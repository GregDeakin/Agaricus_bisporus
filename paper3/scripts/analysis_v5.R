#===============================================================================
#       Load libraries
#===============================================================================

library(limma)
library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
#library(gplot)

#===============================================================================
#       Load one colour Agilent data 
#===============================================================================

# load array location and description
targets<- readTargets("targets_mcb.txt")

# extract descriptions from targets file
colData <- data.frame(Condition=as.factor(sub("[0-9]$","",targets$Condition)),Day=as.factor(sub(".","",targets$Condition)))
levels(colData$Condition) <- c("Treated","Control")
colData$group <- as.factor(paste(colData$Condition,colData$Day,sep="_"))

# function to load one colour data
get_arrays <- function(targets=targets,dpath="../data") {
	RG <- read.maimages(
		targets, 
		path=dpath, 
		columns = list(
			G = "gMedianSignal", 
			Gb = "gBGMedianSignal", 
			R = "gProcessedSignal",
			Rb = "gIsPosAndSignif"
		),
		other.columns=list(
			"gIsPosAndSignif",
			"gIsWellAboveBG",
			"gIsFeatNonUnifOL",
			"gIsBGNonUnifOL",
			"gIsFeatPopnOL",
			"gIsBGPopnOL",
			"IsManualFlag",
			"gIsSaturated"
		),
		annotation = c("Row", "Col","FeatureNum", "ControlType","ProbeName","SystematicName")
	)
	return(RG)
}

# get the array data
RG <- get_arrays(targets)

# backgroung correct data
RG_norm <- backgroundCorrect(RG, method="normexp", offset=50)

# normalize data
RG_norm$G <- normalizeBetweenArrays(RG_norm$G, method="quantile")

# convert expression values to log2 values
RG_norm$G <- log2(RG_norm$G)

#===============================================================================
#	Data filtering
#===============================================================================

# filter function
filter <- function(X) {
	Y <- rowMeans(X)>0.5
	Y <- Y==1
	return(Y)
}

# create filter function for each array "set" 
filter_run <- function(X,n=1) {
	Y=0
	if (n==1){
		Y <- filter(X[,1:4])+ filter(X[,5:8])+ filter(X[,9:12])+ filter(X[,13:16])
	} else {
		Y <- filter(X[,1:4])+ filter(X[,5:8])
	}
	return (Y)
}

#	Y <- filter(X[,1:4])+ filter(X[,5:8])+ filter(X[,9:12])+ filter(X[,13:16])+ filter(X[,17:20])+ filter(X[,21:24])+ filter(X[,25:28])+ filter(X[,29:32])+ filter(X[,33:36])+ filter(X[,37:40])

# apply filter
RG_filt <- RG_norm[filter_run(RG_norm$other$gIsWellAboveBG)>0,]

# remove control genes
RG_filt <- RG_filt[RG_filt$genes$ControlType==0,]

# create new MA list object
RG_ess <- new("MAList", list(targets=RG_filt$targets, genes=RG_filt$genes, source=RG_filt$source,  A=RG_filt$G))

#===============================================================================
#	    Technical Replicates
#===============================================================================

# filter replicates for median values (or set to median if too different)
ff <- function(X) {
	while ( (max(X) - min(X)) > 2 ) {
		if ((max(X)-median(X))>0.5){			
			X[which.max(X)]<-median(X)
		}else{ 
			X[which.min(X)] <- median(X)
		}
	}
	return(median(X))
}

# aggregate by gene name and apply replicates filter (updated to use data tables, it's much faster for this process)
#temp <- aggregate(RG_ess$A,by=list(RG_ess$genes$SystematicName),ff)
#t1 <- data.table(RG_ess$A,genes=RG_ess$genes$SystematicName)
temp <- data.table(RG_ess$A,genes=RG_ess$genes$SystematicName)[,lapply(.SD, ff),by=genes]

# create new MA list
E.avg <- new("MAList", list(targets=RG_filt$targets, genes=temp$genes,source=RG_filt$source,  A=temp[,2:17]))

# remove temp file
rm(temp)

# rename genes
E.avg$genes <- sub("EKV.*_","",E.avg$genes)

#===============================================================================
#	    Annotations
#===============================================================================

# read annotations file
annotations <- read.csv("phd_genes.txt",header=TRUE)

# remove duplicates
annotations <- annotations[!duplicated(annotations$Protein.ID),]

# get the gene names
gene_names <- read.table("gene_names.txt",sep="\t",header=TRUE) 

# join genes to annotations
annotations <- left_join(gene_names,annotations, by=c("JGI_ID"="Protein.ID"))

# annotations <- merge(gene_names,annotations, by.x="JGI_ID",by.y="Protein.ID",all.x=TRUE)

# gene names
genes <- inner_join(data_frame(ENA_ID=E.avg$genes),annotations[,1:2])

#===============================================================================
#	    Statistical Model
#===============================================================================

# add contrast type to colData
contrasts(colData$Condition) <- contr.sum(2)
contrasts(colData$Day) <- contr.sum(2)

# full 2x factor  design model 
design <- model.matrix(~colData$Day*colData$Condition)

# add the design to the gene expression data
fit <- lmFit(E.avg$A, design)

# define the contrast of interest
contrast.matrix <- cbind(
	#Intercept=c(4,0,0,0),
	Day_main_effect=c(0,2,0,0),
	Treatment_main_effect=c(0,0,2,0),
	Day1_Treatment_effect=c(0,0,2,2),
	Day2_Treatment_effect=c(0,0,2,-2),
	Day_Treatment_effect=c(0,2,0,2),
	Day_Control_effect=c(0,2,0,-2),
	Interaction=c(0,0,0,4))

# add contrast to fit
fit2 <- contrasts.fit(fit, contrast.matrix)

# fit the model to the data
fit2 <- eBayes(fit2)

# number of differentially expressed genes
res <- data.table(decideTests(fit2))
res$Day1_Treatment_AND_Day2_Treatment <- (res[,3]&res[,4])*res[,3]#((res[,3]&res[,4] ) - (res[,3]&res[,4] )*res[,3])/2 + ((res[,3]&res[,4] ) + (res[,3]&res[,4] )*res[,3])/2
res$Day1_Treatment_NOT_Day2_Treatment <- (res[,3]&!res[,4])*res[,3] #((res[,3]&!res[,4]) - (res[,3]&!res[,4])*res[,3])/2 + ((res[,3]&!res[,4]) + (res[,3]&!res[,4])*res[,3])/2
res$Day2_Treatment_NOT_Day1_Treatment <- (res[,4]&!res[,3])*res[,4] #(((res[,4]&!res[,3]) - (res[,4]&!res[,3]))*res[,4])/2 + ((res[,4]&!res[,3]) + (res[,4]&!res[,3])*res[,4])/2
res$Treatment_NOT_Day <-(res[,2]&!res[,1])*res[,2]

df <- data.frame(
	"increased_expression"=apply(res*(res==1),2,sum),
	"decreased_expression"=apply(res*(res==-1)*-1,2,sum),
	"no_change"=mapply(sum,res*(res==1)*-1,res*(res==-1),res*(res==0)+1))

reslist <- list(
	Day_main_effect=topTable(fit2, adjust="BH", coef=1, genelist=E.avg$genes, number=length(E.avg$genes)),
	Treatment_main_effect= topTable(fit2, adjust="BH", coef=2, genelist=E.avg$genes, number=length(E.avg$genes)),
	Day1_Treatment_effect=topTable(fit2, adjust="BH", coef=3, genelist=E.avg$genes, number=length(E.avg$genes)),
	Day2_Treatment_effect= topTable(fit2, adjust="BH", coef=4, genelist=E.avg$genes, number=length(E.avg$genes)),
	Treatment_effect=topTable(fit2, adjust="BH", coef=5, genelist=E.avg$genes, number=length(E.avg$genes)),
	Control_effect=topTable(fit2, adjust="BH", coef=6, genelist=E.avg$genes, number=length(E.avg$genes)),
	Interaction=topTable(fit2, adjust="BH", coef=7, genelist=E.avg$genes, number=length(E.avg$genes))
)


#===============================================================================
#	    Data analysis
#===============================================================================

### Heirachical clustering
gene_exprs <- E.avg$A
colnames(gene_exprs) <- paste(sub("C","Control_",sub("A","Treated_",targets$Condition)),seq(1,4),sep="_")
d <- dist(gene_exprs,method="euclidean")	  
h <- hclust(d)
gene_exprs$ID <- E.avg$genes
gene_exprs$order <- h$order
gene_exprs$cuts <- cutree(h, k=100)
#gene_exprs$genes <- row.names(gene_exprs)
#gene_exprs <- left_join(gene_exprs,annotations,by=c("genes"="ENA_ID"))
#gene_exprs <- gene_exprs[order(gene_exprs$order),]
#########

# put all the results (the interesting columns) in a single table
results <- Reduce(function(...) 
	inner_join(...,by="ID"),
		      lapply(reslist,'[',c(1,2,6)
			    ))
# give the columns some names
colnames(results)[-1] <-  rbind(paste(names(reslist),"Log2FC",sep="_"),paste(names(reslist),"adjP",sep="_"))
# add expression values (and hc order)
results <- inner_join(results,gene_exprs)
# add annotations
results <- left_join(results,annotations,by=c("ID"="ENA_ID"))
# write out the table
write.table(results,"new_results.txt",sep="\t",row.names=F,quote=F,na="")


#===============================================================================
#	    Plots
#===============================================================================

##### pca plot using metabarcoding plotOrd function ####
#colData <-sample_info[,2:3]
#colnames(colData) <- c("Condition","Batch")
#colData[colData$Condition=="compost_control",1] <- "C"
#colData[colData$Condition=="compost_mvx",1] <- "A"
#colData[colData$Batch=="A",2] <- "V4"
#colData[colData$Batch=="B",2] <- "V5"

rownames(colData) <- colnames(E.avg$A)
	     
mypca <- prcomp(t(E.avg$A))
mypca$percentVar <- mypca$sdev^2/sum(mypca$sdev^2)
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

ggsave("PCA_NEW_COLOUR.pdf",plotOrd(df,colData,design="Condition",shape="Day",alpha=0.75))

g_pca <- plotOrd(df,colData,design="Sample",alpha=0.75,textSize=16,legend="bottom",pointSize=1.5)#ylims=c(-8,16))

# g <- plotOrd(df,colData,dimx=1,dimy=2,design="Condition",xlabel="PC1",ylabel="PC2") 
# g_pca <- g + theme_classic_thin(16) %+replace% theme(legend.position = "bottom",axis.text = element_text(colour="grey20"))

#### ma plot using function below ####
plot_ma(mvx_effect,xlim=c(-4,4))
plot_ma(day1_effect,xlim=c(-4,4))
plot_ma(day2_effect,xlim=c(-4,4))
plot_ma(day2_effect,legend=T)

#### virus plot ####
# get list of viruses	     
viruses <- data.table(contig=E.avg$genes[grep("^C\\d.*",E.avg$genes,perl=T)],E.avg$A[grep("^C\\d.*",E.avg$genes,perl=T),])
# get new virus names
vnames <- read.table("viruses.txt",header=T,sep="\t")
# add sample names  
colnames(viruses)[-1] <- paste(sub("C","Control_",sub("A","Treated_",targets$Condition)),seq(1,4),sep="_")
#rownames(viruses) <- vnames$virus
# subset viruses (get rid of rev)
viruses <- left_join(vnames,viruses)
#colnames(viruses) <- paste(targets_mcb$Condition,seq(1,4),sep="_")
 

# melt the data
test2 <- melt(viruses)
# set virus factor order
test2$virus <- factor(test2$virus, levels = as.factor(viruses[hclust(dist((viruses[,c(-1,-2)])))$order,2]))	
# remove contig column
test2 <- test2[-1]
#colnames(test2)[3] <- "Scale"
# mean centre the data
test2$Scale <- scale(test2[,3],scale=F)

pdf("virus_plot_3.pdf")	   
g <- ggplot(test2,aes(x=Var2,y=Var1,fill=Scale))
g<- g+ geom_tile(colour = "black")
g <- g + scale_fill_gradient2(mid="orange", low = "red",high = "yellow", na.value = "black")
g <- g + labs(x=NULL,y=NULL)
g + labs(title=" |------Day1------|   |------Day2------|")+ theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size =16),plot.title = element_text(size=16))
dev.off()

# V2
library(grid)
library(gridExtra)
g2 <- g + labs(title=" |------Day1------|   |------Day2------|")
g2 <- g2 + theme(axis.text=element_text(colour = "grey20"),axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size =16),plot.title = element_text(size=16,colour="grey20"))
g2 <- g2 + annotation_custom(textGrob("Treated",gp = gpar(fontsize = 15,col="grey20")),xmin=2.5,xmax=2.5,ymin=0,ymax=0)
g2 <- g2 + annotation_custom(textGrob("Control",gp = gpar(fontsize = 15,col="grey20")),xmin=6.5,xmax=6.5,ymin=0,ymax=0)
g2 <- g2 + annotation_custom(textGrob("Treated",gp = gpar(fontsize = 15,col="grey20")),xmin=10.5,xmax=10.5,ymin=0,ymax=0)
g2 <- g2 + annotation_custom(textGrob("Control",gp = gpar(fontsize = 15,col="grey20")),xmin=14.5,xmax=14.5,ymin=0,ymax=0)
g2 <- ggplot_gtable(ggplot_build(g2))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.arrange(g2)
dev.off()

# V3
# preprepared virus table
viruses <- read.table("virus2.txt",header=T,sep="\t",row.names=1)

# get probalities
padj <- read.table("virus_padj.txt",header=T,sep="\t")
padj <- melt(padj)

# get expression means
v3 <- sapply(seq(1,4),function(x) rowMeans(viruses[,(4*x-3):(4*x)]))
colnames(v3) <- c("T1","C1","T2","C2")

# melt expression values
test2 <- melt(as.matrix(v3))

# join p values to mean expression values
test2 <- left_join(test2,padj,by=c("Var1"="ID","Var2"="variable"))
test2$value.y[is.na(test2$value.y)] <- 1
test2$label <- sapply(test2$value.y,function(x) if(as.numeric(x)<=0.05) {return("\u2020")}else{return("")})
test2$Var1 <- factor(test2$Var1, levels = as.factor(row.names(viruses)[ hclust(dist((viruses)))$order]))
# scale expression values
test2$Scale <- as.vector(scale(test2[,3],scale=F))
		      
# plot sprintf("%0.2f", round(Scale, digits = 1))
g <- ggplot(test2,aes(x=Var2,y=Var1,fill=Scale))
g<- g+ geom_tile(colour = "black") + geom_text(label=test2$label,size=2.8,hjust = -3,vjust=-0.1)
g <- g + scale_fill_gradient2(mid="orange", low = "red",high = "yellow", na.value = "black")
g <- g + labs(x=NULL,y=NULL)
g <- g + theme(text = element_text(size =16),legend.position = "bottom",legend.text = element_text(size = 10))
ggsave("virus_plot_3.pdf",g,device=cairo_pdf,width=4)
g_virus <- g

#### anti-viral heatplot ####
# get genes
av <- read.table("anti_genes.txt",header=T,sep="\t",)
av$JGI_ID <- as.character(av$JGI_ID) # could do this in the read.table

# get probalities
padj <- av
padj <- inner_join(padj,results[,c(36,7,9)])
colnames(padj)[3:4] <- c("T1","T2")
padj <- melt(padj[,c(2,3,4)])

# get expression values
av <- inner_join(av,annotations[,c(1,4)])
av <- inner_join(av[,2:3],gene_exprs[,1:17],by=c("ENA_ID"="ID"))
rownames(av) <- av$Description
av <- av[,c(-1,-2)]
av <- sapply(seq(1,4),function(x) rowMeans(av[,(4*x-3):(4*x)] ))
# melt and scale expression values
test2 <- t(apply(av,1,scale,scale=F))
colnames(test2) <-c("T1","C1","T2","C2")
test2 <- melt(as.matrix(test2))
#test2$Scale <- scale(test2[,3],scale=F)
test2$Scale <- test2$value

# join p values to mean expression values
test2 <- left_join(test2,padj,by=c("Var1"="Description","Var2"="variable"))
test2$value.y[is.na(test2$value.y)] <- 1
test2$label <- sapply(test2$value.y,function(x) if(as.numeric(x)<=0.05){return("\u2020")}else{return("")})
#test2$Var2 <-  as.factor(test2$Var2)
		      
# plot  
g <- ggplot(test2,aes(x=Var2,y=Var1,fill=Scale))
g<- g+ geom_tile(colour = "black") + geom_text(label=test2$label,size=2.8,hjust = -3,vjust=0.1)
g <- g + scale_fill_gradient2(mid="orange", low = "red",high = "yellow", na.value = "black")
g <- g + labs(x=NULL,y=NULL)
g <- g + theme(text = element_text(size =16),legend.position = "bottom")
#ggsave("av_plot_1.pdf",g,device=cairo_pdf,width=4)
g_anti <- g

# V2
#library(grid)
#library(gridExtra)
g2 <- g + labs(title=" |------Day1------|   |------Day2------|")
g2 <- g2 + theme(axis.text=element_text(colour = "grey20"),axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size =16),plot.title = element_text(size=16,colour="grey20"))
g2 <- g2 + annotation_custom(textGrob("Treated",gp = gpar(fontsize = 15,col="grey20")),xmin=2.5,xmax=2.5,ymin=0,ymax=0)
g2 <- g2 + annotation_custom(textGrob("Control",gp = gpar(fontsize = 15,col="grey20")),xmin=6.5,xmax=6.5,ymin=0,ymax=0)
g2 <- g2 + annotation_custom(textGrob("Treated",gp = gpar(fontsize = 15,col="grey20")),xmin=10.5,xmax=10.5,ymin=0,ymax=0)
g2 <- g2 + annotation_custom(textGrob("Control",gp = gpar(fontsize = 15,col="grey20")),xmin=14.5,xmax=14.5,ymin=0,ymax=0)
g2 <- ggplot_gtable(ggplot_build(g2))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.arrange(g2)
dev.off()

#### combine PCA, Virus (v3) and anti_viral plots ####
title.a <- textGrob(label = "a",x = unit(0, "lines"),y = unit(0, "lines"),hjust = -0.5, vjust = 0,gp = gpar(fontsize = 20,face="bold"))
title.b <- textGrob(label = "b",x = unit(0, "lines"),y = unit(0, "lines"),hjust = -0.5, vjust = 0,gp = gpar(fontsize = 20,face="bold"))
title.c <- textGrob(label = "c",x = unit(0, "lines"),y = unit(0, "lines"),hjust = -0.5, vjust = 0,gp = gpar(fontsize = 20,face="bold"))
g2 <- arrangeGrob(g_pca+theme_classic_thin(base_size=14), top = title.a)
g3 <- arrangeGrob(g_virus, top = title.b)
g4 <- arrangeGrob(g_anti, top = title.c)
g <- grid.arrange(g2,g3,g4,layout_matrix=rbind(c(1,1),c(2,3),c(2,3)))
		      
#g <- grid.arrange(arrangeGrob(g2),arrangeGrob(g3,g4,ncol=2),heights=unit.c(unit(1, "npc") - unit(22, "lines"),unit(22, "lines")),nrow=2)
ggsave("Figure_4_v3.pdf",g,device=cairo_pdf,width=8,height=9)

#### Enzymes ####
enzymes <- read.table("enzymes.txt",header=T,sep="\t",)
enzymes$JGI_ID <- as.character(enzymes$Protein.ID)

# get probalities
padj <- enzymes
padj <- inner_join(padj[,2:4],results[,c(36,7,9)])
colnames(padj)[4:5] <- c("T1","T2")
padj <- melt(padj[,c(1,2,4,5)])

# get expression values
enzymes <- inner_join(enzymes[,2:4],annotations[,c(1,4)])
enzymes <- inner_join(enzymes,gene_exprs[,1:17],by=c("ENA_ID"="ID"))
enzymes2 <- sapply(seq(1,4),function(x) rowMeans(enzymes[,(4*x+1):(4*x+4)] ))
# melt and scale expression values
test2 <- t(apply(enzymes2,1,scale,scale=F))
colnames(test2) <-c("T1","C1","T2","C2")
test2 <- cbind(enzymes[,1:3],test2)
test2 <- melt(test2)

test2$Scale <- test2$value
test2$group_var <- paste(test2$name,test2$JGI_ID,sep="_")

# join p values to mean expression values
test2 <- left_join(test2,padj)
test2$label <- sapply(test2$T1,function(x) if(as.numeric(x)<=0.05){return("\u2020")}else{return("")})
#test2$Var2 <-  as.factor(test2$Var2)
		      
# plot  
g <- ggplot(test2,aes(x=variable,y=name,fill=Scale))
g<- g+ geom_tile(colour = "black") + geom_text(label=test2$label,size=2.8,hjust = -3,vjust=-0.1)
g <- g + scale_fill_gradient2(mid="orange", low = "red",high = "yellow", na.value = "black")
g <- g + labs(x=NULL,y=NULL)
g <- g + theme(text = element_text(size =16),legend.position = "bottom")
ggsave("enzyme_plot_1.pdf",g,device=cairo_pdf,width=4)

############################################
############################################
		      
		      
# heatmap of log2foldchange for each contrast	

dfe <- Reduce(function(...) 
	inner_join(...,by="ID"),
		      lapply(reslist,'[',c(1,2)
			    ))
dfe <- inner_join(dfe,gene_exprs[,c(17,19)])
test2 <- melt(dfe)
#test2$order <- as.factor(dfe$ID[dfe$cuts]))
=2^abs(value)/(abs(value)/value)
pdf("test_plot_2.pdf")	   
g <- ggplot(test2,aes(x=variable,y=cuts,fill=value))
#g<- g+ geom_tile(colour = "white")
g <- g + geom_raster()
g <- g + scale_fill_gradient2(mid="yellow", low = "blue",high = "red", na.value = "black")
g <- g + scale_fill_gradient2(mid="orange", low = "red",high = "yellow", na.value = "black")
g <- g + labs(x=NULL,y=NULL)
g + theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size =16),axis.text.y=element_blank(),axis.ticks.y=element_blank())
dev.off()
# this is rubbish

# OLD anti-viral mechanisms
anti<-read.table("anti.txt",sep="\t",header=T)
anti[,1]<-as.character(anti[,1])
anti_merged <- left_join(anti[,1:2],annotations[,c(1,4)])

anti<-left_join(anti,A[,c(1,2,6)])
anti<-left_join(anti,C1[,c(1,2,6)],by="JGI_ID")
anti<-left_join(anti,C2[,c(1,2,6)],by="JGI_ID")
anti<-left_join(anti,A1[,c(1,2,6)],by="JGI_ID")
anti<-left_join(anti,A2[,c(1,2,6)],by="JGI_ID")

colnames(anti) <- c("JGI_ID","Description","MVX","p_A","C_Day1","p_C1","C_Day2","p_C2","T_Day1","p_A1","T_Day2","p_A2")
anti$Description<-factor(anti$Description)

test <- melt(anti[,c(2,3,9,11,5,7)])
test2<- melt(anti[,c(2,4,10,12,6,8)])
test$value<-2^(abs(test$value ))*abs(test$value )/test$value 
colnames(test)[2:3]<-c("Sample","Fold_change")
colnames(test2)[2]<-"Sample"

frames<-test2[test2$value<=0.05,1:2]
frames$Sample = as.integer(frames$Sample)
frames$Description = as.integer(frames$Description)

	     
g <- ggplot(test,aes(x=Sample,y=Description,fill=Fold_change))
g<- g+ geom_tile(colour = "white")
#g <- g + geom_raster()
g <- g + scale_fill_gradient2(mid="white", low = "steelblue", high = "red", na.value = "black")
g<-g+geom_rect(data=frames, size=1, fill=NA, colour="black",aes(xmin=Sample - 0.5, xmax=Sample + 0.5, ymin=Description - 0.5, ymax=Description + 0.5))


	     
g <- ggplot(test3,aes(x=Sample,y=Description,fill=value.x,z=hl))
g<- g+geom_rect(size=1,fill=NA,colour="yellow",aes(xmin=z*(Sample-0.5),xmax=z*(Sample+0.5),ymin=z*(Sample-0.5),ymax=z*(Sample+0.5)))



	   

geom_raster(aes(x=Var1, y=Var2, fill=value)) +
     scale_fill_gradient2(low="blue", high="red", na.value="black", name="") +
     geom_rect(data=frames, size=1, fill=NA, colour="black",
       aes(xmin=Var1 - 0.5, xmax=Var1 + 0.5, ymin=Var2 - 0.5, ymax=Var2 + 0.5)) 


#===============================================================================
#	   Functions
#===============================================================================

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_ma <- function
(
	fitObj,
	xlim=c(-6,6),
	textsize=16,
	legend=F,
	crush=T
)

{
	d <- fitObj[,c(2,3,6)]
	colnames(d) <- c("log2FoldChange","baseMean","padj")
	d$group<-1
	d[d$padj<=0.05,4]<-2 
	d[abs(d$log2FoldChange)>1,4]<-3 
	d[(d$padj<=0.05)&(abs(d$log2FoldChange)>1),4]<-4
	d$group<-as.factor(d$group)
	d$shape<-16

	if(crush){
		d[d$log2FoldChange<xlim[1],5]<-25
		d[d$log2FoldChange<xlim[1],1]<-xlim[1]
		d[d$log2FoldChange>xlim[2],5]<-24
		d[d$log2FoldChange>xlim[2],1]<-xlim[2]
	}

	g <- ggplot(data=d,aes(x=log2FoldChange,y=baseMean,colour=group,shape=shape))
	g <- g + theme_bw()
	g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	if(!legend) {
		g <- g+ theme(legend.position="none")
	}
	g <- g + theme(axis.line.x = element_line(size=0.3,colour = "black"),
	     axis.line.y = element_line(size=0.3,colour = "black"),
	     axis.text = element_text(colour = "black"),
	     text=element_text(size=16)
	)
	g <- g + scale_shape_identity() 
	g <- g + geom_point(size=3)
	g <- g + scale_colour_manual(values=cbbPalette)
	g <- g + xlab(expression("Log"[2]*" Fold Change"))
	g <- g + ylab(expression("Log"[2]*" Mean Expression"))
	g <- g + xlim(xlim)
	g <- g + expand_limits(x = xlim[1], y = 5)
	g <- g + coord_flip()
	return(g)
}


#########################################
#########################################
		      
kogclass<-fread("kog_class",sep="\t")
	     
# apply(kogclass,1,function(str) nrow(day1[kogClass %like% str & adj.P.Val <= 0.05 & logFC > 0 ]))

kog_nums <- sapply(seq_along(kogclass),function(i) {
		data.frame(all=	nrow(day1[kogClass %like% kogclass[i]]),
		   day1_up=nrow(day1[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC > 0 ]),
		   day1_down=nrow(day1[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC < 0 ]),
		   day2_up=nrow(day2[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC > 0 ]),
		   day2_down=nrow(day2[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC < 0 ])
		)
	}
)

	     
	     
day_1_2 <- data.table(inner_join(day1,day2[,c(1,2,4)],by="ID"))

colnames(kog_nums) <- t(kogclass)
x<-t(apply(kog_nums,1,function(i) prop.table(as.numeric(i))))
## ignnore
	   ))}
###
	     
colnames(x) <- t(kogclass) 
df <- melt(x, id = row.names)
ggplot(df, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity",colour="white")+theme(legend.position="none")	   


colnames(kog_nums_2) <- t(kogclass)
df <- data.table(kog_nums_2,keep.rownames = TRUE)
df<-melt(df,id="rn")
df$value <- as.numeric(df$value)
   
	     
ggplot(df, aes(x = rn, y = value, fill = variable)) + geom_bar(stat = "identity",colour="white")+theme(legend.position="none")

		      

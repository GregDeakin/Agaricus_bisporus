#===============================================================================
#       Load libraries
#===============================================================================

library(limma)
library(virtualArray)
library(dplyr)
library(affyPLM)
library(reshape2)
library(ggplot2)
library(gplot)

#===============================================================================
#       Load one colour Agilent data 
#===============================================================================

targets_v5<- readTargets("targets_mcb.txt")
targets_v4<- readTargets("targets_phd.txt")

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

RG_v5 <- get_arrays(targets_v5)
RG_v4 <- get_arrays(targets_v4)

RG_norm_v5 <- backgroundCorrect(RG_v5, method="normexp", offset=50)
RG_norm_v4 <- backgroundCorrect(RG_v4, method="normexp", offset=50)

# need to work on whether this should be done now or only later...
RG_norm_v4$G <- normalizeBetweenArrays(RG_norm_v4$G, method="quantile")
RG_norm_v5$G <- normalizeBetweenArrays(RG_norm_v5$G, method="quantile")

RG_norm_v5$G <- log2(RG_norm_v5$G)
RG_norm_v4$G <- log2(RG_norm_v4$G)

#RG_ess <- new("MAList", list(targets=RG_norm$targets, genes=RG_norm$genes, source=RG_norm$source,  A=RG_norm$G))

#===============================================================================
#	Data filtering
#===============================================================================

filter <- function(X) {
	Y <- rowMeans(X)>0.5
	Y <- Y==1
	return(Y)
}

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

RG_filt_v5 <- RG_norm_v5[filter_run(RG_norm_v5$other$gIsWellAboveBG)>0,]
RG_filt_v4 <- RG_norm_v4[filter_run(RG_norm_v4$other$gIsWellAboveBG,2)>0,]

RG_filt_v5 <- RG_filt_v5[RG_filt_v5$genes$ControlType==0,]
RG_filt_v4 <- RG_filt_v4[RG_filt_v4$genes$ControlType==0,]

RG_ess_v5 <- new("MAList", list(targets=RG_filt_v5$targets, genes=RG_filt_v5$genes, source=RG_filt_v5$source,  A=RG_filt_v5$G))
RG_ess_v4 <- new("MAList", list(targets=RG_filt_v4$targets, genes=RG_filt_v4$genes, source=RG_filt_v4$source,  A=RG_filt_v4$G))

#===============================================================================
#	    Technical Replicates
#===============================================================================

# below for v4 microarray only
RG_ess_v4$genes$SystematicName[grep("Agabi_varbisH97_2$",RG_ess_v4$genes$SystematicName)] <- RG_ess_v4$genes$ProbeName[grep("Agabi_varbisH97_2$",RG_ess_v4$genes$SystematicName)]


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

temp <- aggregate(RG_ess_v5$A,by=list(RG_ess_v5$genes$SystematicName),ff)
E.avg.v5 <- new("MAList", list(targets=RG_filt_v5$targets, genes=temp$Group.1,source=RG_filt_v5$source,  A=temp[,2:17]))
# for V5 only
E.avg.v5$genes <- sub(".*_","",E.avg.v5$genes)
			
temp <- aggregate(RG_ess_v4$A,by=list(RG_ess_v4$genes$SystematicName),ff)
E.avg.v4 <- new("MAList", list(targets=RG_filt_v4$targets, genes=temp$Group.1,source=RG_filt_v4$source,  A=temp[,2:9]))
rm(temp)

#===============================================================================
#	    Merge arrays (needs a bit of work to add viruses and filtered from one only)
#===============================================================================

annotations <- read.csv("phd_genes.txt",header=TRUE)
annotations <- annotations[!duplicated(annotations$Protein.ID),]
gene_names <- read.table("gene_names.txt",sep="\t",header=TRUE) 
annotations <- merge(gene_names,annotations, by.x="JGI_ID",by.y="Protein.ID",all.x=TRUE)

genes_v4 <- inner_join(data_frame(Name=E.avg.v4$genes),annotations[,c(1,6)])
genes_v5 <- inner_join(data_frame(ENA_ID=E.avg.v5$genes),annotations[,1:2])

genes_to_keep <- inner_join(genes_v4,genes_v5)

test<- cbind(E.avg.v4$genes,E.avg.v4$A)
colnames(test)[1] <- "Name"
test2 <- inner_join(test,genes_to_keep)
test.v4 <- test2

test<- cbind(E.avg.v5$genes,E.avg.v5$A)
colnames(test)[1] <- "ENA_ID"
test2 <- inner_join(test,genes_to_keep)
test.v5 <- test2

test <- inner_join(test.v4,test.v5)

E.avg <- new("MAList", list(targets=rbind(RG_ess_v4$targets,RG_ess_v5$targets), genes=test$JGI_ID,gi=test$ENA_ID,source=RG_ess_v4$source,  A=test[,c(2:9,12:27)]))

#===============================================================================
#      Normalisation
#===============================================================================

merged <- new("ExpressionSet", exprs = as.matrix(E.avg$A))
merged_norm <- normalize.ExpressionSet.quantiles(merged)
sample_info <- E.avg$targets
colnames(sample_info)[1:2] <-c("Array.name", "Sample.name")
sample_info$Batch <- c(rep("A",8),rep("B",16))			
sampleNames(merged_norm) <- sample_info[, 1]

merged_combat <- merged_norm
exprs(merged_combat) <- virtualArrayComBat(expression_xls = exprs(merged_combat),sample_info_file=sample_info)

pData(merged_combat) <- as.data.frame(sample_info)
merged_combat_norm <- normalize.ExpressionSet.quantiles(merged_combat)
E.avg$M <- exprs(merged_combat_norm)

#===============================================================================
#	    Statistical Analysis
#===============================================================================

# combined
f <- factor(colData$Condition, levels = unique(colData$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(E.avg$M, design)
contrast.matrix <- makeContrasts(
	"A-C",
	"C1-C",
	"C2-C",
	"A1-C",
	"A2-C",
	"A1-C1",
	"A2-C2",
	"A-C1",
	"A-C2",
	"(C1+C2+A)/3-C",
	levels=design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# v4 only
f <- factor(colData$Condition[1:8], levels = unique(colData$Condition[1:8]))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(E.avg.v4$A, design)
contrast.matrix <- makeContrasts(
	"A-C",
	levels=design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2.v4 <- fit2

# v5 only
f <- factor(colData$Condition[9:length(colData$Condition)], levels = unique(colData$Condition[9:length(colData$Condition)]))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(E.avg.v5$A, design)
contrast.matrix <- makeContrasts(
	"A1-C1",
	"A2-C2",
	levels=design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2.v5 <- fit2

#===============================================================================
#	    Data analysis
#===============================================================================

mvx <- compost <- topTable(fit2.v4, adjust="BH", coef="A-C", genelist=E.avg.v4$genes, number=length(E.avg.v4$genes))
day1 <- topTable(fit2.v5, adjust="BH", coef="A1-C1", genelist=E.avg.v5$genes, number=length(E.avg.v5$genes))
day2 <- topTable(fit2.v5, adjust="BH", coef="A2-C2", genelist=E.avg.v5$genes, number=length(E.avg.v5$genes))

mvx <- inner_join(mvx[,c(1:3,6)],annotations,by=c("ID"="Name"))
day1 <- inner_join(day1[,c(1:3,6)],annotations,by=c("ID"="ENA_ID"))
day2 <- inner_join(day2[,c(1:3,6)],annotations,by=c("ID"="ENA_ID"))

mvx_effect <- compost <- topTable(fit2, adjust="BH", coef="A-C", genelist=E.avg$genes, number=length(E.avg$genes))
names(mvx_effect)[1] <- "JGI_ID"

day1_effect <- topTable(fit2, adjust="BH", coef="A1-C1", genelist=E.avg$genes, number=length(E.avg$genes))
day2_effect <- topTable(fit2, adjust="BH", coef="A2-C2", genelist=E.avg$genes, number=length(E.avg$genes))

A <- mvx_effect
C1 <- topTable(fit2, adjust="BH", coef="C1-C", genelist=E.avg$genes, number=length(E.avg$genes))
C2 <- topTable(fit2, adjust="BH", coef="C2-C", genelist=E.avg$genes, number=length(E.avg$genes))
A1 <- topTable(fit2, adjust="BH", coef="A1-C", genelist=E.avg$genes, number=length(E.avg$genes))
A2 <- topTable(fit2, adjust="BH", coef="A2-C", genelist=E.avg$genes, number=length(E.avg$genes))

colnames(A)[1] <- "JGI_ID"
colnames(C1)[1] <- "JGI_ID"
colnames(C2)[1] <- "JGI_ID"
colnames(A1)[1] <- "JGI_ID"
colnames(A2)[1] <- "JGI_ID"

con_mvx <- inner_join(A[,c(1:3,6)],inner_join(C1[,c(1:3,6)],inner_join(C2[,c(1:3,6)],inner_join(A1[,c(1:3,6)],A2[,c(1:3,6)],by="JGI_ID"),by="JGI_ID"),by="JGI_ID"),by="JGI_ID")
con_mvx <- con_mvx[,c(1,2,4,5,7,8,10,11,13,14,16,3)]
colnames(con_mvx) <- c("JGI_ID","FC_A","P_A","FC_C1","P_C1","FC_C2","P_C2","FC_A1","P_A1","FC_A2","P_A2","BaseMean")

mvx_effect_all <- topTable(fit2, adjust="BH", coef="(C1+C2+A)/3-C", genelist=E.avg$genes, number=length(E.avg$genes))
dim(mvx_effect_all[mvx_effect_all$adj.P.Val<=0.05,]) 
colnames(mvx_effect_all)[1] <- "JGI_ID"

#mvx_effect <- compost <- topTable(fit2, adjust="BH", coef="(compost_mvx+C1+C2+A1+A2)-compost_control", genelist=E.avg$genes, number=length(E.avg$genes))

mvx_effect_all <- inner_join(mvx_effect_all,annotations)


kogclass<-fread("kog_class",sep="\t")

# apply(kogclass,1,function(str) nrow(day1[kogClass %like% str & adj.P.Val <= 0.05 & logFC > 0 ]))


kog_nums <- sapply(seq(1,nrow(kogclass)),function(i) {
		data.frame(all=	nrow(day1[kogClass %like% kogclass[i]]),
		   day1_up=nrow(day1[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC > 0 ]),
		   day1_down=nrow(day1[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC < 0 ]),
		   day2_up=nrow(day2[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC > 0 ]),
		   day2_down=nrow(day2[kogClass %like% kogclass[i] & adj.P.Val <= 0.05 & logFC < 0 ])
		)
	}
)

day_1_2 <- data.table(inner_join(day1,day2[,c(1,2,4)],by="ID"))

kog_nums_2 <- sapply(seq(1,nrow(kogclass)),function(i) {
		data.frame(all=	nrow(day1[kogClass %like% kogclass[i]]),
		   early_up=nrow(day_1_2[kogClass %like% kogclass[i] & adj.P.Val.x <= 0.05 & logFC.x > 0 & adj.P.Val.y > 0.05 ]),
		   early_down=nrow(day_1_2[kogClass %like% kogclass[i] & adj.P.Val.x <= 0.05 & logFC.x < 0 & adj.P.Val.y > 0.05]),
		   late_up=nrow(day_1_2[kogClass %like% kogclass[i] & adj.P.Val.y <= 0.05 & logFC.y > 0 & adj.P.Val.x > 0.05]),
		   late_down=nrow(day_1_2[kogClass %like% kogclass[i] & adj.P.Val.y <= 0.05 & logFC.y < 0& adj.P.Val.x > 0.05 ]),
		   long_up=nrow(day_1_2[kogClass %like% kogclass[i] & adj.P.Val.y <= 0.05 & logFC.y > 0 & adj.P.Val.x <= 0.05]),
		   long_down=nrow(day_1_2[kogClass %like% kogclass[i] & adj.P.Val.y <= 0.05 & logFC.y < 0 & adj.P.Val.x <= 0.05])			   
		)
	}
)


colnames(kog_nums) <- t(kogclass)
x<-t(apply(kog_nums,1,function(i) prop.table(as.numeric(i))))
colnames(x) <- t(kogclass) 
df <- melt(x, id = row.names)
ggplot(df, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity",colour="white")+theme(legend.position="none")	   

colnames(kog_nums_2) <- t(kogclass)
df <- data.table(kog_nums_2,keep.rownames = TRUE)
df<-melt(df,id="rn")
df$value <- as.numeric(df$value)
   
ggplot(df, aes(x = rn, y = value, fill = variable)) + geom_bar(stat = "identity",colour="white")+theme(legend.position="none")
	   
	   
#===============================================================================
#	    Plots
#===============================================================================

# pca plot using metabarcoding plotOrd function
colData <-sample_info[,2:3]
colnames(colData) <- c("Condition","Batch")
colData[colData$Condition=="compost_control",1] <- "C"
colData[colData$Condition=="compost_mvx",1] <- "A"
colData[colData$Batch=="A",2] <- "V4"
colData[colData$Batch=="B",2] <- "V5"

mypca <- prcomp(t(E.avg$M))
mypca$percentVar <- mypca$sdev^2/sum(mypca$sdev^2)
df <- t(data.frame(t(mypca$x)*mypca$percentVar))
plotOrd(df,colData,dimx=1,dimy=2,design="Condition",xlabel="PC1",ylabel="PC2")

# ma plot using function below
plot_ma(mvx_effect,xlim=c(-4,4))
plot_ma(day1_effect,xlim=c(-4,4))
plot_ma(day2_effect,xlim=c(-4,4))
plot_ma(day2_effect,legend=T)


anti<-read.table("anti.txt",sep="\t",header=T)
anti[,1]<-as.factor(anti[,1])

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
g

g <- ggplot(test3,aes(x=Sample,y=Description,fill=value.x,z=hl))
g<- g+geom_rect(size=1,fill=NA,colour="yellow",aes(xmin=z*(Sample-0.5),xmax=z*(Sample+0.5),ymin=z*(Sample-0.5),ymax=z*(Sample+0.5)))

## virus plot
viruses <- E.avg.v5$A[grep("^C\\d.*",E.avg.v5$genes,perl=T),]
vnames <- read.table("viruses.txt",header=T,sep="\t")
rownames(viruses) <- vnames$virus
#colnames(viruses) <- paste(targets_mcb$Condition,seq(1,4),sep="_")
   
colnames(viruses) <- paste(sub("C","Control_",sub("A","Treated_",targets_mcb$Condition)),seq(1,4),sep="_")
#test2 <- melt(as.matrix(viruses))

test3 <- apply(viruses,2, scale,scale=F)
test3 <- t(apply(viruses,1, scale,scale=F))
colnames(test3) <- colnames(viruses)	   
test2 <- melt(as.matrix(test3))
test2$Var1 <- factor(test2$Var1, levels = as.factor(row.names(viruses)[ hclust(dist((viruses)))$order]))	   
colnames(test2)[3] <- "Scale"

pdf("virus_plot_2.pdf")	   
g <- ggplot(test2,aes(x=Var2,y=Var1,fill=Scale))
g<- g+ geom_tile(colour = "black")
g <- g + scale_fill_gradient2(mid="orange", low = "red",high = "yellow", na.value = "black")
g <- g + labs(x=NULL,y=NULL)
g + labs(title="  |------Day1------|   |------Day2------|")+ theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size =16),plot.title = element_text(size=16))
dev.off()
   
	   

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

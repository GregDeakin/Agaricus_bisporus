#===============================================================================
#       Load libraries
#===============================================================================
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(vegan)
library(lmPerm)
library(ape)
library(gridExtra)
library(grid)
library(fpc)
library("factoextra")
library(clustertend)
library(cluster)
library(mixtools)
library(ggdendro)


load_all("~/pipelines/metabarcoding/scripts/myfunctions")

#===============================================================================
#       Load data
#===============================================================================

countData <- fread("countData")
colData <- fread("colData")

DF <- countData[,4:21]

list_DT <- list(both=DF,exp1=countData[Experiment==1,4:21],exp2=countData[Experiment==2,4:21])


#===============================================================================
#       Cluster analysis
#===============================================================================

# Hopskin statistic
set.seed(sum(utf8ToInt("Kerry Burton")))
lapply(list_DT,get_clust_tendency, 17, graph = F)

$both$hopkins_stat
[1] 0.2385218
$exp1$hopkins_stat
[1] 0.2975166
$exp2$hopkins_stat
[1] 0.3145393

# Ward Hierarchical Clustering
d <- lapply(list_DT,function(DF) {dist(t(DF), method = "euclidean")}) # distance matrix
fit <- lapply(d,hclust,method="ward.D2")
#lapply(fit,plot)  # display dendogram
g_hi <- lapply(fit,function(fit) {ggdendrogram(fit,rotate=F)+ theme_blank(base_size=16) %+replace%
 theme(axis.text.x = element_text(angle = 45, vjust = 1.1,hjust = 1),
 axis.ticks=element_blank(),
 axis.title=element_blank(),
 plot.title = element_text(hjust = -0.06),
 plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))})

# Compute the gap statistic
set.seed(sum(utf8ToInt("Kerry Burton")))
gap_stat_p <- lapply(list_DT,function(DF) {clusGap(t(DF), FUN = pam,K.max = 8, B = 500)})

#######  FIGURE S1 - no. clusters ########

g1 <- lapply(gap_stat_p,fviz_gap_stat,maxSE=list(method="globalSEmax", SE.factor = 3))
g2 <- lapply(list_DT,function(DF) {fviz_nbclust(t(DF), pam, method = "wss",k.max = 8)})
g3 <- lapply(list_DT,function(DF) {fviz_nbclust(t(DF), pam, method = "silhouette",k.max = 8)})

g1 <- lapply(g1,function(g) {g +ylab("") +theme_classic_thin(base_size=14) %+replace% theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())})
g2 <- lapply(g2,function(g) {g +ylab("") +theme_classic_thin(base_size=14) %+replace% theme(plot.title = element_blank(),axis.title.x = element_blank())})
g3 <- lapply(g3,function(g) {g +ylab("") +theme_classic_thin(base_size=14) %+replace% theme(plot.title = element_blank())})

g1[[1]] <- g1[[1]] + ggtitle("Both")
g1[[2]] <- g1[[2]] + ggtitle("Experiment 1") + ylab("Gap statistic (k)")
g1[[3]] <- g1[[3]] + ggtitle("Experiment 2")
g2[[2]] <- g2[[2]] + ylab("Total Within SS")
g3[[2]] <- g3[[2]] + ylab("Avg silhouette width")
g3[[3]] <- g3[[3]] + xlab("Number of clusters")

 ggsave("Figure_S1_NEW_ANALYSIS.pdf",grid.arrange(g1[[2]],g1[[3]],g1[[1]],
			g2[[2]],g2[[3]],g2[[1]],
			g3[[2]],g3[[3]],g3[[1]],ncol=3,nrow=3))


set.seed(sum(utf8ToInt("Kerry Burton")))
Clusters <- lapply(list_DT,function(DF) {pam(t(DF),4)[[3]]})

colData <- as.data.table(inner_join(colData,as.data.frame(Clusters$both)%>% mutate(Sample=rownames(as.data.frame(Clusters$both)))))
colData <- colData[order(colData[,3],Group),]
colData$Clusters <- c(rep(2,6),rep(3,5),rep(1,5),4,4)
colData <- colData[order(Clusters,Group,Sample),]
colData$hclust <- as.factor(c(rep(1,6),rep(2,5),rep(3,5),rep(4,2)))

#===============================================================================
#       PCA
#===============================================================================

# covariance PCA (probably best as data is already on a log scale)
DF <- DF[,colData$Sample,with=F]

mypca.cov <- prcomp(t(DF))

mypca.cov$percentVar <- mypca.cov$sdev^2/sum(mypca.cov$sdev^2)

d <- t(data.frame(t(mypca.cov$x)*mypca.cov$percentVar))

# figure 1 B (PCA with medoids)
g1 <- plotOrd(d,colData,design="Clusters",cbPalette=T,pointSize=1.5,axes=c(2,3),alpha=0.75,labels=T,sublabels=c(seq(1,18))[-6])+ stat_ellipse(type="norm",geom="polygon", level=0.85, alpha=0.2)
ggsave("Figure_1_B.pdf",g1)


# figure 2 B - correlation + kmedoids
g2 <- plotOrd(d,colData,design="Group",shapes="Clusters",cbPalette=T,pointSize=1.5,axes=c(2,3),alpha=0.75,label="Sample",sublabels=c(seq(1,18))[-6])+ stat_ellipse(type="norm",geom="polygon", level=0.85, alpha=0.2)

ggsave("Figure_2_B.pdf",g2)

# PCA 1 vs 2 with
ggsave("Figure_S4.pdf",plotOrd(d,colData,design="Group",shape="Clusters",cbPalette=T,pointSize=1.5,axes=c(1,2),alpha=0.75,labels=T,sublabels=c(seq(1,18))[c(-2,-6,-11)])+ stat_ellipse(type="norm",geom="polygon", level=0.85, alpha=0.2))

#===============================================================================
#       NEW FIGURE 1 - Clustering
#===============================================================================

#ggsave("Figure_S1_B.pdf",plot(fit))

#ggsave("Figure_S1_C.pdf",g)

# Figure 1
d <- dist(t(DF), method = "euclidean") # distance matrix
set.seed(sum(utf8ToInt("Kerry Burton")));
fit <- hclust(d, method="ward.D2")
#g_hi <- fviz_dend(fit,main=NULL,ggtheme = theme_classic(base_size=16),lwd=0.5)# I had to hack fviz_dend and set offset_labels to -1
library(ggdendro)


g_hi <- ggdendrogram(fit,rotate=F)+ ggtitle("A") + theme_blank(base_size=16) %+replace%
theme(axis.text.x = element_text(angle = 45, vjust = 1.1,hjust = 1),
 axis.ticks=element_blank(),
 axis.title=element_blank(),
 plot.title = element_text(hjust = -0.06),
 plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))

#set.seed(sum(utf8ToInt("Kerry Burton")))
#pam(t(DF),4,metric="euclidean",stand = F)
#colData$Cluster <- as.factor(c(2,2,3,3,3,3,3,2,2,2,1,4,4,1,1,2,1,1))

#d <- t(data.frame(t(mypca.cov$x)*mypca.cov$percentVar))
g_pca <-plotOrd(obj=d,textSize=14,colData=colData,design="Clusters",cbPalette=T,pointSize=1.5,axes=c(2,3),alpha=0.75,labels=T,sublabels=c(seq(1,18))[c(-6)])+
stat_ellipse(type="norm",geom="polygon", level=0.75, alpha=0.2) +
ggtitle("B") + theme_classic_thin(base_size=14) %+replace% theme(plot.title = element_text(hjust = -0.11),plot.margin = unit(c(-5,5.5,-5,5.5), "pt"))

### FIGURE 1 ####
ggsave("Figure_1.pdf",grid.arrange(g_hi,g_pca,nrow=2,padding = unit(-0.5, "line")),width=7,height=7)

### FIGURE S3
set.seed(sum(utf8ToInt("Kerry Burton")));
fit_clus <- c(2,2,3,3,3,3,3,2,2,2,1,4,4,1,1,1,1,1)
sil1 <- silhouette(fit_clus, dist(scale(t(DF),scale=F)))
rownames(sil1) <- colData$Sample
g <- fviz_silhouette(sil1,label=F,palette=c("#000000", "#E69F00", "#56B4E9", "#009E73"))
g <- g + theme_classic_thin(base_size=14) %+replace%
	theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
	axis.ticks=element_blank(),
	plot.title = element_text(hjust = -0.07),
	axis.title.y=element_blank(),
	plot.margin = unit(c(-5,5.5,-5,5.5), "pt"))
g_sil1 <- g +	scale_y_continuous(expand = c(0, 0), limits = c(-0,0.75)) + ggtitle("A")

sil2 <- silhouette(as.number(colData$Cluster), dist(scale(t(DF),scale=F)))
rownames(sil2) <- colData$Sample
g <- fviz_silhouette(sil2,label=F,palette=c("#000000", "#E69F00", "#56B4E9", "#009E73"))
g <- g + theme_classic_thin(base_size=14) %+replace%
	theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
	axis.ticks=element_blank(),
	plot.title = element_text(hjust = -0.07),
	axis.title.y=element_blank(),
	plot.margin = unit(c(5.5,5.5,-5,5.5), "pt"))
g_sil2 <- g +	scale_y_continuous(expand = c(0, 0), limits = c(-0,0.75)) + ggtitle("B")

sil3 <- silhouette(as.numeric((as.factor(colData$Group))), dist(scale(t(DF),scale=F)))
rownames(sil3) <- colData$Sample
g <- fviz_silhouette(sil3,label=F,palette=c("#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442"))
g <- g + theme_classic_thin(base_size=14) %+replace%
	theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
	axis.ticks=element_blank(),
	plot.title = element_text(hjust = -0.07),
	axis.title.y=element_blank(),
	plot.margin = unit(c(-5,5.5,5.5,5.5), "pt"))
g_sil3 <- g +	scale_y_continuous(expand = c(0, 0), limits = c(-0,0.75)) + ggtitle("C")

ggsave("Figure_S3_NEW.pdf",grid.arrange(sil1,sil2,sil3,nrow=3),width=7,height=9)

#===============================================================================
#      Correlation heatmap (or correlation matrix)
#===============================================================================

cormat <- cor(DF)

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

set.seed(sum(utf8ToInt("Kerry Burton")))
cormat <- reorder_cormat(cormat)

cormat[lower.tri(cormat,diag=F)] <- NA

melted_cormat <- melt(cormat, na.rm = TRUE)

g <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))
g <- g + geom_tile(color = "white")
g <- g + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab",  name="Pearson\nCorrelation")
g <- g + coord_fixed() #+ geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 2.5)
g <- g + geom_hline(yintercept=c(1.5,3.5,8.5,13.5,18.5), linetype="dashed", color = "black")
g <- g + annotate("text", y = c(3,6,11,16), x=1.5, label = c("Group IV", "Group III", "Group II", "Group I"),size=3)
g3 <- g + theme_minimal(base_size=11) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),
axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank())
ggsave("Fig_2A.pdf",g3)

#===============================================================================
#      Figure 2 - Viral correlation
#===============================================================================

#g2_1 <- g3 + annotation_custom(grob=textGrob(label="A",hjust = 0, gp = gpar(cex = 1.5)),ymin=18,ymax=18,xmin=19,xmax=19)
#g1_1 <- ggplotGrob(g2+annotate("text",label=paste("B"), x=-20, y=3,size=6))

g2_1 <- g3 + annotation_custom(grob=textGrob(label="A",hjust = 0, gp = gpar(cex = 1.5)),ymin=18,ymax=18,xmin=-4,xmax=-4)
gt2 <- ggplot_gtable(ggplot_build(g2_1))
gt2$layout$clip[gt2$layout$name == "panel"] <- "off"
#grid.draw(gt2)

g1_1 <- g2 + annotation_custom(grob=textGrob(label="B",hjust = 0, gp = gpar(cex = 1.5)),ymin=5.5,ymax=5.5,xmin=-13.5,xmax=-13.5)
gt1 <- ggplot_gtable(ggplot_build(g1_1))
gt1$layout$clip[gt1$layout$name == "panel"] <- "off"
#grid.draw(gt1)

layout_matrix <- cbind(c(1,1,2),c(1,1,2))

ggsave("Figure_2.pdf",grid.arrange(gt2,gt1,layout_matrix=layout_matrix),width=8,height=9)

#===============================================================================
#       
#	Statistical analysis of virus dCt experiment and clusterd/viruses			    
#			    
#	ANOVA - variance not equal nor are the distributions normal
#
#
#===============================================================================
			    
library(car)

# convert data from wide to long format
anova_data <- melt(countData[,c(-2,-3)])

# add an experiment column to the data
anova_data$experiment <- "E1"
anova_data[Sample %like% "D1"]$experiment <- "E2"

# join data with metadata
anova_data <- left_join(anova_data,colData, by=c("variable"="Sample"))

# rename virus column
names(anova_data)[2] <- "virus"

# remove rows with missing values (shouldn't be any)
anova_data <- anova_data[complete.cases(anova_data),]

# build linear model first and do some test
model <- lm(value~experiment*virus,data=anova_data)
grid.arrange(arrangeGrob(grobs=plot_lm(model),se=F))
plot_lm(model)

# permutation ANOVA of ordinal abundance (removes problem of unequal variance and non normality)


# levene (or  Brown-Forsythe) test for equivelence of variance
model <- leveneTest(value~experiment*Clusters,data=anova_data,center=median)

# then post hoc tests

anova_data <- anova_data %>% group_by(virus,experiment,Clusters) %>% mutate(dat.med = ifelse(value,median(value, na.rm=TRUE), ifelse(value==NA, NA)))


anova_data$dat.med.res<-abs(anova_data$value-anova_data$dat.med)

# Then we run an ANOVA, and post-hoc if necessary:
levene.dat.aov<-aov(dat.med.res~experiment*Clusters,anova_data)
summary(levene.dat.aov)
TukeyHSD(levene.dat.aov)


# run ANOVA model (unbalanced design - using type III anova (interaction is significant))

library(lsmeans)
library(emmeans)

options(contrasts = c("contr.sum", "contr.poly"))

model = lm(dat.med~virus*experiment,anova_data)

Anova(model,type="III")
#TukeyHSD(Anova(lm(value~experiment*Cluster,anova_data),type="III"))
lsmeans(model,pairwise~virus*experiment,adjust="Tukey")

marginal = emmeans(model,~experiment*Cluster)


#fviz_nbclust(t(DF), kmeans, method = "wss") +   geom_vline(xintercept = 4, linetype = 2)
#fviz_nbclust(t(DF), hcut, method = "wss")   +   geom_vline(xintercept = 4, linetype = 2)


#fviz_nbclust(t(DF), kmeans, method = "silhouette") +   geom_vline(xintercept = 4, linetype = 2,colour=blue)
#fviz_nbclust(t(DF), hcut, method = "silhouette",hc_method = "ward.D2")   +   geom_vline(xintercept = 4, linetype = 2)


library(NbClust)
ind <- c("kl", "ch", "hartigan", "cindex", "db", "silhouette", "duda", "pseudot2", "beale",
        "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
        "gamma", "gplus", "tau", "dunn", "hubert", "sdindex",
        "dindex", "sdbw")

nb <- sapply(ind,function(ind) {NbClust(t(DF), distance = "euclidean", min.nc = 2, max.nc = 10, method = "centroid", index =ind)})


#===============================================================================
#      Figure 3 - Mean/median abundance and variance of viral populations
#===============================================================================
library(matrixStats)

# this is actually sd_means or sd)medians
var_means_exp <- as.data.frame(data.table(Sample=colnames(DF),var_1=apply(DF[1:40,],2,sd),var_2=apply(DF[41:89,],2,sd),mean_1=colMeans(40-DF[1:40,]),mean_2=colMeans(40-DF[41:89,])))
# qf <- function(X){x=1/quantile(scale(X))[[4]];mad(X,constant=x)}
#var_means_exp <- as.data.frame(data.table(Sample=colnames(DF),var_1=apply(DF[1:40,],2,qf),var_2=apply(DF[41:89,],2,qf),mean_1=colMedians(40-as.matrix(DF[1:40,])),mean_2=colMedians(40-as.matrix(DF[41:89,]))))

rownames(var_means_exp) <- var_means_exp$Sample
var_means_exp <- var_means_exp[rownames(cormat)[18:1],]
var_means_exp <- var_means_exp[c(1:5,18,6:17),]


var_means <- data.table(Virus=as.factor(c(var_means_exp$Sample,var_means_exp$Sample)),
						Experiment=as.factor(c(rep("Experiment 1",18),rep("Experiment 2",18))),
						Mean=c(var_means_exp$mean_1,var_means_exp$mean_2),
						Variance=c(var_means_exp$var_1,var_means_exp$var_2))
	     
var_means$Virus <- factor(var_means$Virus,levels=colData$Sample)
#var_means$Virus <- factor(var_means$Virus,levels(var_means$Virus)[c(15,14,9,16,17,13,3,10,11,2,1,5:8,18,4,12)])

melted_var_means <- melt(var_means,id.vars=c("Virus","Experiment"))

melted_var_means[variable=="Variance",y_min:=0]
melted_var_means[variable=="Variance",y_max:=62]
melted_var_means[variable=="Mean",y_min:=0]
melted_var_means[variable=="Mean",y_max:=30.1]

melted_means <- melted_var_means[variable=="Mean",]
melted_vars <- melted_var_means[variable=="Variance",]

test <- melted_means[melted_vars,on=c("Virus","Experiment"),]
test$Virus <- factor(test$Virus,levels=rownames(cormat)[18:1]) # this is not correct
test$Virus <- factor(test$Virus,levels=c(
  "ORFan2","ORFan3","ORFan5","ORFan7","MBV","AbV2","AbSV","AbV10","AbV12","AbV6_RNA1","AbV6_RNA2",
  "AbV16_RNA1","AbV16_RNA2","AbV16_RNA3","AbV16_RNA4","ORFAN8","AbV14","AbV9"))


g <- ggplot(data=test,aes(x=Virus,y=value,fill=variable,g=Experiment))
g <- g + geom_bar(stat="identity",colour="white",width=0.7,position = position_dodge(width=0.75) )
g <- g + scale_fill_manual(values = c("Mean" = "orange"))
g <- g + geom_errorbar(aes(ymin=value-i.value, ymax=value+i.value),width=.2,position=position_dodge(.9))
#g <- g + geom_vline(xintercept=c(5.5,6.5,11.5,16.5), linetype="dashed", color = "red") # this is for correlation groups
g <- g + geom_vline(xintercept=c(6,11.5,16.5), linetype="dashed", color = "red")
g <- g + ylab(expression(40 - Delta*"Ct"))
g <- g + facet_wrap(~Experiment,ncol=1,scales="free_y")
g <- g + expand_limits(y=30)
g <- g + theme_blank(base_size=12) %+replace%
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5,linetype=1),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
	axis.ticks=element_blank())

ggsave("NEW_FIG3_MEDIANS.pdf",g)
ggsave("NEW_FIG3_MEANS.pdf",g)

# box plots
DF$Experiment <- c(rep("Experiment 1",40),rep("Experiment 2",49))
t2 <- melt(DF,id.vars="Experiment")
t2$value <- 40 - t2$value
t3 <- colData[t2,on=c("Sample"="variable")]
t3$Sample <- factor(t3$Sample,levels=c(
  "ORFan2","ORFan3","ORFan5","ORFan7","MBV","AbV2","AbSV","AbV10","AbV12","AbV6_RNA1","AbV6_RNA2",
  "AbV16_RNA1","AbV16_RNA2","AbV16_RNA3","AbV16_RNA4","ORFAN8","AbV14","AbV9"))
g <- ggplot(data=t3,aes(x=Sample,y=value))
g <- g + geom_boxplot()
g <- g + geom_vline(xintercept=c(6,11.5,16.5), linetype="dashed", color = "red")
g <- g + ylab(expression(40 - Delta*"Ct"))
g <- g + facet_wrap(~Experiment,ncol=1,scales="free_y")
g <- g + theme_blank(base_size=12) %+replace%
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5,linetype=1),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
	axis.ticks=element_blank())
ggsave("NEW_FIG3_BOXPLOTS.pdf",g)

# other plots???
g <- ggplot(data=melted_means,aes(x=Virus,y=value,fill=variable,g=Experiment))
g <- g + geom_bar(stat="identity",colour="white",width=0.7,position = position_dodge(width=0.75) )
g <- g + scale_fill_manual(values = c("Mean" = "black", "Variance" = "orange"))
#g <- g + geom_vline(xintercept=c(5.5,6.5,11.5,16.5), linetype="dashed", color = "red") # this is for correlation groups
g <- g + geom_vline(xintercept=c(6,11.5,16.5), linetype="dashed", color = "red")
g <- g + ylab(expression(40 - Delta*"Ct"))
g <- g + facet_wrap(~Experiment+variable, nrow = 2,ncol=1,scales="free_y")
g <- g + expand_limits(y=30)
g <- g + theme_blank(base_size=12) %+replace%
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5,linetype=1),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
	axis.ticks=element_blank())

# remove space at bottom of grapha and retain space at top
g <- g + scale_y_continuous(expand = expand_scale(mult = c(0, .075)))

# set different y_limits for each facet
g <- g + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
g1 <- g


g <- ggplot(data=melted_vars,aes(x=Virus,y=value,fill=variable,g=Experiment))
g <- g + geom_bar(stat="identity",colour="white",width=0.7,position = position_dodge(width=0.75) )
g <- g + scale_fill_manual(values = c("Mean" = "black", "Variance" = "orange"))
g <- g + geom_vline(xintercept=c(6,11.5,16.5), linetype="dashed", color = "red")
g <- g + ylab("")
#g <- g + ylab(expression(40 - Delta*"Ct"))
g <- g + facet_wrap(~Experiment+variable, nrow = 2,ncol=1,scales="free_y")
g <- g + expand_limits(y=30)
g <- g + theme_blank(base_size=12) %+replace%
	theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5,linetype=1),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),
	axis.ticks=element_blank())

# remove space at bottom of grapha and retain space at top
g <- g + scale_y_continuous(expand = expand_scale(mult = c(0, .075)))

# set different y_limits for each facet
g <- g + geom_blank(aes(y = y_min)) + geom_blank(aes(y = y_max))
g2 <- g

ggsave("Figure_3_Mean_Abundance.pdf",grid.arrange(g1,g2,ncol=2,nrow=1),width=8,height=8)


	     
#===============================================================================
#      Figure 5 - Cap colour and viral abundance scatter plots
#===============================================================================

# melt the imported data
colour <- melt(countData[,-1],id.vars=c("Experiment","Log(?E)"))

# convert experiment from a numeric to a factor
colour$Experiment <- as.factor(colour$Experiment)

# add colnames
colnames(colour) <- c("Experiment", "Log_dE", "Virus","dCT")
	     
# calculate correlation probablities for both experiments
colData$col_prob_1 <- sapply(seq(4,21),function(i) cor.test(unlist(countData[Experiment=="1",3]),unlist(countData[Experiment=="1",..i]))[[3]])
colData$col_prob_2 <- sapply(seq(4,21),function(i) cor.test(unlist(countData[Experiment=="2",3]),unlist(countData[Experiment=="2",..i]))[[3]])
			     
# subset colour based on significant correlation
col_small <- colour[colour$Virus%in%colData[colData$col_prob_2<=0.05,Sample],]

# convert dCt to 40-
col_small$dCT <- 40-col_small$dCT

# data frame for probabilities
p1 <- c("italic('p')==~0.028",
		"italic('p')==~0.014",
		"italic('p')==~0.022",
		"italic('p')==~0.106",
		"italic('p')==~0.023",
		"italic('p')==~0.878"
)


p2 <- c("italic('p')<~0.001",
		"italic('p')<~0.001",
		"italic('p')<~0.001",
		"italic('p')<~0.001",
		"italic('p')<~0.001",
		"italic('p')<~0.001"
)

dat <- data.frame(x = rep(10, 6), y = rep(1.3, 6),Virus=colData$Sample[levels(colour$Virus)%in%colData[colData$col_prob_2<=0.05,Sample]],p1=p1,p2=p2)

# plot
g <- ggplot(data=col_small,aes(x=dCT,y=Log_dE,colour=Experiment))
g <- g + geom_point(size=2,na.rm = TRUE)+ scale_colour_manual(values=c("black","orange"))

# make expression for greek letter axes labels
g <- g + xlab(expression(40 - Delta*"Ct"))
g <- g + ylab(expression(Delta*"E*"))
g <- g + theme_classic_thin() %+replace% theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))

# add regression lines
g <- g + stat_smooth(method="lm", se=FALSE)

# facet the plot
g <- g + facet_wrap(~Virus, nrow = 3,ncol=2,scales="free_y")

g <- g + geom_text(aes(x, y, label=p1, group=NULL),data=dat,inherit.aes=F,size=2.5,parse = T) +
geom_text(aes(x, y-0.05, label=p2, group=NULL),data=dat,inherit.aes=F,colour="orange",size=2.5,parse = T)

ggsave("Figure_4.pdf",g)


### PER CLUSTER ###

# combine colour with colData	     
colour2 <- as.data.table(left_join(colour,colData,by=c("Virus"="Sample")))
			     
# probabilities per cluster
cor.test(colour2$Log_dE[colour2$Experiment=="1"&colour2$Clusters==1],colour2$dCT[colour2$Experiment=="1"&colour2$Clusters==1])[[3]]
cor.test(colour2$Log_dE[colour2$Experiment=="1"&colour2$Clusters==2],colour2$dCT[colour2$Experiment=="1"&colour2$Clusters==2])[[3]]
cor.test(colour2$Log_dE[colour2$Experiment=="1"&colour2$Clusters==3],colour2$dCT[colour2$Experiment=="1"&colour2$Clusters==3])[[3]]
cor.test(colour2$Log_dE[colour2$Experiment=="1"&colour2$Clusters==4],colour2$dCT[colour2$Experiment=="1"&colour2$Clusters==4])[[3]]
			     
cor.test(colour2$Log_dE[colour2$Experiment=="2"&colour2$Clusters==1],colour2$dCT[colour2$Experiment=="2"&colour2$Clusters==1])[[3]]
cor.test(colour2$Log_dE[colour2$Experiment=="2"&colour2$Clusters==2],colour2$dCT[colour2$Experiment=="2"&colour2$Clusters==2])[[3]]
cor.test(colour2$Log_dE[colour2$Experiment=="2"&colour2$Clusters==3],colour2$dCT[colour2$Experiment=="2"&colour2$Clusters==3])[[3]]
cor.test(colour2$Log_dE[colour2$Experiment=="2"&colour2$Clusters==4],colour2$dCT[colour2$Experiment=="2"&colour2$Clusters==4])[[3]]

#cor.test(colour2$Log_dE[colour2$Experiment=="1"&colour2$hclust==1],colour2$dCT[colour2$Experiment=="1"&colour2$hclust==1])[[3]]
#cor.test(colour2$Log_dE[colour2$Experiment=="2"&colour2$hclust==1],colour2$dCT[colour2$Experiment=="2"&colour2$hclust==1])[[3]]
#cor.test(colour2$Log_dE[colour2$Experiment=="1"&colour2$hclust==2],colour2$dCT[colour2$Experiment=="1"&colour2$hclust==2])[[3]]
#cor.test(colour2$Log_dE[colour2$Experiment=="2"&colour2$hclust==2],colour2$dCT[colour2$Experiment=="2"&colour2$hclust==2])[[3]]

p1 <- c("italic('p')==~0.863",
	"italic('p')==~0.427",
	"italic('p')==~0.044",
	"italic('p')==~0.779"
)

p2 <- c("italic('p')<~0.001",
	"italic('p')==~0.165",
	"italic('p')<~0.001",
	"italic('p')==~0.081"
)
dat <- data.frame(x = rep(10, 4), y = rep(1.2, 4),Clusters=c(1:4),p1=p1,p2=p2)
dat[4,1] <- c(30)

# convert dCt to 40-
colour2$dCT <- 40-colour2$dCT

# make a labeler (saves changing the Cluster factor levels and remaking dat)
relabel <- c(
	`1` = "Cluster 1",
	`2` = "Cluster 2",
	`3` = "Cluster 3",
	`4` = "Cluster 4"
)
# plot
g <- ggplot(data=colour2,aes(x=dCT,y=Log_dE,colour=Experiment))
#g <- g + geom_point(size=2,na.rm = TRUE)+ scale_colour_manual(values=c("black","orange"))
g <- g+ scale_colour_manual(values=c("black","orange"))

# make expression for greek letter axes labels
g <- g + xlab(expression(40 - Delta*"Ct"))
g <- g + ylab(expression(Delta*"E*"))
g <- g + theme_classic_thin() %+replace% theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))

# add regression lines
g <- g + stat_smooth(method="lm", se=T)

# facet the plot and relabel 
g <- g + facet_wrap(~Clusters, nrow = 3,ncol=2,scales="free_y",labeller = as_labeller(relabel))

# 
g <- g + geom_text(aes(x, y, label=p1, group=NULL),data=dat,inherit.aes=F,size=2.5,parse = T) +
geom_text(aes(x, y-0.03, label=p2, group=NULL),data=dat,inherit.aes=F,colour="orange",size=2.5,parse = T)

ggsave("Figure_5_V2.pdf",g)

summary(lm(deltaE~ORFan2+ORFan3+ORFan5+ORFan7+MBV,data=countData[Experiment=="1",]))
summary(lm(deltaE~ORFan2+ORFan3+ORFan5+ORFan7+MBV,data=countData[Experiment=="2",]))
summary(lm(deltaE~AbV2+AbSV+AbV10+AbV12+AbV6_RNA1+AbV6_RNA2,data=countData[Experiment=="1",]))			     
summary(lm(deltaE~AbV2+AbSV+AbV10+AbV12+AbV6_RNA1+AbV6_RNA2,data=countData[Experiment=="2",]))					     
summary(lm(deltaE~AbV16_RNA1+AbV16_RNA2+AbV16_RNA3+AbV16_RNA4+ORFAN8,data=countData[Experiment=="1",]))
summary(lm(deltaE~AbV16_RNA1+AbV16_RNA2+AbV16_RNA3+AbV16_RNA4+ORFAN8,data=countData[Experiment=="2",]))
summary(lm(deltaE~AbV14+AbV9,data=countData[Experiment=="1",]))			     
summary(lm(deltaE~AbV14+AbV9,data=countData[Experiment=="2",]))			     
			     
#===============================================================================
#      Density plots
#===============================================================================

dense <- melt(countData[,c(2,4:21)],id.vars=c("Experiment"))
dense$Experiment <- as.factor(dense$Experiment)
dense$value <- 40 - dense$value
dense$variable <- factor(dense$variable,levels=colData$Sample)

# cut dense into exp1 and exp2
l1 <- lapply(levels(dense[["variable"]]),function(s) dense[variable==s&Experiment=="1",value])
l2 <- lapply(levels(dense[["variable"]]),function(s) dense[variable==s&Experiment=="2",value])

# get bootstrap probability that is double compared to single normal
sink("/dev/null")
 set.seed(sum(utf8ToInt("Kerry Burton")));
 booted_1 <- lapply(l1,boot.comp,mix.type = "normalmix",max.comp=1,B=10000)
 booted_2 <- lapply(l2,boot.comp,mix.type = "normalmix",max.comp=1,B=10000)
sink()

# table S3
data.frame(row.names=levels(dense[["variable"]]),log.lik=sapply(booted_1,function(l) l$obs.log.lik), p=sapply(booted_1,function(l) l$p.values),log.lik=sapply(booted_2,function(l) l$obs.log.lik), p=sapply(booted_2,function(l) l$p.values))

colData$p_1 <- sapply(booted_1,'[[',1)
colData$p_2 <- sapply(booted_2,'[[',1)


#===============================================================================
#      Figure 4 - Best fit models for viral RNAs for each experiment
#===============================================================================

# set various labels
ylabel  <- textGrob("Density", rot = 90, vjust = 0.5,gp = gpar(fontsize = 16))
xlabel  <- textGrob(expression(40 - Delta*"Ct"),hjust=-11,,gp = gpar(fontsize = 14))
exp1_label   <- textGrob("Experiment 1",gp = gpar(fontsize = 16))
exp2_label   <- textGrob("Experiment 2",gp = gpar(fontsize = 16))
ylabel_empty <- textGrob(" ", rot = 90, vjust = 0.5)

# ggplot mixture function
gg.mixEM <- function(EM) {
  require(ggplot2)
  x       <- with(EM,seq(min(x),max(x),len=1000))
  pars    <- with(EM,data.frame(Components=colnames(posterior), mu, sigma,lambda))
  em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
  em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sigma))
  em.df$n <- with(em.df,dnorm(x,mean=mean(EM$x),sd=sd(EM$x)))
  em.df <- rbind(em.df[,c(1,2,6)],data.frame(x=em.df$x,Components="norm",y=em.df$n))
  g <- ggplot(data.frame(x=EM$x),aes(x,y=..density..))
  g <- g + geom_histogram(fill=NA,color="lightgrey")
  g <- g + geom_line(data=em.df,aes(x,y=y,colour=Components),inherit.aes=F) + scale_colour_manual(values=c( "#E69F00", "#56B4E9","#000000"))
  #g <- g + geom_polygon(data=em.df,aes(x,y,fill=comp),colour="grey50", alpha=0.25)
#  g <- g + geom_line(data=em.df,aes(x,n),colour="red")
  g <- g + scale_fill_discrete("Component\nMeans",labels=format(em.df$mu,digits=3))
  g <- g + xlab(NULL)+ ylab(NULL)
  g + theme_minimal() %+replace% theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),panel.grid= element_blank(),legend.position="none")
}

# Experiment 1
set.seed(sum(utf8ToInt("Kerry Burton")));
ll <- lapply(levels(dense[["variable"]]),function(s)  {normalmixEM(dense[variable==s&Experiment=="1",value])})
g <- lapply(seq_along(ll), function(i) gg.mixEM(ll[[i]])+ggtitle(levels(dense[["variable"]])[i])+theme(plot.title = element_text(size = 12)))
# set limits for AbV14
g[[8]] <- g[[8]] + scale_y_continuous(limits = c(0,0.2))
g[[17]] <- g[[17]] + scale_y_continuous(limits = c(0,0.2))

# add red borders for double normal p <=0.05
g[which(colData$p_1<=0.05)] <- lapply(g[which(colData$p_1<=0.05)],function(g) g + theme(panel.border = element_rect(colour = "red", fill=NA, size=1)))
# g[which(sapply(booted_1,function(l) l$p.values<=0.05)==T)] <- lapply(g[which(sapply(booted_1,function(l) l$p.values<=0.05)==T)],function(g) g + theme(panel.border = element_rect(colour = "red", fill=NA, size=1)))

g1       <- g
# g1[c(4)] <- lapply(g1[c(4)],function(g) g + theme(panel.border = element_rect(colour = "black",size=0.5)))
g1       <- lapply(g1,function(g) g  + theme(text = element_text(size=14)))


# Experiment 2
set.seed(sum(utf8ToInt("Kerry Burton")))
ll <- lapply(levels(dense[["variable"]]),function(s)  {normalmixEM(dense[variable==s&Experiment=="2",value])})
g <- lapply(seq_along(ll), function(i) gg.mixEM(ll[[i]])+ggtitle(levels(dense[["variable"]])[i])+theme(plot.title = element_text(size = 12)))

g[[3]] <- g[[3]] + scale_y_continuous(limits = c(0,0.5))
g[[4]] <- g[[4]] + scale_y_continuous(limits = c(0,1.5))
g[[5]] <- g[[5]] + scale_y_continuous(limits = c(0,0.5))

# add red borders for double normal p <=0.05
g[which(colData$p_2<=0.05)] <- lapply(g[which(colData$p_2<=0.05)],function(g) g + theme(panel.border = element_rect(colour = "red", fill=NA, size=1)))

#g[which(sapply(booted_2,function(l) l$p.values<=0.05)==T)] <- lapply(g[which(sapply(booted_2,function(l) l$p.values<=0.05)==T)],function(g) g + theme(panel.border = element_rect(colour = "red", fill=NA, size=1)))
g2       <- g
#g2[c(1)] <- lapply(g2[c(1)],function(g) g + theme(panel.border = element_rect(colour = "red",size=0.5)))
g2       <- lapply(g2,function(g) g  + theme(text = element_text(size=14)))

# get legend
mlegend <- myfunction::get_legend(g1[[1]]+theme_blank(base_size=14) %+replace% theme(legend.position = c(6, 0.5),legend.direction="horizontal"))

# plot graph
gxt <- grid.arrange(
  ylabel_empty, exp1_label,                   ylabel_empty, exp2_label,
  ylabel,       arrangeGrob(grobs=g1,ncol=4), ylabel_empty, arrangeGrob(grobs=g2,ncol=4),
  arrangeGrob(xlabel,mlegend,nrow=2),
  widths=unit.c(unit(2, "lines"), unit(0.5, "npc") - unit(2, "lines"),unit(2, "lines"), unit(0.5, "npc") - unit(2, "lines")),
  heights=unit.c(unit(2, "lines"),unit(1, "npc") - unit(5, "lines"),unit(3, "lines")),
  nrow=3)
dev.off()
ggsave("Figure_4.pdf",grid.draw(gxt),width=16,height=9)


#===============================================================================
#      Figure S1
#===============================================================================

#g <- DF %>% gather() %>% ggplot(aes(value)) +  geom_density()
#g <- g + facet_wrap(~ key, scales = "free",ncol=4)
#g <- g + xlab(expression(Delta*"Ct")) + ylab("Density")
#g <- g + theme_minimal() %+replace% theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),panel.grid= element_blank(),)
#ggsave("density_plots.pdf",g,width=7,height=7)

g <- ggplot(dense,aes(value,fill=Experiment,variable=variable))
g <- g + geom_density(alpha=0.25,kernel="bi")  + scale_x_continuous(labels = function (x) round(x,0)) #+ scale_fill_manual(values=c("black","orange"))
g <- g + facet_wrap(~ variable, scales = "free",ncol=4)
g <- g + xlab(expression(40 - Delta*"Ct")) + ylab("Density")
g <- g + theme_minimal() %+replace% theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),panel.grid= element_blank(),legend.position="right",plot.title = element_text(hjust = 0))
ggsave("Figure_S1_V2.pdf",g,width=7,height=7)


# mixture models
#sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}
#xMix = normalmixEM(dense$value[dense$variable=="AbV10"])




# old stuff

### all graphs ###



gx <- grid.arrange(ylabel,
  arrangeGrob(grobs=g,ncol=4),
  arrangeGrob(xlabel,mlegend,nrow=2),
  widths=unit.c(unit(3, "lines"), unit(1, "npc") - unit(3, "lines")),
  heights=unit.c(unit(1, "npc") - unit(3, "lines"),unit(3, "lines")),
  nrow=2)
dev.off()

ggsave("Figure S3 Experiment 1.pdf",grid.draw(gx))
gx <- grid.arrange(ylabel,
  arrangeGrob(grobs=g,ncol=4),
  arrangeGrob(xlabel,mlegend,nrow=2),
  widths=unit.c(unit(3, "lines"), unit(1, "npc") - unit(3, "lines")),
  heights=unit.c(unit(1, "npc") - unit(3, "lines"),unit(3, "lines")),
  nrow=2)
dev.off()

ggsave("Figure S3 Experiment 2.pdf",grid.draw(gx))

#c2 <- apply(dense,2,function(i)

ll <- lapply(levels(dense[["variable"]]),function(s) lapply(levels(dense[["Experiment"]]), function(ss) {set.seed(sum(utf8ToInt("Kerry Burton")));normalmixEM(dense[variable==s&Experiment==ss,value])}))
names(ll) <- levels(dense[["variable"]])
fwrite(as.data.table((lapply(ll,function(l) data.table(comp1.means=l[3],comp1.sd=l[4],comp2.means=l[12],comp2.sd=l[13])))),"norm.out.txt",sep="\t")

EM.x.1    <- as.data.table(sapply(ll,function(l) l[1]))
EM.x.2    <- as.data.table(sapply(ll,function(l) l[10]))
EM.pars.1 <- as.data.table(sapply(ll,function(l) data.frame(l[1])) with(EM,data.frame(Components=colnames(posterior), mu, sigma,lambda))
EM.pars.2 <- as.data.table(sapply(ll,function(l) l[1]))

set.seed(sum(utf8ToInt("Kerry Burton")))
g1 <- gg.mixEM(normalmixEM(dense[variable=="ORFan2"&Experiment=="1",value],lambda = .5)) + ggtitle("ORFan2, Experiment 1")
set.seed(sum(utf8ToInt("Kerry Burton")))
g2 <- gg.mixEM(normalmixEM(dense[variable=="ORFan2"&Experiment=="2",value])) + ggtitle("ORFan2, Experiment 2")
set.seed(sum(utf8ToInt("Kerry Burton")))
g3 <- gg.mixEM(normalmixEM(dense[variable=="AbV10"&Experiment=="1",value])) + ggtitle("AbV10, Experiment 1")
set.seed(sum(utf8ToInt("Kerry Burton")))
g4 <- gg.mixEM(normalmixEM(dense[variable=="AbV10"&Experiment=="2",value])) + ggtitle("AbV10, Experiment 2")
set.seed(sum(utf8ToInt("Kerry Burton")))
g5 <- gg.mixEM(normalmixEM(dense[variable=="AbV16_RNA1"&Experiment==1,value])) + ggtitle("AbV16_RNA1, Experiment 1")
set.seed(sum(utf8ToInt("Kerry Burton")))
g6 <- gg.mixEM(normalmixEM(dense[variable=="AbV16_RNA1"&Experiment==2,value])) + ggtitle("AbV16_RNA1, Experiment 2")


ylabel  <- textGrob("Density", rot = 90, vjust = 0.5)
#xlabel  <- textGrob(expression(40 - Delta*"Ct"),hjust=0)
xlabel  <- textGrob(expression(40 - Delta*"Ct"),hjust=-5.5)
mlegend <- myfunction::get_legend(g[[1]]+theme_blank() %+replace% theme(legend.position = c(4, 0.5),legend.direction="horizontal"))

gt <- grid.arrange(ylabel,
 arrangeGrob(g1,g2,g3,g4,g5,g6,ncol=2),
 arrangeGrob(xlabel,mlegend,nrow=2),
 widths=unit.c(unit(3, "lines"), unit(1, "npc") - unit(3, "lines")),
 heights=unit.c(unit(1, "npc") - unit(3, "lines"),unit(3, "lines")),
 nrow=2)
ggsave("Figure S3 V2.pdf",gt)


booted_2 <- list(
booted_2[[15]],
booted_2[[14]],
booted_2[[16]],
booted_2[[18]],
booted_2[[17]],
booted_2[[11]],
booted_2[[10]],
booted_2[[1]],
booted_2[[2]],
booted_2[[9]],
booted_2[[8]],
booted_2[[4]],
booted_2[[3]],
booted_2[[5]],
booted_2[[7]],
booted_2[[6]],
booted_2[[13]],
booted_2[[12]]
)

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
librabry(cluster)

load_all("~/pipelines/metabarcoding/scripts/myfunctions")

#===============================================================================
#       Load data
#===============================================================================

countData <- fread("countData")
colData <- fread("colData")

DF <- countData[,4:21]

#===============================================================================
#       ANOVA
#===============================================================================

# convert data from wide to long format
anova_data <- melt(DF)

# add an experiment column to the data
anova_data$experiment <- "E1"
anova_data[Sample %like% "D1"]$experiment <- "E2"

# join data with metadata
anova_data <- left_join(anova_data,colData, by=c("variable"="Sample"))

# rename virus column
names(anova_data)[2] <- "virus"

# remove rows with missing values
anova_data <- anova_data[complete.cases(anova_data),]

# remove MBV as it is not a group
anova_data_minus_MBV <- anova_data[anova_data$group!="X",]

# run ANOVA model
aov(value~experiment*group,anova_data_minus_MBV)
summary(aov(value~experiment*group,anova_data_minus_MBV))
TukeyHSD(aov(value~experiment*group,anova_data_minus_MBV))

#===============================================================================
#       Cluster analysis
#===============================================================================

mylen <- nrow(t(DF)) - 1

km_list <- lapply(seq(1,10),function(i) kmeans(t(DF),centers=i,iter.max=1000)

#Plot the original dataset
# plot(DF$x,DF$y,main="Original Dataset")

#Scree plot to deterine the number of clusters
wss <- (mylen)*sum(apply(t(DF),2,var))
sapply(seq(2,mylen),function(i) wss[i] <<- sum(kmeans(t(DF),centers=i)$withinss))
plot(1:17, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# Ward Hierarchical Clustering
d <- dist(t(DF), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D2") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

#Silhouette analysis for determining the number of clusters
asw <- sapply(seq(2,mylen ),function(i) pam(t(DF), i)$silinfo$avg.width)

k.best <- which.max(asw)
cat("silhouette-optimal number of clusters:", k.best, "\n")
pdf("silhouette.pdf", height=6,width=8)
plot(pam(t(DF), 4))
plot(pam(d, 5))
dev.off()

      
# K-Means Cluster Analysis
fit <- kmeans(t(DF),4)
# get cluster means 
aggregate(t(DF),by=list(fit$cluster),FUN=mean)
# append cluster assignment
DF <- data.frame(DF, clusterid=fit$cluster)
plot(DF$x,DF$y, col = fit$cluster, main="K-means Clustering results")

plot(mypca.cov$x[,2],mypca.cov$x[,3],col = fit$cluster, main="K-means Clustering results")


# NEW ANALYSIS
res <- get_clust_tendency(t(DF), 17, graph = F)
# Hopskin statistic
res$hopkins_stat
0.2440934

library("cluster")
# Compute the gap statistic
gap_stat <- clusGap(DF, FUN = kmeans, nstart = 25,K.max = 10, B = 500) 
# plot  (comes out as 3 - 6)
fviz_gap_stat(gap_stat)

set.seed(sum(utf8ToInt("Kerry Burton")))
#km.res <- lapply( seq(3,15,by=4),function(i) kmeans(t(DF), 5, nstart = i))

km.res <- lapply( seq(3,6)),function(i) kmeans(t(DF), i, nstart = 25))

g <- lapply(km.res,fviz_cluster,t(DF),geom = c("point","text"),repel=T,main=NULL)
lapply(g,function(o) o+coord_fixed())

#===============================================================================
#       PCA
#===============================================================================

# covariance PCA (probably best as data is already on a log scale)
mypca.cov <- prcomp(t(DF))

mypca.cov$percentVar <- mypca.cov$sdev^2/sum(mypca.cov$sdev^2)

set.seed(sum(utf8ToInt("Kerry Burton")))
km <- kmeans(t(DF),4,nstart=25)
#km <- mkeans(mkypca.cov$x,centers=4,iter.max=1000,nstart=25)
colData$Cluster <- as.factor(c(2,2,3,3,3,3,3,2,2,2,1,4,4,1,1,2,1,1))

d <- t(data.frame(t(mypca.cov$x)*mypca.cov$percentVar))


g1 <-  plotOrd(d,colData,design="Group",shape="Cluster",cbPalette=T,pointSize=1.5,axes=c(2,3),alpha=0.75,labels=T,sublabels=c(seq(1,18))[-16])+ stat_ellipse(type="norm",geom="polygon", level=0.95, alpha=0.2)

ggsave("Figure_1_B.pdf",g1)

d <- mypca.cov$x

ggsave("pca.pdf",plotOrd(d,colData,design="group",pointSize=1.5,axes=c(1,2),alpha=0.75))  + stat_ellipse(geom="polygon", level=cluster, alpha=0.2)
ggsave("nogroup.pca.pdf",plotOrd(d,colData,pointSize=1.5,axes=c(1,2),alpha=0.75,cluster=0.95,centers=4))


# correlation PCA
mypca.cor <- prcomp(t(countData[,4:21]),scale=T)
mypca.cor$percentVar <- mypca.cor$sdev^2/sum(mypca.cor$sdev^2)

d <-t(data.frame(t(mypca.cor$x)*mypca.cor$percentVar))
ggsave("corrected.cor_pca.pdf",plotOrd(d,colData,design="group",pointSize=1.5,axes=c(1,2),alpha=0.75,cluster=0.95,centers=4))
ggsave("corrected.cor_pca.2_3.pdf",plotOrd(d,colData,design="group",pointSize=1.5,axes=c(2,3),alpha=0.75,cluster=0.95,centers=4))

d <- mypca.cor$x
ggsave("cor_pca.pdf",plotOrd(d,colData,design="group",pointSize=1.5,axes=c(1,2),alpha=0.75,cluster=0.95,centers=4))

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
g2 <- g + theme_minimal(base_size=11) %+replace% theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),
axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank())
ggsave("Fig_1A.pdf",g2)

#===============================================================================
#      Figure 1
#===============================================================================

g2_1 <- g2 + annotation_custom(grob=textGrob(label="A",hjust = 0, gp = gpar(cex = 1.5)),ymin=18,ymax=18,xmin=19,xmax=19)
g1_1 <- ggplotGrob(g1+annotate("text",label=paste("B"), x=-20, y=3,size=6))

g2_1 <- g2 + annotation_custom(grob=textGrob(label="A",hjust = 0, gp = gpar(cex = 1.5)),ymin=18,ymax=18,xmin=-4,xmax=-4)
gt2 <- ggplot_gtable(ggplot_build(g2_1))
gt2$layout$clip[gt2$layout$name == "panel"] <- "off"
#grid.draw(gt2) 

g1_1 <- g1 + annotation_custom(grob=textGrob(label="B",hjust = 0, gp = gpar(cex = 1.5)),ymin=3.5,ymax=3.5,xmin=-12.5,xmax=-12.5)
gt1 <- ggplot_gtable(ggplot_build(g1_1))
gt1$layout$clip[gt1$layout$name == "panel"] <- "off"
#grid.draw(gt1) 

layout_matrix <- cbind(c(1,1,2),c(1,1,2))

ggsave("Figure_1_NEW.pdf",grid.arrange(gt2,gt1,layout_matrix=layout_matrix),width=8,height=9)


#===============================================================================
#      Figure 2 
#===============================================================================

var_means_exp <- as.data.frame(data.table(Sample=colnames(DF),var_1=apply(DF[1:40,],2,var),var_2=apply(DF[41:89,],2,var),mean_1=colMeans(40-DF[1:40,]),mean_2=colMeans(40-DF[41:89,])))
rownames(var_means_exp) <- var_means_exp$Sample
var_means_exp <- var_means_exp[rownames(cormat)[18:1],]
var_means_exp <- var_means_exp[c(1:5,18,6:17),]

var_means <- data.table(Virus=as.factor(c(var_means_exp$Sample,var_means_exp$Sample)),
						Experiment=as.factor(c(rep("Experiment 1",18),rep("Experiment 2",18))),
						Mean=c(var_means_exp$mean_1,var_means_exp$mean_2),
						Variance=c(var_means_exp$var_1,var_means_exp$var_2))
var_means$Virus <- factor(var_means$Virus,levels(var_means$Virus)[c(15,14,9,16,17,13,3,10,11,2,1,5:8,18,4,12)])						

melted_var_means <- melt(var_means,id.vars=c("Virus","Experiment"))

melted_var_means[variable=="Variance",y_min:=0]
melted_var_means[variable=="Variance",y_max:=62]
melted_var_means[variable=="Mean",y_min:=0]
melted_var_means[variable=="Mean",y_max:=30.1]

melted_means <- melted_var_means[variable=="Mean",]
melted_vars <- melted_var_means[variable=="Variance",]


g <- ggplot(data=melted_means,aes(x=Virus,y=value,fill=variable,g=Experiment))
g <- g + geom_bar(stat="identity",colour="white",width=0.7,position = position_dodge(width=0.75) )
g <- g + scale_fill_manual(values = c("Mean" = "black", "Variance" = "orange"))
g <- g + geom_vline(xintercept=c(5.5,6.5,11.5,16.5), linetype="dashed", color = "red") 
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
g <- g + geom_vline(xintercept=c(5.5,6.5,11.5,16.5), linetype="dashed", color = "red") 
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
     
ggsave("Figure_2_NEW.pdf",grid.arrange(g1,g2,ncol=2,nrow=1),width=8,height=8)





#===============================================================================
#      Figure 3 
#===============================================================================

# melt the imported data
colour <- melt(countData[,-1],id.vars=c("Experiment","Log(?E)"))))

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
			     
ggsave("Figure_3.pdf",g)

#===============================================================================
#      Figure S1 
#===============================================================================
ggsave("Figure_S1_A.pdf",plotOrd(d,colData,design="Group",shape="Cluster",cbPalette=T,pointSize=1.5,axes=c(1,2),alpha=0.75,labels=T,sublabels=c(seq(1,18))[c(-2,-15,-16)])+ stat_ellipse(type="norm",geom="polygon", level=0.85, alpha=0.2))

ggsave("Figure_S1_B.pdf",plot(fit))
			     
sil <- silhouette(as.number(colData$Cluster), dist(scale(t(DF),scale=F)))
rownames(sil) <- colData$Sample
g <- fviz_silhouette(sil,label=T,palette=c("#000000", "#E69F00", "#56B4E9", "#009E73"))
g <- g + theme_classic_thin(base_size=16) %+replace% 
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1),
	axis.ticks=element_blank()) + scale_y_continuous(expand = c(0, 0), limits = c(-0,0.75))

ggsave("Figure_S1_C.pdf",g)

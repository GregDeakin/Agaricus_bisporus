#===============================================================================
#       Load libraries
#===============================================================================

library(limma)
library(virtualArray)
library(dplyr)
library(affyPLM)


#===============================================================================
#       Load one colour Agilent data 
#===============================================================================

targets_mcb<- readTargets("targets_mcb.txt")
targets_phd<- readTargets("targets_phd.txt")

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

RG_mcb <- get_arrays(targets_mcb)
RG_phd <- get_arrays(targets_phd)

RG_norm_mcb <- backgroundCorrect(RG_mcb, method="normexp", offset=50)
RG_norm_phd <- backgroundCorrect(RG_phd, method="normexp", offset=50)

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


f1 <-  filter_run(RG_norm$other$gIsWellAboveBG) 

filt <- f1>0


RG_filt <- RG_norm[filt,]
RG_filt <- RG_filt[RG_filt$genes$ControlType==0,]
RG_ess <- new("MAList", list(targets=RG_filt$targets, genes=RG_filt$genes, source=RG_filt$source,  A=RG_filt$G))




#===============================================================================
#	    Technical Replicates
#===============================================================================

# below for v4 microarray only
RG_ess$genes$SystematicName[grep("Agabi_varbisH97_2$",RG_ess$genes$SystematicName)] <- RG_ess$genes$ProbeName[grep("Agabi_varbisH97_2$",RG_ess$genes$SystematicName)]


ff <- function(X) {
	while ( (max(X) - min(X)) > 2 ) 
		if ((max(X)-median(X))>0.5) 
			X[which.max(X)]<-median(X) 
		else 
			X[which.min(X)] <- median(X)
	
	return(median(X))
}


test <- aggregate(RG_ess$A,by=list(RG_ess$genes$SystematicName),ff)

E.avg <- new("MAList", list(targets=RG_filt$targets, genes=test$Group.1,source=RG_filt$source,  A=test[,2:31]))

# for V5 only
E.avg$genes <- sub(".*_","",E.avg$genes)

E.avg.v4 <- E.avg
E.avg.v5 <- E.avg

#===============================================================================
#	    Merge arrays
#===============================================================================


annotations <- read.csv("phd_genes.txt",header=TRUE)
annotations <- annotations[!duplicated(annotations$Protein.ID),]
gene_names <- read.table("gene_names.txt",sep="\t",header=TRUE) 

annotations <- merge(gene_names,annotations, by.x="JGI_ID",by.y="Protein.ID",all.x=TRUE)

genes_v4 <- inner_join(data_frame(Name=E.avg.v4$genes),annotations[,c(1,6)])
genes_v5 <- inner_join(data_frame(ENA_ID=E.avg.v5$genes),annotations[,1:2])

genes_to_keep <- inner_join(genes_v4,genes_v5)

test<- cbind(E.avg.v4$genes,E.avg,v4$A)
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
sampleNames(merged_norm) <- sample_info[, 1]

merged_combat <- merged_norm
exprs(merged_combat) <- virtualArrayComBat(expression_xls = exprs(merged_combat),sample_info_file=sample_info)

pData(merged_combat) <- as.data.frame(sample_info)
merged_combat_norm <- normalize.ExpressionSet.quantiles(merged_combat)
E.avg$t2 <- exprs(merged_combat_norm)



#===============================================================================
#       Load libraries
#===============================================================================

library(limma)

#===============================================================================
#       Load one colour Agilent data 
#===============================================================================

targets<- readTargets() # will read in targets.txt by default
array_path = "../data"

RG <- read.maimages(
	targets, 
	path=array_path, 
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

#===============================================================================
#      Normalisation
#===============================================================================

RG_norm <- backgroundCorrect(RG, method="normexp", offset=50)

RG_norm$G <- normalizeBetweenArrays(RG_norm$G, method="quantile")

RG_norm$G <- log2(RG_norm$G)

#===============================================================================
#	Data filtering
#===============================================================================

# describes how the data is going to be filtered. 
filter <- function(X) {
	Y <- rowMeans(X)>0.5 # more than half of the samples must be well above background
	Y <- Y==1 # doen't see the point of this, Y is already logical?
	return(Y)
}

# pass biological replicates to the filter seperately
filter_fun <- function(X) {
	Y <- filter(X[,1:4])+ filter(X[,5:8])+ filter(X[,9:12])+ filter(X[,13:16])+ filter(X[,17:20])+ filter(X[,21:24])
	return (Y)
}

f1 <-  filter_fun(RG_norm$other$gIsWellAboveBG) 

filt <- f1>0 # test that at least one of the contrasts/experimental conditions passes the filter


RG_filt <- RG_norm[filt,]
RG_filt <- RG_filt[RG_filt$genes$ControlType==0,]
RG_ess <- new("MAList", list(targets=RG_filt$targets, genes=RG_filt$genes, source=RG_filt$source,  A=RG_filt$G))

#===============================================================================
#	No filter
#===============================================================================

#RG_ess <- new("MAList", list(targets=RG_norm$targets, genes=RG_norm$genes, source=RG_norm$source,  A=RG_norm$G))

#===============================================================================
#	    Technical Replicates
#===============================================================================

ff <- function(X) {
	while ( (max(X) - min(X)) > 2 ) 
		if ((max(X)-median(X))>0.5) 
			X[which.max(X)]<-median(X) 
		else 
			X[which.min(X)] <- median(X)
	
	return(median(X))
}

test <- aggregate(RG_ess$A,by=list(RG_ess$genes$SystematicName),ff)

E.avg <- new("MAList", list(targets=RG_filt$targets, genes=test$Group.1,source=RG_filt$source,  A=test[,2:25]))
E.avg$genes <- sub(".*_","",E.avg$genes)

#===============================================================================
#	    Statistical Analysis
#===============================================================================

f <- factor(targets$Condition, levels = unique(targets$Condition))

design <- model.matrix(~0 + f)

colnames(design) <- levels(f)

fit <- lmFit(E.avg$A, design)

contrast.matrix <- makeContrasts(
	"cond1-control",
	"cond2-control",
	levels=design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


#===============================================================================
#	    Data analysis
#===============================================================================

cond_1 <-topTable(fit2, adjust="BH", coef="cond1-control", genelist=E.avg$genes, number=nrow(RG_filt))
cond_2 <-topTable(fit2, adjust="BH", coef="cond2-control", genelist=E.avg$genes, number=nrow(RG_filt))

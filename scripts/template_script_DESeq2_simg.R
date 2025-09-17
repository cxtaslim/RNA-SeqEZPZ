################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### March 20th, 2018
### designed to be executed with SARTools 1.6.3
###
### modified by CT for run inside singularity image
### v2: this version uses library Polychrome to auto generate colors
### ran with:
### singularity shell --bind project1:/mnt --bind scripts:/scripts centos7-ct_4h.sif
### salloc --partition=himem --ntasks=48
### sartools_env
### export OMP_NUM_THREADS=48
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################

# removing user library path
if(length(.libPaths())>1) .libPaths(.libPaths()[-1])
options(echo=TRUE)

# getting arguments from command line
args = commandArgs(trailingOnly=TRUE)

# parse command line arguments
for(i in 1:length(args)){
  tmp.args <- unlist(strsplit(args[i],"="))
  switch(tmp.args[1],
  padj={padj.cutoff=as.numeric(tmp.args[2])},
  email={email=tmp.args[2]},
  batch_adjust={batch_adjust=tmp.args[2]}
 )}

workDir <- "/mnt/outputs/diff_analysis_rslt"
dir.create(file.path(workDir), showWarnings = FALSE)
scriptDir <- "/scripts/internal_R"
setwd(workDir)

### loading libraries
library(SARTools)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(kableExtra)
library(ggplot2)
forceCairoGraph <- FALSE
if (forceCairoGraph) options(bitmapType="cairo")
library(fdrtool)
library(tibble)
library(EnsDb.Hsapiens.v86)
library(Polychrome)
library(scales) # hue_pal()

# using modified exportResults.DESeq2 to output gene names
source(file.path(scriptDir,'exportResults.DESeq2.R'))
assignInNamespace("exportResults.DESeq2",exportResults.DESeq2,ns="SARTools")

# using modified pairwiseScatterPlots to plot correlation matrices
source(file.path(scriptDir,'pairwiseScatterPlots.R'))
assignInNamespace("pairwiseScatterPlots",pairwiseScatterPlots,ns="SARTools")

# using modified loadTargetFile to allow for punctuation in groupname
source(file.path(scriptDir,'loadTargetFile.R'))
assignInNamespace("loadTargetFile",loadTargetFile,ns="SARTools")

# the following might be my own new function that doesn't exist in SARTools
source(file.path(scriptDir,'getDiffTop.R'))
source(file.path(scriptDir,'calcDiffTotal.R'))
source(file.path(scriptDir,'diffPlot.R'))
source(file.path(scriptDir,'run.DESeq2.R'))
assignInNamespace("run.DESeq2",run.DESeq2,ns="SARTools")
source(file.path(scriptDir,'rawpHist.R'))
assignInNamespace("rawpHist",rawpHist,ns="SARTools")
source(file.path(scriptDir,'summarizeResults.DESeq2.R'))
assignInNamespace("summarizeResults.DESeq2",summarizeResults.DESeq2,ns="SARTools")
source(file.path(scriptDir,'descriptionPlots.R'))
assignInNamespace("descriptionPlots",descriptionPlots,ns="SARTools")
# put the code below as there is error "could not find function'diagSizeFactorsPlots'" when
# calling summarizeResults
environment(summarizeResults.DESeq2)<-asNamespace('SARTools')
source(file.path(scriptDir,'exploreCounts.R'))
assignInNamespace("exploreCounts",exploreCounts,ns="SARTools")
source(file.path(scriptDir,'PCAPlot.R'))
assignInNamespace("PCAPlot",PCAPlot,ns="SARTools")
source(file.path(scriptDir,'MAPlot.R'))
assignInNamespace("MAPlot",MAPlot,ns="SARTools")
source(file.path(scriptDir,'volcanoPlot.R'))
assignInNamespace("volcanoPlot",volcanoPlot,ns="SARTools")
source(file.path(scriptDir,'writeReport.DESeq2.R'))
assignInNamespace("writeReport.DESeq2",writeReport.DESeq2,ns="SARTools")

projectName <- "RNA-seq_differential_analysis" # name of the project
# figures folder name will be appended with folder.name
folder.name="Differential_analysis_result"
author <- email                     # author of the statistical analysis/report

### creating file with experiment info
dir <- file.path('../STAR_2pass','Pass2')
tmp <- read.table("../../samples.txt",comment.char="#")
groupname=tmp[,1]
controlname=tmp[,2]
repname=tmp[,3]
label=paste0(groupname,repname)
files <- paste0(groupname,"_",repname,"_counts.txt")
# use idx.filt to filter out certain samples
idx.filt=1:length(label)
label <- label[idx.filt]
files = files[idx.filt]

# specify subset of analysis to run
subset=FALSE;
# specify the group you want to compare
# sub.names=c('*EWS_ETV4*|*iEF_empty*')
# skip iEF EWS ETV4
sub.names=c('*empty*|*FUS*|*wtEF714*')
## only include R2 thru R4 but skip iluc empty r2
#sub.names="^iEF.*R[2|3|4]|^iLuc_empty.*R[3|4]"
## only include R2 thru R4 iEF_* samples
#sub.names="^iEF.*R[2|3|4]"

#sub.names=c('iEF_empty197_R2|iEF_empty197_R3|iEF_empty197_R4|iEF_EF85_R2|iEF_EF85_R3|iEF_EF85_R4')
#iEF_empty197|iEF_wtEF714|iLuc_empty197')
if(subset){
  files = files[grep(sub.names,files)]
  label = label[grep(sub.names,label)]
}

group <- groupname
rep <- repname
write.table(data.frame(label,files,group,rep),file=file.path('target.txt'),
            quote=FALSE,row.names = FALSE,sep="\t")
targetFile <- file.path("target.txt")     # path to the design/target file
rawDir <- dir # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

varInt <- "group"                                    # factor of interest
# reference biological condition
# condRef is just used to determine reference for calculating fold-change
# changed to get all controname not only the first one.
condRef <-  controlname[!is.na(controlname)]
# blocking factor: NULL (default) or "batch" for example
if (tolower(batch_adjust) == "yes"){
 batch <- "rep"
}else{
 batch <- NULL
 batchRem=FALSE
}
fitType <- "local"                                   # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- padj.cutoff                                # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"
typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors


################################################################################
###                             running script                               ###
################################################################################

# loading target file
# NOTE: SARTools does not allow group names to have underscores.
# so I modified loadTargetFile to allow for underscores
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# vector of colors for each group or condition
ngroup=length(levels(target[,varInt]));
if(ngroup<=5){
	#colors<-c("#ADD8E6","#FFA500","#C71585")
	colors <- scales::hue_pal()(ngroup)
}else if (ngroup<=12){
	# library(colorspace)
	# colors<-qualitative_hcl(ngroup,palette="Dark 3")
	# library(Polychrome)
	# colors=as.vector(dark.colors(ngroup))
	colors=brewer.pal(ngroup,"Paired")

}else{
	# library(Polychrome)
	set.seed(701138) 
	colors <-as.vector(createPalette(40, c("#ADD8E6","#FFA500","#C71585"),M=10000))
}

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
if(!is.null(batch)){
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors,
              batch=batch,varInt=varInt,batchRem=TRUE)
}
# plot PCA without removing batch
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors,
 batch=batch,varInt=varInt,batchRem=FALSE)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot,
# diffPlot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)
# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
 majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
 targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
 condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
 independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
 typeTrans=typeTrans, locfunc=locfunc, colors=colors)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# # Renamed figures folder
# file.rename("figures",paste0("figures_",folder.name))
# file.rename("tables",paste0("tables_",folder.name))

## Author details
## Script name: Script_LNB_pipeline_v1.R
## Purpose of the script: To analyze the samples using minfi and ChAMP
##
## Author(s): Shumaila Sayyab, Ph.D.
## Date Created: 2022-08-16
## Date Modified: 2022-08-16
## Copyright statement: MIT open license
## Contact information: shumaila.sayyab@liu.se
## Please cite: @article{,
##  title={LNB_Study},
##  author={},
##  journal={},
##  year={2022},
##  publisher={Cold Spring Harbor Laboratory Press}}
## Notes: The user needs to have the IDAT files.
##
## github link:


####################################################################

## Pipeline for DNA methylation analysis

# 1. Processing of IDAT files with minfi
# 2. Preparing the data and Pheno files and Annotations
# 3. Minfi QC step: Normalization, quality check, filtering (poor Quality probes/samples, chr XY, probes over SNPs, cross reactive probes)
# 4. SVD analysis using ChAMP
# 5. Batch affect correction using Combat
# 6. Density plots, MDS plots, PCA plots to visualize the data (before and after QC/correction)
# 7. Co-variates / confounder identification
# 8. Cell type adjustments and analysis
# 9. DMC/DMGs identification
#10. Volcano, Heatmaps and Pathway Enrichment (ORA)/Network analysis
#11.Estimated genomic inflation and bais using BACON package in R.


##
setwd("~/OneDrive - Linköpings universitet/LNB_Anna_2022/data/BEA20P070_ML")
pack_R <- c("minfi", "gmqn", "ChAMP", "bnstruct","missMethyl","limma","RColorBrewer","stringr","ggpubr","ggrepel","EnhancedVolcano","calibrate","FlowSorted.Blood.EPIC")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
col1=c("LNB"="maroon","non-LNB"="royalblue")
set.seed(777)

#####################################################################
# 1. Processing of IDAT files with minfi                           #
#####################################################################

#--------------------------- loading the rgSet object or read in the raw dataset
# set the data directory
dataDir <- setwd(getwd())
list.files(dataDir)
baseDir<-"."
targets<-read.metharray.sheet(baseDir,pattern="BEA20P070_ML_sample_sheet.csv")
# read in the raw data from the IDAT files; warnings can be ignored.
rgSet <- read.metharray.exp(targets = targets)
sampleNames(rgSet) <- targets$Sample_Name
getManifest(rgSet)
colnames(rgSet)
#save(rgSet,file="../../Analysis_220816/rgSet_LNB.RData")
## or load the data already run through above steps
#load("rgSet_LNB.RData")
#
############################################################
## 2. Preparing the data and Pheno files, Annotations Hg38
#############################################################

setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")
## subset dataset
dim(rgSet)
colnames(rgSet)
rgSet_LNB<-rgSet[,c(19:32)] ## remove unknown samples (not part of LNB project)
dim(rgSet_LNB)
colnames(rgSet_LNB)

sampleNames(rgSet_LNB)
pData(rgSet_LNB)
MET_pd <- read.table("../Samples/Clinical data standardized 220707_R.csv",header=T,sep=",")

## ################################

#-------------------------- Prepare the Pheno file

## Matching the row names of excel file to the column names of methylation data in rgSet
match(colnames(rgSet_LNB),MET_pd$SampleName)
## renaming two of the samples mislabeled in the samplesheet files
colnames(rgSet_LNB)[10:11]<-c("F22a","F28a")
match(colnames(rgSet_LNB),MET_pd$SampleName)
samples_info<-MET_pd[match(colnames(rgSet_LNB),MET_pd$SampleName),]
match(colnames(rgSet_LNB),samples_info$SampleName)
# done

rgSet_LNB$Age<-samples_info$Age
rgSet_LNB$Sex<-samples_info$Sex
rgSet_LNB$Sample_Group<-samples_info$Sample_Group

samples_info$Slide<-rgSet_LNB$Slide
samples_info$Array<-rgSet_LNB$Array

##################### SAVE Working data files
save(rgSet_LNB,file="rgSet_LNB_processed.RData")
write.table(samples_info,"samples_info_LNB_220816.csv",sep=",",quote = F,row.names = F)

#####################  START HERE ******************

setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")
load("rgSet_LNB.RData")
rgSet<-NULL
#rgSet<-rgSet_LNB
samples<-read.table("samples_info_LNB_220816.csv",header=T,sep=",")
getManifest(rgSet)
pData(rgSet_LNB)
colnames(rgSet_LNB)

## ##############

annotate_hg38<-function(rgSet)
{
  #BiocManager::install("devtools")
  library(devtools)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
  #rgSet@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
  # get the 850k annotation data
  ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
 
} ##end of annotation

############################################################
## 3. Pipeline Minfi QC steps
#############################################################
minfi_pipeline<-function(rgSet,samples)
{
  mycell<-"LNB"
  #-------------- Quality Control of the dataset
  
  detP <- detectionP(rgSet)
  qcReport(rgSet, sampNames=samples$Sample_Name, sampGroups=samples$Group,
           pdf="qcReport_1.pdf")
  keep <- colMeans(detP) < 0.05
  
  rgSet <- rgSet[,keep]
  table(keep)
  ## remove poor quality samples from targets data
  samples <- samples[keep,]
  dim(samples)
  samples[,1:5]
  
  ## remove poor quality samples from detection p-value table
  detP <- detP[,keep]
  table(keep)
  
  # Normalization
  mSetRaw <- preprocessRaw(rgSet)
  mSetSq <- preprocessSWAN(rgSet) ## 1. SWAN normalization, used here
  #mSetSq2 <- preprocessQuantile(rgSet) ## 2. Quantile normalization generated same results
  
  # Visualise what the data looks like before and after normalisation
  
  ##----------------------Plot density
  
  NormPlot <- function(normPlot){
    
    par(mfrow=c(1,2))
    densityPlot(rgSet, sampGroups=samples$Sample_Group,main="Raw", legend=FALSE,pal = col1)
    legend("top", legend = levels(factor(samples$Sample_Group)),
           text.col=col1)
    densityPlot(getBeta(mSetSq), sampGroups=samples$Sample_Group,
                main="SWAN.Normalized", legend=FALSE,pal=col1)
    legend("top", legend = levels(factor(samples$Sample_Group)),
           text.col=col1)
  }
  name<-paste0(mycell,"_Raw_Density_1.png")
  png(name,height = 10, width = 10, units = "in", res = 300)
  NormPlot()
  dev.off()
  dir()
  ##---------------------Plot MDS customized function
  MDSplot <-function(MDSplot){
    par(mfrow=c(1,1))
    #col1=c("HC"="mediumseagreen","CT"="orange","ATB"="red")
    mdsPlot(getBeta(mSetSq), numPositions = 1000, sampNames = samples$SampleName, sampGroups = samples$Sample_Group,
            pch = 19, pal = col1, legendPos = "bottomleft",
            main = "Beta MDS")
    
  }
  
  name<-paste0(mycell,"_Raw_mds_swanNorm_IDs.png")
  jpeg(name,width = 589,height = 426)
  #par(mfrow=c(1,1))
  MDSplot()
  dev.off()
  
  dim(mSetSq)

  ##-------------------- Filtering of the data
  detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
  
  ## remove any probes that have failed in one or more samples
  keep <- rowSums(detP < 0.01) == ncol(mSetSq)
  table(keep)

  mSetSqFlt <- mSetSq[keep,]
  dim(mSetSqFlt)
  ## [1] 863164     14
  
  ## if your data includes males and females, remove probes on the sex chromosomes
  keep <- !(featureNames(mSetSqFlt) %in% ann850k$Name[ann850k$chr %in%
                                                        c("chrX","chrY")])
  table(keep)
  ##
  # keep
  # FALSE   TRUE
  #  19176 843988 
  mSetSqFlt <- mSetSqFlt[keep,]
  dim(mSetSqFlt)
  ## [1] 843988     14
  
  ## remove probes with SNPs at CpG site
  gset <- mapToGenome(mSetSqFlt)
  mSetSwFlt <- dropLociWithSnps(gset)
  dim(mSetSwFlt)
  ##[1] 835693     14
  dim(samples)
  ##[1] 14 52
  
  ## exclude cross reactive probes
 
  nsp <- read.table("~/Linkoping University/Data/48639-non-specific-probes-Illumina450k.csv",sep=",",header=T)
  keep <- !(featureNames(mSetSwFlt) %in% nsp$TargetID )
  table(keep)
  ##

  mSetSwFlt <- mSetSwFlt[keep,]
  dim(mSetSwFlt)

  ##
  # Calculation of M values and beta values
  ## calculate M-values for statistical analysis
  mVals <- getM(mSetSqFlt)
  bVals <- getBeta(mSetSqFlt)
  colnames(bVals)<-samples$SampleName
  match(colnames(bVals),samples$SampleName)
  
  ## MDS cleaned
  
  name<-paste0(mycell,"_After_Quality_MDS.png")
  jpeg(name,width = 589,height = 426)
  #par(mfrow=c(1,1))
  mdsPlot(dat =bVals , numPositions = 1000, sampNames = samples$SampleName, sampGroups = samples$Sample_Group,
          pch = 19, pal = col1,
          legendPos ="topright",
          main = "Normalized beta matrix After QC")
  dev.off()
  
  ## Denisty plot for filtered bVals and mVals
  
  ##---------- save the R objects
   save(bVals,file="lnb_bVals_220816.RData")
  # save(mVals,file="MS_mVals_15032022.RData")
   save(samples,file="samples_info_LNB_220816.RData")
  
} ##end of minfi

############################################################
## 4. SVD analysis using ChAMP
## 5. Batch affect correction using Combat
## 6. Confounders and cov. identification
#############################################################

SVDplot.champ<-function(bVals,samples)
{
  
  library(ChAMP)
  library("minfi")
  mycell<-"LNB"
  col1=c("LNB"="red","non-LNB"="blue")
  
  
  setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")
  load("lnb_bVals_220816.RData")
  load("samples_info_LNB_220816.RData")
  #load("MS_bVals.corrected_15032022.RData")
  samples<-read.table("samples_info_LNB_220816.csv",header=T,sep=",")
  ##-------------------- SVD and batch effects
  names(samples)
  names(samples)[50:51]<-c("tickBite","Borreliosis")
  tmp_samples<-samples[,c(1:3,5,49:52)]
  tmp_samples<-tmp_samples[,-c(5,6,7,8)]
  colnames(tmp_samples)
  head(tmp_samples)
  colnames(samples)[5]<-"Age_cat"
  colnames(tmp_samples)[4]<-"Age"
  
  table(tmp_samples$Age_cat)
  match(tmp_samples$SampleName,colnames(bVals))
  
  ##-----------
  
  # Calculate first 10 SVD components (champ.SVD)
  source("calculate_svd.R")
  
  pdf(file = "svd_bVals.corrected.pdf", 
      width = 10, height = 8)
  calculate_svd(counts = bVals.corrected, 
                metadata = tmp_samples, 
                components = 10, 
                Rplot = TRUE)
  dev.off()
  ##-----------
  #champ.SVD(beta=bVals,pd=tmp_samples)
  bVals.corrected<-champ.runCombat(beta=bVals,pd=tmp_samples)

  champ.SVD(beta=bVals.corrected,pd=tmp_samples,resultsDir = "SVD_bVals.corrected_230224_")
  save(bVals.corrected,file="LNB_bVals.corrected_220816.RData") ##***
  save(samples,file="samples_info_LNB_220816.RData") ##***
  
  ## plot MDS of the batch corrected beta values
  
  name<-paste0(mycell,"_BatchCorrected_MDS_ID.png")
  jpeg(name,width = 589,height = 426)
  
  
  mdsPlot(dat =bVals.corrected, numPositions = 1000, sampNames = samples$SampleName, sampGroups = samples$Sample_Group,
          pch = 19, pal = col1,
          legendPos ="topright",
          main = "MDS of Norm. corrected beta matrix")
  dev.off()
  
  ##------------- check all clinical info
  library("ggplot2")
  setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")
  load("LNB_bVals.corrected_220816.RData") ##***
  load("samples_info_LNB_220816.RData") ##***
  colnames(samples)
  colnames(samples)[50:51]<-c("tickBite","Borreliosis")
  colnames(samples)[5]<-"Age_cat"
  colnames(samples)[8:11]<-c("Borrelia_IgG.serum","Borrelia_IgM.serum","Borrelia_IgG.CFS","Borrelia_IgM.CFS")
  colnames(samples)[c(23,25,27,29,31,33,35,37,39,41,43,45,47)]<-c("Headache","Fatigue","Fever","Neck.pain","Poor.apetite","Nausea","Weight.loss","Vertigo","Concentration.diff","Radiculitis","Myalgia.Arthralgia","Paresthesia","Facial.palsy")
  tmp_samples<-samples[,c(1:3,5,50,51,8:11,23,25,27,29,31,33,35,37,39,41,43,45,47,52,53)]
  tmp_samples[,c(2:25)] <- sapply(tmp_samples[,c(2:25)],as.character)
  sapply(tmp_samples, class)
  #tmp_samples2<-tmp_samples[,c(1,3,4,5,6,7)]
  library("ChAMP")
  champ.SVD(beta=bVals.corrected,pd=tmp_samples,resultsDir = "SVD_bVals.corrected_230309_")
  save(tmp_samples,file="tmp_samples_info_LNB_230309.RData") ##***
  
  
} ## end of batch correction and SVD

#####################################################################################################
## 6. Density plots, MDS plots, PCA plots to visualize the data (before and after QC/correction) ##
######################################################################################################

##----------------------  PCAPlot.merge
PCAPlot.merge <-function(bVals.corrected,tmp_samples,p1,p2)
{
  library(FactoMineR)
  library(factoextra)
  library(gridExtra)
  library(devtools)
  library(jcolors)
  library(ggpubr)
  library(ggplot2)
  library(pander)
  library(gridExtra)
  
  data1_top <- as.matrix(bVals.corrected)
  betaPCAt<-as.data.frame(t(data1_top))
  betaPCAt$sample<-tmp_samples$Sample_Group
  betaPCAfig <- betaPCAt[,1:length(betaPCAt)-1]
  dim(betaPCAfig)
  PCA <- PCA(betaPCAfig, graph = F,ncp=10) # beta PCA
  #c("ATB"="red","CT"="orange","HC"="mediumseagreen")
  betaPCAt$sample <- factor(betaPCAt$sample, levels = c("LNB","non_LNB"))
  #betaPCAt$sample <- factor(betaPCAt$sample, levels = c("1","2"))
  p <- fviz_pca_ind(PCA, repel = T,axes = c(p1,p2),
                    #geom.ind ="point",
                    geom.ind =c("point","text"),
                    pointsize = 5,
                    col.ind = betaPCAt$sample,
                    mean.point = FALSE,
                    palette = c("LNB"="red","non_LNB"="blue"),
                    #palette = c("1"="red","2"="blue"),
                    legend.title = "Group",
                    addEllipses = T,ellipse.level = 0.70)+
    scale_shape_manual(values=c(19,19))
  #+labs(title ="", x = "PC1", y = "PC2")
  
  
  plot1<-p+
    theme(text = element_text(family = "Arial", size = 28),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(size = 2))
  
  name<-paste0("PCA_BatchCorrected_LNB_230224.png")
  png(name,height = 10, width = 10, units = "in", res = 300 )
  grid.arrange(plot1, ncol=1, nrow =1)
  dev.off()
  
}

PCAPlot.merge(bVals.corrected = bVals.corrected,tmp_samples = samples,2,3)


###############################################################################

# 8.A Cell type adjustments and analysis
############### 1. Cell deconvolution: Estimate of cell counts using Houseman et al.

###############################################################################

run_houseman <- function(houseman){
  
  counts_hm<-FlowSorted.Blood.EPIC::estimateCellCounts2(rgSet = rgSet,compositeCellType = "Blood",
                                                        processMethod = "preprocessQuantile",
                                                        referencePlatform ="IlluminaHumanMethylationEPIC",
                                                        returnAll = TRUE,meanPlot = FALSE,
                                                        probeSelect = "IDOL",
                                                        IDOLOptimizedCpGs =IDOLOptimizedCpGs)
  
  counts_cellprop<-counts_hm$counts
  write.table(counts_cellprop,
              "HouseMan_counts_cellprop.csv",
              sep=",",quote=F)
  ## perform the stats on different cell types
  ## to get the signf. difference between groups
  
}

###############################################################################

# 8.B Cell type adjustments and analysis
############### 2. EPIDISH : est. cell prop. for Epi, Fib, Fat and IC.
## used the robust partial correction method

###############################################################################


Epidish_analysis<-function()
{
  #BiocManager::install("EpiDISH")
  library(EpiDISH)
  library(gridExtra)
  data(centEpiFibIC.m)
  
  Frac.m <- epidish(beta.m = bVals.corrected, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
  
  name<-paste0("LNB","_Cell_type_propotions.png")
  png(name,height = 10, width = 10, units = "in", res = 300 )
  boxplot(Frac.m)
  dev.off()
  
  write.table(Frac.m, "Epidish_LNB_Frac.m.csv",
              sep = ",",
              quote = F,
              row.names = T)
  
  ## Compare the cell types between the Groups
  Cell_type_Plot()
  
}

######################################################################

## 9.DMC idenfication using linear model (Limma package in R)

######################################################################

DMC_identification<-function()
{
  library(limma)
  library(minfi)
  # load("samples_info_MS_15032022.RData")
  # load("MS_bVals.corrected_15032022.RData")
  
  load("LNB_bVals.corrected_220816.RData") ##***
  load("samples_info_LNB_220816.RData") ##***
  
  ## removing the outlier identified based on the PCA. This will be added back in the heatmap later.
  bVals.corrected_wo = bVals.corrected[,-12]
  targets = samples[-12,]
  #<-samples
  match(colnames(bVals.corrected_wo),targets$SampleName)
  group <- factor(targets$Sample_Group)
  sex <- factor(targets$Sex)
  
  #design_2<-model.matrix(~0+group,data=targets)
  design_3<-model.matrix(~0+group+sex,data=targets)

  colnames(design_3) <- c(levels(group),"sex")

  contMatrix <- makeContrasts(LNB-non_LNB,levels=design_3)
  
  dim(design_2)
  
  fit<-lmFit(bVals.corrected_wo,design_3)
  colnames(fit$coefficients)
  rownames(contMatrix)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  summary(decideTests(fit2))
  
  ##--- Annotate the DMCs
  
  ann850kSub_combat <- ann850k[match(rownames(bVals.corrected_wo),ann850k$Name),
                               c(1:ncol(ann850k))]
  
  
  res1_CovSex<-topTable(fit2,coef = 1,number = Inf,genelist = ann850kSub_combat)
  
  
 
  save(res1_CovSex,file="res1_CovSex_LNB_nonLNB_limma.RData")
 
  load("res1_CovSex_LNB_nonLNB_limma.RData")

  
  ## No DMC passed the FDR cutoff
  ##----------- Top DMC identified based on the nominal p values
  topset3_sex<-subset(res1_CovSex, res1_CovSex$P.Value<0.05 & abs(res1_CovSex$logFC)>0.2)
  write.table(topset3_sex,"DMC_LNB_nominal_p0.05_MMD0.2_SEXCov.csv",sep = ",",quote=F, row.names = F)
  dim(topset3_sex)

  ## total DMCs : 428
  
  #Hypermethylated DMCs
  topset3_sex_hyper<-subset(res1_CovSex, res1_CovSex$P.Value<0.05 & res1_CovSex$logFC>0.2)
  dim(topset3_sex_hyper)
  ## 256
  
  #Hypomethylated DMCs
  topset3_sex_hypo<-subset(res1_CovSex, res1_CovSex$P.Value<0.05 & res1_CovSex$logFC < -0.2)
  dim(topset3_sex_hypo)
  ## 172
  CHeatmapPlot_220818(betaPCA = bVals.corrected,tp1 = topset1_sex,samples = samples,title1 = "Heatmap_full_pval_0.05_MMD0.3_230309_SEXCov")
  CHeatmapPlot_220818(betaPCA = bVals.corrected,tp1 = topset3_sex,samples = samples,title1 = "Heatmap_full_pval_0.05_MMD0.2_230309_SEXCov")
  
  ##------------ Super relaxed list
  #relaxed_sex<-subset(res1_CovSex, res1_CovSex$P.Value<0.05 & abs(res1_CovSex$logFC)>0.1)
  #dim(relaxed_sex)
  #write.table(relaxed_sex,"DMC_LNB_nominal_p0.05_MMD0.1_SEXCov_relaxed.csv",sep = ",",quote=F, row.names = F)
  
  
}

## Find DMR 23-05-17
DMR_identification<-function()
{

  library("ChAMP")
  library(limma)
  library(minfi)
  library(missMethyl)
  library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
  library(DMRcate)
  
  load("Normalization/bVals.corrected.RData") ##***
  load("Normalization/samples.RData") ##***
  
  colnames(samples)
  ## removing the outlier
  match(samples$Sample_Name,colnames(bVals.corrected))
  targets<-samples
  match(colnames(bVals.corrected),targets$SampleName)
  
  # Region level analysis 
  
  group <- factor(targets$Sample_Group)
  sex <- factor(targets$Sex)

  design<-model.matrix(~group+sex,data=targets)
  #colnames(design) <- c(levels(group),"sex")
  contMatrix 

  myAnnotation <- cpg.annotate(object = bVals.corrected, datatype = "array", what = "Beta", 
                               arraytype = c("EPIC"), 
                               analysis.type = "differential", design = design, 
                               fdr=0.05,
                               coef = 1)
 
  
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, betacutoff=0.2)
  
  results.ranges <- extractRanges(DMRs)
  results.ranges

  # Plot DMR
  cols <- c(1,2)[factor(samples$Sample_Group)]
  names(cols) <-factor(samples$Sample_Group)
  
  par(mfrow=c(1,1))
  # dmr_region selects which region to plot
  dmr_region <- 5
  DMR.plot(ranges=results.ranges, dmr=dmr_region, CpGs=bVals.corrected, phen.col=col1, 
           what="Beta", arraytype="EPIC", genome="hg38")

  ##--- Annotate the DMRs
  #Eight DMR found
  
}
######################################################################

## 10.Volcano, Heatmaps and Pathway Enrichment (ORA)/Network analysis

## Using the DMCs identified above, you can merge them or use them individually
## generate the heatmaps and volcano

######################################################################

## updated 220818
CHeatmapPlot_220818<-function(betaPCA,tp1,samples,title1){
  
  library(ComplexHeatmap)
  library(dendsort)
  library("RColorBrewer")
  library("viridis")
  top_betaPCA<-betaPCA[rownames(betaPCA) %in% tp1$Name,]
  dim(samples)
  dim(top_betaPCA)
  ## load data end
  match(colnames(top_betaPCA),samples$SampleName)
  data1_top <- as.matrix(top_betaPCA) ##*
  betaPCAt<-as.data.frame(t(data1_top))
  #betaPCAt$fatigue<-samples$Fatigue.1.yes..0.no
  betaPCAt$group<-samples$Sample_Group
  betaPCAt$sex<-samples$Sex
  #betaPCAt$fatigue<- factor(betaPCAt$fatigue, levels = c("0", "1"))
  betaPCAt$group <- factor(betaPCAt$group, levels = c("LNB", "non_LNB"))
  
  betaPCAt$sex <- factor(betaPCAt$sex, levels = c("1", "2"))
  
  col = list( 
    #fatigue= c("0"="grey","1"="black"),
              group= c("LNB"="maroon","non_LNB"="royalblue")
              ,sex= c("2"="pink","1"="lightskyblue1")
  )
  ha <- HeatmapAnnotation(
    group= betaPCAt$group,
    #fatigue= betaPCAt$fatigue,
    sex=betaPCAt$sex,
    col = col,
    
    annotation_label = gt_render(c("Group","Sex"),
                                 gp = gpar(col = "black", fontsize=22)),
    annotation_name_side = "right",
    simple_anno_size =unit(1, "cm"),
    
    #Sample_Group=anno_text(gp = gpar(fontsize = 22)),
    annotation_legend_param = list(
      group = list(
        title = "Group",
        ncol = 1, by_row = FALSE,
        title_gp = gpar(col = "black", fontsize = 22),
        grid_height = unit(1, "cm"), 
        grid_width = unit(1, "cm"),
        labels_gp = gpar(col = "black", fontsize = 22)),
      # fatigue = list(
      #   title = "Fatigue",
      #   title_gp = gpar(col = "black", fontsize = 22),
      #   grid_height = unit(1, "cm"), 
      #   grid_width = unit(1, "cm"),
      #   labels_gp = gpar(col = "black", fontsize = 22)),
      sex = list(
        title = "Sex",
        title_gp = gpar(col = "black", fontsize = 22),
        grid_height = unit(1, "cm"), 
        grid_width = unit(1, "cm"),
        labels_gp = gpar(col = "black", fontsize = 22))
      
    )
  )
  
  library(viridis)
  row_dend = dendsort(hclust(dist(data1_top)))
  col_dend = dendsort(hclust(dist(t(data1_top))))
  
  # Heatmap(data1_top,top_annotation = ha, name = "Beta value", cluster_rows = row_dend, cluster_columns = col_dend,
  #        column_title = "",row_dend_reorder = TRUE,show_row_names=FALSE,show_column_names =FALSE,col = viridis(100),show_row_dend = FALSE)
  ## to be checked ..
  ht=Heatmap(data1_top,
             top_annotation = ha,
             name = "Beta value",
             cluster_rows = row_dend,
             #cluster_row_slices = FALSE,
             cluster_column_slices = FALSE,
             cluster_columns = col_dend,
             #col=viridis(100),
             heatmap_legend_param = list(
               title = "Beta value", at = c(0, 0.5,1),
               #direction = "horizontal",
               #column_gap = unit(6, "cm"),
               #x = unit(1, "npc"), y = unit(1, "npc"),
               #just = c("right", "top"),
               grid_width = unit(1, "cm"),
               legend_height = unit(5, "cm"),
               labels_gp = gpar(col = "black", fontsize = 21),
               title_gp = gpar(col = "black", fontsize = 21)
               
             ),
             #width = unit(8, "cm"), height = unit(8, "cm"),
             column_title = "",
             row_dend_reorder = TRUE,
             column_dend_reorder = TRUE,
             
             show_row_names=FALSE,
             show_column_names =TRUE,
             ##--- *** --- ##
             
             #show_column_dend = TRUE,
             show_row_dend = TRUE)
  ## adjust global parameters in heatmap:Chapter 4.13
  ht_opt(legend_border = "black",legend_labels_gp=gpar(lwd=3),
         legend_gap=unit(7,"mm"),
         #ANNOTATION_LEGEND_PADDING=unit(7,"mm"),
         HEATMAP_LEGEND_PADDING=unit(3,"mm"))
  ##
  #name<-paste0(title1,"_CHeatmap_230309.pdf")
  #png(name,height = 10, width = 10, units = "in", res = 300)
  pdf(file = paste0(title1,"_CHeatmap_230309.pdf"))
  
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}
CHeatmapPlot_220818(betaPCA = bVals.corrected,tp1 = topset1_sex,samples = samples,title1 = "Heatmap_full_pval_0.05_MMD0.3_230309_SEXCov")
CHeatmapPlot_220818(betaPCA = bVals.corrected,tp1 = topset3_sex,samples = samples,title1 = "Heatmap_full_pval_0.05_MMD0.2_230309_SEXCov")

######################################################################

## 11.Estimated genomic inflation and bais using BACON package in R.

######################################################################

##----------------
run_bacon<-function(){
  library("bacon")
  bc1<-bacon(res1$t)  ##
  estimates(bc1)
  inflation(bc1)
  bias(bc1)
  bc<-bc1
  plot(bc1, type="qq")
}

##--------------- end of pipeline


### Section 2 for plotting
MDS_plots<-function()
{
  library("minfi")
  LNB_mdsPlot_v5_220816 <- function(dat, numPositions = 1000, sampNames = NULL,
                                    sampGroups = NULL, Vacc = NULL, xlim, ylim, pch = 1,
                                    pal = brewer.pal(8, "Dark2"),
                                    #pal = "black",
                                    legendPos = "bottomleft",
                                    legendNCol, main = NULL, sub=NULL, bty="L",
                                    bg_col=brewer.pal(8, "Dark2"),lwd=lwd.igra) {
    # Check inputs
    if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
      b <- getBeta(dat)
    } else if (is(dat, "matrix")) {
      b <- dat
    } else {
      stop("dat must be an 'MethylSet', 'RGChannelSet', or 'matrix'.")
    }
    if (is.null(main)) {
      main <- sprintf(
        "Beta MDS\n%d most variable positions",
        numPositions)
    }
    
    o <- order(rowVars(b), decreasing = TRUE)[seq_len(numPositions)]
    d <- dist(t(b[o, ]))
    fit <- cmdscale(d)
    if (missing(xlim)) xlim <- range(fit[, 1]) * 1.2
    if (missing(ylim)) ylim <- range(fit[, 2]) * 1.2
    if (is.null(sampGroups)) sampGroups <- rep(1, numPositions)
    sampGroups <- as.factor(sampGroups)
    col <- pal[Vacc]
    Vacc <- as.factor(Vacc)
    bg<-bg_col[sampGroups]
    lwd<-lwd[Vacc]
    if (is.null(sampNames)) {
      par(mar = c(5, 4.7, 4, 10),
          xpd = TRUE,lwd=3)
      plot(
        x = fit[, 1],
        y = fit[, 2],
        col = col,
        bg=bg,
        pch = pch,
        xlim = xlim,
        ylim = ylim,
        xlab = "Dimension 1",
        ylab = "Dimension 2",
        cex=2.5,
        las=1,
        type="p",
        #lwd=5,
        #lty=1,
        cex.main=2,
        cex.axis=2,
        cex.lab=1.7,
        sub=sub,
        bty=bty,
        lwd=lwd,
        main = main)
      # axis(side=1,lwd = 3,cex.axis=2)
      #axis(side=2,lwd = 3, las=1, cex.axis=2)
      #box(lwd=3)
      
    } else {
      par(mar = c(5, 4, 4, 8),
          xpd = TRUE)
      plot(
        x = 0,
        y = 0,
        type = "L",
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = "",
        main = main)
      #axis(side=1,lwd = 3,cex.axis=1.5)
      #axis(side=2,lwd = 3,cex.axis=1.5)
      #box(lwd=3)
      text(x = fit[, 1], y = fit[, 2], sampNames, col = col)
      
      
    }
    
  }
  
  
  ## ----------- corrected bVals 
  load("LNB_bVals.corrected_220816.RData") ##***
  load("samples_info_LNB_220816.RData") ##***
  
  col.vacc=c("LNB"="black","non_LNB"="black")
  lwd.vacc=c("LNB"="1","non_LNB"="1")
  group= c("LNB"="maroon","non_LNB"="royalblue")
  ticks<-samples$Have.noticed.a.tick.bite.1.yes..0.no
  
  par(lwd=3)
  #png("mdsPlot_bVals_corrected_v1.png",height = 150, width = 170, unit = "mm", res = 300)
  pdf(file = "mdsPlot_bVals_corrected_v1.pdf",width = 7,height = 7)
  
  
  LNB_mdsPlot_v5_220816(dat = bVals.corrected, numPositions = 1000, 
                        sampNames = NULL, 
                        sampGroups = samples$Sample_Group,
                        Vacc = samples$Sample_Group,
                        pch = 21,
                        pal = col.vacc,
                        #pal="transparent",
                        bg_col=group,
                        legendPos ="topright",lwd=lwd.vacc,
                        main = "", sub = "")
  legend("topright",title="Group",legend=c("non-LNB","LNB"),col=c("royalblue","maroon"),pch=19,bty="n",cex=1.5, inset = c(- 0.50, 0.20))
  #legend("topright",title="Ticks",legend=levels(factor(ticks)),col=col.vacc[seq_len(length(levels(factor(ticks))))],pch=21,bty="n",cex=2, inset = c(- 0.42, 0),pt.lwd =lwd.vacc)
  
  dev.off()
  
}

#########################---------------- boxplot for the cell types
Cell_type_Plot<-function()
{
  # library
  library(ggpubr)
  library(ggplot2)
  library(pander)
  library(gridExtra)
  #install.packages("tidyr")
  library(tidyr)
  
  load("samples_info_LNB_220816.RData")
  Frac.m<-read.table("Epidish_LNB_Frac.m.csv",
                     sep = ",",
                     header=T)
  counts<-Frac.m[,-7]
  
  Cell_type<-rep(c("B","NK","CD4T","CD8T","Mono","Neutro"),each=14)
  Group=samples$Sample_Group
  match(rownames(counts),samples$SampleName)
  counts<-data.frame(counts)
  data<-counts %>% gather(cellType, cellProp)
  data1=data.frame(Cell_type,Group,data[,2])
  Prop<-data[,2]
  
  Group<-factor(samples$Sample_Group, levels = c("LNB","non-LNB"))
  # grouped boxplot
  p1<-ggplot(data1, aes(x=Cell_type, y=Prop, fill=Group,col=Group)) + 
    #geom_boxplot()
    geom_boxplot(alpha=0.9,outlier.shape = 20,lwd=2)+
    # #geom_text(check_overlap = TRUE,size =8,
    # #         position=position_jitter(width=0.15))+
    # geom_jitter(position = position_jitter(0.1), size =6)+
    #theme_pubr(base_size = 22, margin = FALSE,legend="right")+
     theme(plot.title = element_text(color="black", size=26),
           #axis.text = element_text(color="black", size=26),
            
           panel.background = element_blank(),
           axis.line = element_line(colour = "black",size = 1.5),
           axis.text.x = element_text(size=26,color="black"),
           axis.text.y = element_text(size=26,color="black"),
           plot.margin = margin(0.8,0.8,0.1,0.1, "cm"),
           plot.tag.position = "top",
           
           legend.key.size = unit(2, 'cm'),
           legend.background = element_blank(),
           legend.text = element_text(size=26),
           legend.title = element_text(size=26))+
    scale_fill_manual(values = c("maroon","royalblue"))+
    scale_color_manual(values = c("maroon","royalblue"))+
    
    #stat_compare_means(method = "t.test", label = "p.signif", label.y = 0.8)+
    #scale_shape_manual(values=pch.cohort)+
    #geom_point(aes(colour=factor(betaPCAt$sample)), shape=21)+
    ##xlab()+
    #labs(y=beta~"value", x=NULL)
    #labs(y="Proportion", x="cell type")
  png("CellType_proportion_LNB_1.png",height = 10, width = 10, units = "in", res = 300)
  #png("CellType_proportion_LNB.png",height = 150, width = 170, unit = "mm", res = 300)
  #par(lwd=3)
  p1
dev.off()   
  
}

## 220902

Fig_Dotplot<-function(data=data)
{
  # Install and load the package (already installed in the Docker version)
  #BiocManager::install("clusterProfiler")
  #library(clusterProfiler)
  library(ggpubr)
  library(ggplot2)
  library(pander)
  library(gridExtra)
  library(scales)
  library(rasterVis)
  library(dplyr)
  library(tidyverse)
  #library(httr)
  head(data)
  data<-data.frame(data)
  dim(data)
  #data_sort= order(data, by= "Gene.Ratio")
  colnames(data) <- c("PID", "Pathway", "link","Gene.Count", "Overlap", "Expected", "Gene.Ratio", "p.value","FDR","OL_ID","UserID")
  dim(data)
  #data <- data[1:10,]
  data
  
  dim(data)
  #data_sort<-data[order(as.numeric(data$Gene.Ratio)),]
  # Dot plot
  attach(data)
  p <- ggplot() +
    geom_point(data=data,
               aes(x= Gene.Ratio, y= reorder(Pathway,Gene.Ratio), size= Gene.Count, fill= p.value),
               shape=21, alpha=0.8, stroke=1) +
    xlab("Enrichment Ratio")+
    ylab("")+
    
    #xlim(0,3)+ ##Fig 2B, 
    xlim(1,16)+ ##Fig 2A,
    
    scale_size_continuous(range = c(4, 8), breaks=pretty_breaks(4)) +
    scale_fill_distiller(direction = 1, palette="RdBu", breaks=pretty_breaks(4)) +
    guides(fill = guide_colorbar(order = 2), size = guide_legend(order = 1))
  
  detach(data)
  # ***## Main pathway figs 4A and 5A
  plot1<-p+theme_pubr(base_size = 13,legend="right")
  #+theme(plot.title = element_text(color="black", size=24, face="bold"),
  #axis.text = element_text(color="black", size=24, face="bold"),
  
  #panel.background = element_blank(),
  #axis.line = element_line(colour = "black",size = 1.5),
  #plot.margin = margin(0.8,0.8,0.1,0.1, "cm"),
  #plot.tag.position = "top",
  #legend.key.size = unit(1.5, 'cm'),
  #legend.title = element_text(size=22, face="bold")
  #  )
  
  #plot1<-p+labs_pubr() ## Suppl. Fig S5A, S5B
  plot1
  
}

##MMD 0.3 and nominal pval<0.05

# 1) LNB vs non-LNB
setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")
data <- read.delim("Pathways/MMD0.3_p0.05/enrichment_results_LNB_vs_non_ReactomePathways.txt", sep="\t", header = T)
png("Pathway_Reactome_LNB_MMD0.3_p0.05_220902.png",height = 125, width = 170, unit = "mm", res = 300)
LNB_pathways<-Fig_Dotplot(data = data)
png("Pathway_Reactome_LNB_MMD0.3_p0.05_220902.png",height = 10, width = 10, unit = "cm", res = 300)
LNB_pathways
dev.off()



MDS_plots_SupplementaryFigure<-function()
{
  library("minfi")
  
  source("LNB_mdsPlot_v5_220816.R") ## normal MDS
  source("Scripts/LNB_mdsPlot_v5_230310_studygroup.R")
  
  ## ----------- corrected bVals 
  load("LNB_bVals.corrected_220816.RData") ##***
  load("tmp_samples_info_LNB_230309.RData") ##****
  #load("samples_info_LNB_220816.RData")
  
  name<-paste0("Supplementary_plots_MDS/","ALL_supplementary_MDS_230310_withShapes.pdf")
  pdf(file = name,width = 7,height = 7)
  
  par(lwd=3)
  tmp_samples$Sample_Group
  group<-factor(tmp_samples$Sample_Group)
  pch.group=c("LNB"=24,"non_LNB"=22)
  col.vacc=c("LNB"="black","non_LNB"="black")
  lwd.vacc=c("LNB"="1","non_LNB"="1")
  
  ##---------- MDS with shapes
  for(i in 3:length(names(tmp_samples)))
  {
    #name<-paste0("Supplementary_plots_MDS/",names(tmp_samples)[i],"_supplementary_230309.pdf")
    #pdf(file = name,width = 7,height = 7)
    temp<-tmp_samples[,i]
    group.color <- as.integer(factor(temp))
    col<-levels(factor(group.color))
    match(colnames(bVals.corrected),tmp_samples$SampleName)
    LNB_mdsPlot_v5_230310_studygroup(dat = bVals.corrected, numPositions = 1000, 
                          sampNames = NULL, 
                          sampGroups = temp,IGRA = temp,cohort = group, 
                          pch = pch.group,
                          pal = col.vacc,
                          #pal="transparent",
                          bg_col=col,
                          legendPos ="top",lwd=lwd.vacc,
                          main = "", sub = "", bty = "L")
    legend("topright",title=names(tmp_samples)[i],legend=levels(factor(temp)),col=col,pch=19,bty="n",cex=1.5, inset = c(- 0.50, 0.30))
    #legend("topright",title="IGRA",legend=levels(factor(samples_wo_CC19$Vaccinated)),col=col.igra[seq_len(length(levels(factor(samples_wo_CC19$Vaccinated))))],pch=21,bty="n",cex=1.5, inset = c(- 0.4, 0),pt.lwd =lwd.igra)
    legend("topright",title="Group",legend=levels(factor(tmp_samples$Sample_Group)),pch=c(24,22),bty="n",cex=1.5, inset = c(- 0.50, 0))
    
    
  } 
  dev.off()
  
  
  ##--------- MDS without shapes
  for(i in 3:length(names(tmp_samples)))
        {
          #name<-paste0("Supplementary_plots_MDS/",names(tmp_samples)[i],"_supplementary_230309.pdf")
          #pdf(file = name,width = 7,height = 7)
            temp<-tmp_samples[,i]
            group.color <- as.integer(factor(temp))
            col<-levels(factor(group.color))
            match(colnames(bVals.corrected),tmp_samples$SampleName)
            LNB_mdsPlot_v5_220816(dat = bVals.corrected, numPositions = 1000, 
                                  sampNames = NULL, 
                                  sampGroups = temp,
                                  Vacc = temp,
                                  pch = 21,
                                  pal = col.vacc,
                                  #pal="transparent",
                                  bg_col=col,
                                  legendPos ="topright",lwd=lwd.vacc,
                                  main = "", sub = "")
            legend("topright",title=names(tmp_samples)[i],legend=levels(factor(temp)),col=col,pch=19,bty="n",cex=1.5, inset = c(- 0.50, 0.20))
           
            #legend("topright",title="Ticks",legend=levels(factor(ticks)),col=group.color[levels(factor(group.color))],pch=21,bty="n",cex=2, inset = c(- 0.42, 0),pt.lwd =lwd.vacc)
          #dev.off() 
            
        }
  dev.off()
 
  
}

load("res1_CovSex_LNB_nonLNB_limma.RData")

setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")


plot_volcano_hyper_hypo<-function(res){
  res<-res1_CovSex
  library(EnhancedVolcano)
  keyvals.colour <- ifelse(
    res$P.Value<=0.05 & as.numeric(res$logFC) <= -0.2 , 'blue',
    ifelse(res$P.Value<=0.05 & as.numeric(res$logFC) >= 0.2, 'red',
           'grey22'))
  keyvals.colour[is.na(keyvals.colour)] <- 'grey22'
  names(keyvals.colour)[keyvals.colour == 'red'] <- 'Hyper'
  
  names(keyvals.colour)[keyvals.colour == 'grey22'] <- 'Not Sign.'
  names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Hypo'
  
  #jpeg("VolcanoPlot_15022022_yb_res1_PC19_vs_Con.jpeg",width = 700,height = 600)
  p1<-EnhancedVolcano(res,
                      lab = '',
                      caption = '',
                      x = 'logFC',
                      y = 'P.Value',
                      #caption = 'cutoff: p-value 0.01; MMD 0.2',
                      selectLab = rownames(res)[which(names(keyvals.colour) %in% c('Hyper', 'Hypo'))],
                      xlab = 'Mean methylation difference (MMD)',
                      ylab = '-log10(P value)',
                      title = '',
                      subtitle = '',
                      pCutoff = 0.05,
                      FCcutoff = 0.2,
                      pointSize = 4,
                      xlim = c(-0.6,0.6),
                      ylim= c(0, -log10(10e-13)),
                      shape=19,
                      #labSize = 22.0,
                      #legendLabels=c('Not sig.','MMD','Adj.p-value','Adj.p-value & MMD'),
                      #shapeCustom = keyvals.shape,
                      colCustom = keyvals.colour,
                      colAlpha = 1,
                      legendPosition = 'right',
                      legendLabSize = 22,
                      
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      colConnectors = 'grey50',
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0.7
                      #border = 'full',
                      #borderWidth = 1.0,
                      #borderColour = 'black'
  )
  #png("Volcano_MMD_0.2_PVal_0.05_230310.png",height = 160, width = 190, unit = "mm", res = 300)
  svg("Volcano_MMD_0.2_PVal_0.05_230310.svg")
  #pdf(file="Volcano_MMD_0.2_PVal_0.05_230309.pdf")
  p1 + theme(
    text = element_text(colour = "black"),
    axis.text = element_text(size = 22, face = "plain",colour = "black")
    ,legend.text = element_text(size = 22)
    ,axis.title = element_text(size=22)
    ,axis.line = element_line(size = 1.25)
    #,legend.text.align = 1
  )
  
  ## OR # p1+theme_pubr(base_size = 22,legend = "right")
  dev.off()
  
}
plot_volcano_hyper_hypo(res = res1_CovSexAge)

##------------ Date 230223 Boc plots for LNB, HLA-DR cpG sites
Box_plots_HLA<-function()
{
  library("minfi")
  
  #source("LNB_mdsPlot_v5_220816.R") ## normal MDS
  #source("Scripts/LNB_mdsPlot_v5_230310_studygroup.R")
  
  ## ----------- corrected bVals 
  load("LNB_bVals.corrected_220816.RData") ##***
  load("tmp_samples_info_LNB_230309.RData") ##****
  #load("samples_info_LNB_220816.RData")
  
  tmp_samples$Sample_Name<-tmp_samples$SampleName
 setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")

 Top_Boxplots<-function(betaPCA,samples)
 {
   library(ggpubr)
   library(ggplot2)
   library(pander)
   library(gridExtra)
  
     betaPCA<-bVals.corrected
     samples<-tmp_samples
     
     match(colnames(betaPCA),samples$Sample_Name)
     
     HLA_Gene_list<-read.table("HLA_DMC.txt",sep="\t",header=T)
     
     CpG_listz<-HLA_Gene_list$CpG
     Gene_listz<-HLA_Gene_list$GeneName
     
     betaPCA_top<-betaPCA[rownames(betaPCA) %in% HLA_Gene_list$CpG,]
     
     x=HLA_Gene_list$CpG
     
     
     data1_top <- as.matrix( betaPCA_top)
     betaPCAt<-as.data.frame(t(data1_top))
     betaPCAt$sample<-samples$Sample_Group
     betaPCAt$sample<-factor(betaPCAt$sample, levels = c("LNB","non_LNB"))
     
     plot_list = list()
     
     
     ##---------- voilin grid
    
     #file = paste0(main_dir,output_dir, "BOXPLOT_COPE_230316_Top8.pdf")
     pdf(file = "Boxplot_LNB_DMC_in_HLA_230323_names.pdf")
     
     for (i in 1:length(x)) {
       
       dat<-betaPCAt[,c(x[i],"sample")]
       gene<- Gene_listz[i]
       col1=c("LNB"="maroon","non_LNB"="royalblue")
       group<-betaPCAt$sample
       ##print(dat)
       IDs<-rownames(betaPCAt)
       colnames(dat)<-c("beta","sample")
       ymin<-min(dat[,1])
       ymax<-max(dat[,1])
       p<-ggplot(dat, aes(x=sample, y=beta,labCol="",fill=group,col=group,label = rownames(betaPCAt))) +
         geom_boxplot(alpha=0.5)+
         geom_jitter(position = position_jitter(0.1), size =6)+
         theme_pubr(base_size = 18, margin = FALSE,legend="right")+
         theme(plot.title = element_text(color="black", size=24, face="bold"),
               axis.text = element_text(color="black", size=24, face="bold"),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black",size = 1.5),
               plot.margin = margin(0.8,0.8,0.1,0.1, "cm"),
               plot.tag.position = "top",
               legend.key.size = unit(2, 'cm'),
               legend.title = element_text(size=24, face="bold"))+
         ##c("royalblue","#7570B3","#E7298A")
         scale_fill_manual(values = c("LNB"="maroon","non_LNB"="royalblue"))+
         scale_color_manual(values = c("LNB"="maroon","non_LNB"="royalblue"))+
         
         #scale_fill_manual(values = c("royalblue","red"))+
         #scale_color_manual(values = c("royalblue","red"))+
         stat_compare_means(method = "wilcox.test", label = "p.signif", label.y = 0.8)+
         #ggtitle(x[i])
       ggtitle(paste0(gene,"_",x[i]))
       #scale_shape_manual(values=pch.cohort)+
       #geom_point(aes(colour=factor(betaPCAt$sample)), shape=21)+
       #xlab()+
       #labs(y=beta~"value", x=NULL)
       labs(y=NULL, x=NULL)
       
       plot_list[[i]] = p+ geom_text()
       plot(p)
       #plot(p+ geom_text(hjust=0, vjust=0))
     }
     dev.off()
     ## ---------------------- plot voilin
     # jpeg("~/PredictMe/PatentCpG_in_CD3_PACS_TB_Controls_list_39.jpeg",width = 1000,height = 1000)
     # 
     #  ggarrange(plotlist = plot_list,
     #           #labels = "",
     #           font.label = list(size=16,face="bold.italic"),
     #           hjust = -1.3,
     #           vjust = 1.9,
     #           ncol=5,nrow =8
     #           ,common.legend = TRUE
     #           ,legend = "none"
     # )
     # #+theme(plot.margin = margin(0.1,4,0.1,0.1, "cm"))
     # 
     # dev.off()
     
   
 }
 
 
 
}


##------Date 230517 GEO submission of processed datasets

##----------- 230424 Data uploaded in the GEO
GEO_submission<-function()
{
  setwd("/Users/Shumaila/OneDrive - Link?pings universitet/LNB_Anna_2022/Analysis_220816/")
  
  ##-------- used the dataset with Public IDS
  load("LNB_bVals.corrected_220816.RData") ##***
  load("samples_info_LNB_220816.RData") ##***
  
  match(samples$SampleName,colnames(bVals.corrected))
  samples$publicID<-c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14")
  ## S12 is an outlier corresponding to F32a.
  
  colnames(bVals.corrected) <- samples$publicID
  bVals_df <- cbind(as.data.frame(rownames(bVals.corrected)), as.data.frame(bVals.corrected))
  colnames(bVals_df)[1] <- "Probe_ID"
  rownames(bVals_df) <- NULL
  write.table(bVals_df,"~/Linkoping University/GEO_submissions/LNB_GEO/Preparation_GEO/LNB_beta_matrix/BetaVals_LNB_850K.txt",sep="\t",row.names = F,quote = F)
  
  save(samples,file="samples_info_LNB_GEO_Submission_publicIDs_230517.RData")
}

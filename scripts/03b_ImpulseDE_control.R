# Aim is to run ImpulseDE2 on control timecourse to compare to gradient tool results

# Philippa Borrill 
# 29.10.2018
# ImpulseDE is currently a preprint  http://dx.doi.org/10.1101/113548 
# but it has a tutorial which I will follow https://bioconductor.org/packages/release/bioc/vignettes/ImpulseDE2/inst/doc/ImpulseDE2_Tutorial.html 

# install ImpulseDE2
source("https://bioconductor.org/biocLite.R")
biocLite("ImpulseDE2")

library(ImpulseDE2)

out_dir <- "Y:\\PB_AFLF\\RefSeqv1.1_masked\\control_timecourse\\impulseDE2\\"

# this is the example from the tutorial using a simulated dataset
lsSimulatedData <- simulateDataSetImpulseDE2(
  vecTimePointsA   = rep(seq(1,8),3),
  vecTimePointsB   = NULL,
  vecBatchesA      = NULL,
  vecBatchesB      = NULL,
  scaNConst        = 30,
  scaNImp          = 10,
  scaNLin          = 10,
  scaNSig          = 10,
  scaMuBatchEffect = NULL,
  scaSDBatchEffect = NULL,
  dirOutSimulation = NULL)

# the count data is in the form:
head(lsSimulatedData$matObservedCounts)

# the annotation is in the form:
head(lsSimulatedData$dfAnnotation)

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = lsSimulatedData$matObservedCounts, 
  dfAnnotation    = lsSimulatedData$dfAnnotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

head(objectImpulseDE2$dfImpulseDE2Results)


### now want to do it with my control dataset
# should use counts as input
# only want to use genes >0.5 tpm

## do it for control first
count.data <- read.csv(file="Y:\\PB_AFLF\\RefSeqv1.1_masked\\expression_per_gene\\control_timecourse_count.tsv", sep = "\t")
head(count.data)
dim(count.data)

counts <- count.data
head(counts)
print("dimensions of count data before filtering")
dim(counts)

# now select only genes with > 0.5 tpm in at least 1 timepoint FLB
tpmData_FLB <- read.csv(file="Y:\\PB_AFLF\\RefSeqv1.1_masked\\expression_per_gene\\control_timecourse_tpm.tsv", header=T, sep = "\t")
head(tpmData_FLB)
colnames(tpmData_FLB)

# just use FLB data
tpmData_FLB <- tpmData_FLB[,1:30]
head(tpmData_FLB)

# average per timepoint
tpmData_FLB$T3 <- (tpmData_FLB[,1] + tpmData_FLB[,2] + tpmData_FLB[,3]) / 3
tpmData_FLB$T7 <- (tpmData_FLB[,4] + tpmData_FLB[,5] + tpmData_FLB[,6]) / 3
tpmData_FLB$T10 <- (tpmData_FLB[,7] + tpmData_FLB[,8] + tpmData_FLB[,9]) / 3
tpmData_FLB$T13 <- (tpmData_FLB[,10] + tpmData_FLB[,11] + tpmData_FLB[,12]) / 3
tpmData_FLB$T15 <- (tpmData_FLB[,13] + tpmData_FLB[,14] + tpmData_FLB[,15]) / 3
tpmData_FLB$T17 <- (tpmData_FLB[,16] + tpmData_FLB[,17] + tpmData_FLB[,18]) / 3
tpmData_FLB$T19 <- (tpmData_FLB[,19] + tpmData_FLB[,20] + tpmData_FLB[,21]) / 3
tpmData_FLB$T21 <- (tpmData_FLB[,22] + tpmData_FLB[,23] + tpmData_FLB[,24]) / 3
tpmData_FLB$T23 <- (tpmData_FLB[,25] + tpmData_FLB[,26] + tpmData_FLB[,27]) / 3
tpmData_FLB$T26 <- (tpmData_FLB[,28] + tpmData_FLB[,29] + tpmData_FLB[,30]) / 3

# just keep the average per timepoint
#head(tpmData_FLB)
colnames(tpmData_FLB)[31:40]
tpmData_FLB_av <- tpmData_FLB[,31:40]

tpmData_FLB_av$maxtpm <- apply(tpmData_FLB_av[,1:10],1,max)
head(tpmData_FLB_av)

# clean up workspace to remove unnecessary dataframes
rm(tpmData_FLB)

# merge together counts_conf  and tpmData_FLB_av
counts_max_FLB <- merge(counts, tpmData_FLB_av, by.x = 0, by.y = 0)
head(counts_max_FLB)
dim(counts_max_FLB)
# select only rows with a maxtpm >0.5
counts_FLB <- counts_max_FLB[which(counts_max_FLB$maxtpm>0.5),]


## remove LC genes
head(counts_FLB)
counts_FLB <- counts_FLB[!grepl("LC",counts_FLB$Row.names), ]
head(counts_FLB)
dim(counts_FLB)

# make rownames correct
rownames(counts_FLB) <- counts_FLB[,1]
counts_FLB <- counts_FLB[,-1]
head(counts_FLB)
#head(row.names(counts_FLB))

# remove tpm and maxtpm columns
counts_FLB <- counts_FLB[,1:30]
head(counts_FLB)
dim(counts_FLB)

# re-organise
counts_FLB_ordered <- counts_FLB[,c(1,4,7,10,13,16,19,22,25,28,2,5,8,11,14,17,20,23,26,29,3,6,9,12,15,18,21,24,27,30)]
head(counts_FLB_ordered)

### now make annotation dataframe with sample details. Need 4 columns "Sample", "Condition", "Time" and "Batch")
annotationdf <- data.frame("Sample" = colnames(counts_FLB_ordered), 
                           "Condition" = rep("case",30), 
                           "Time" = rep(c(3,7,10,13,15,17,19,21,23,26),3),
                           "Batch" = rep("B_NULL",30))
head(annotationdf)
annotationdf

# now let's run ImpulseDE2 on just 5 genes as a test

test_counts_FLB_ordered <- as.matrix(round(counts_FLB_ordered[1:5,]))
test_counts_FLB_ordered

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = test_counts_FLB_ordered, 
  dfAnnotation    = annotationdf,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

head(objectImpulseDE2$dfImpulseDE2Results)

# let's plot the gene expression
library(ggplot2)
lsgplotsGenes <- plotGenes(
  vecGeneIDs       = NULL,
  scaNTopIDs       = 3,
  objectImpulseDE2 = objectImpulseDE2,
  boolCaseCtrl     = FALSE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL,
  boolSimplePlot=FALSE)
print(lsgplotsGenes[[1]])


### now run for real on full set of ~52,905 genes
head(counts_FLB_ordered)
dim(counts_FLB_ordered)

counts_FLB_ordered <- as.matrix(round(counts_FLB_ordered))
head(counts_FLB_ordered)

real_objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts_FLB_ordered, 
  dfAnnotation    = annotationdf,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

head(real_objectImpulseDE2$dfImpulseDE2Results)

# save output
write.csv(file=paste0(out_dir,"control_ImpulseDE2_results.csv"), real_objectImpulseDE2$dfImpulseDE2Results, row.names = F)


result_df <- real_objectImpulseDE2$dfImpulseDE2Results
head(result_df)

head(result_df[result_df$padj <0.05,])
nrow(result_df[result_df$padj <0.05,])

nrow(result_df[result_df$padj <0.01,])

nrow(result_df[result_df$padj <0.001,])

## for some reason making plots doesn't show the dots of the normalised counts ###
library(ggplot2)
lsgplotsGenes <- plotGenes(
  vecGeneIDs       = NULL,
  scaNTopIDs       = 5,
  objectImpulseDE2 = real_objectImpulseDE2,
  boolCaseCtrl     = FALSE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL,
  boolSimplePlot=FALSE)
  
print(lsgplotsGenes[[1]]) 


# let's compare ImpulseDE2 to gradient tool results

##### now let's compare ImpulseDE2 to gradient tool results using gradient tool tpm normalised (my preferred one) ####
# get gradient tool results
gradienttool_res <- read.csv(file="Y:\\PB_AFLF\\RefSeqv1.1_masked\\control_timecourse\\gradient_tool\\control_gradient_tool_tpm_normalised\\sorted_up_down_data_HC_only.csv")
head(gradienttool_res)
dim(gradienttool_res)

gradienttool_res[gradienttool_res$X =="TraesCS4B02G311100",]

# now select only DE genes from gradient tool
gradienttool_res_DE <- gradienttool_res[gradienttool_res$pattern != "0_0_0_0_0_0_0_0_0_0",]
head(gradienttool_res_DE)
dim(gradienttool_res_DE)


# now select only DE genes from ImpulseDE2 padj0.001
impulse_DE <- result_df[result_df$padj <0.001,c(1,3)]
head(impulse_DE)
dim(impulse_DE)

# merge together and find overlap
merged_res <- merge(gradienttool_res_DE, impulse_DE, by.x="X", by.y="Gene", all = T)
head(merged_res)

dim(merged_res)
sum(complete.cases(merged_res))


write.csv(file=paste0(out_dir,"merged_impuseDE2_padj0.001_gradient_tool_tpm_normalised.csv"),merged_res,row.names = F)


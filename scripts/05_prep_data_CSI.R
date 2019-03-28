# prep data for CSI
# Philippa Borrill
# 1-6-2018


library(tidyr)
library(dplyr)
library(ggplot2)

base_dir <- "Y:\\PB_AFLF\\RefSeqv1.1_masked\\"

# load genes DE in ImpulseDE2 + gradient tool WT
WT_DE <- read.csv(file=paste0(base_dir,"control_timecourse\\impulseDE2\\merged_impuseDE2_padj0.001_gradient_tool_tpm_normalised.csv"))
head(WT_DE)
dim(WT_DE)
# now keep only DE in gradient tool and impulseDE2
WT_DE <- WT_DE[!is.na(WT_DE$pattern) & !is.na(WT_DE$padj),]
dim(WT_DE)
head(WT_DE)

### now want to add TF information
TF_v1.1 <- read.csv(file=paste0(base_dir, "control_timecourse\\gradient_tool\\control_gradient_tool_tpm_normalised\\","TFs_v1.1.csv"))
head(TF_v1.1)

## add TF information
all_DE_with_TF <- merge(WT_DE, TF_v1.1, by.x="X", by.y="gene", all.x=T)
head(all_DE_with_TF)

# how many genes are DE in one or more conditions?
nrow(all_DE_with_TF)

## how many TFs are DE 
nrow(all_DE_with_TF[!is.na(all_DE_with_TF$TF_family),])


# keep only DE in WT
DE_TFs <- all_DE_with_TF[!is.na(all_DE_with_TF$TF_family),]
dim(DE_TFs)
head(DE_TFs)

# how many NAC, MYB and WRKYs are there?
DE_TFs %>%
  filter(TF_family == "NAC"|
           TF_family=="MYB"|
           TF_family=="WRKY" ) %>%
  dim()

# how many non NAC, MYB and WRKYs are there?
DE_TFs %>%
  filter(TF_family != "NAC"&
           TF_family!="MYB"&
           TF_family!="WRKY" ) %>%
  dim()

## let's load data WT tpms and see how many TF would be included at different tpm thresholds
data_dir <- "Y:\\PB_AFLF\\RefSeqv1.1_masked\\expression_per_gene\\"

## read in data and only select genes with > 0.5 tpm ##
control.tpm.data <- read.csv(file=paste0(data_dir,"control_timecourse_tpm.tsv"), sep = "\t")
head(control.tpm.data)

# calculate max tpm for each gene
# just use FLB data
control.tpm.data <- control.tpm.data[,1:30]
head(control.tpm.data)

# average per timepoint
control.tpm.data$T3 <- (control.tpm.data[,1] + control.tpm.data[,2] + control.tpm.data[,3]) / 3
control.tpm.data$T7 <- (control.tpm.data[,4] + control.tpm.data[,5] + control.tpm.data[,6]) / 3
control.tpm.data$T10 <- (control.tpm.data[,7] + control.tpm.data[,8] + control.tpm.data[,9]) / 3
control.tpm.data$T13 <- (control.tpm.data[,10] + control.tpm.data[,11] + control.tpm.data[,12]) / 3
control.tpm.data$T15 <- (control.tpm.data[,13] + control.tpm.data[,14] + control.tpm.data[,15]) / 3
control.tpm.data$T17 <- (control.tpm.data[,16] + control.tpm.data[,17] + control.tpm.data[,18]) / 3
control.tpm.data$T19 <- (control.tpm.data[,19] + control.tpm.data[,20] + control.tpm.data[,21]) / 3
control.tpm.data$T21 <- (control.tpm.data[,22] + control.tpm.data[,23] + control.tpm.data[,24]) / 3
control.tpm.data$T23 <- (control.tpm.data[,25] + control.tpm.data[,26] + control.tpm.data[,27]) / 3
control.tpm.data$T26 <- (control.tpm.data[,28] + control.tpm.data[,29] + control.tpm.data[,30]) / 3

# just keep the average per timepoint
colnames(control.tpm.data)[31:40]
tpmData_av_con <- control.tpm.data[,31:40]
head(tpmData_av_con)

tpmData_av_con$maxtpm <- apply(tpmData_av_con[,1:10],1,max)
head(tpmData_av_con)

# now keep only tpms for DE TFs:
head(DE_TFs)
DE_TF_max_tpm <- merge(tpmData_av_con, DE_TFs, by.x=0, by.y= "Gene", all.y=T)
head(DE_TF_max_tpm)

# rename max column 
colnames(DE_TF_max_tpm)[12] <- "max"

head(DE_TF_max_tpm)
min(DE_TF_max_tpm$max)
is.na(DE_TF_max_tpm$max)

# how many TFs over tpm thresholds from 1 - 10 tpm?
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 0.5,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 1,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 2,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 3,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 4,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 5,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 6,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 7,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 8,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 9,])
nrow(DE_TF_max_tpm[DE_TF_max_tpm$max > 10,])

TFs_to_use5tpm <- DE_TF_max_tpm[DE_TF_max_tpm$max > 0.5,] # 0.5 tpm threshold
head(TFs_to_use0.5tpm)
dim(TFs_to_use0.5tpm)
colnames(TFs_to_use0.5tpm)[1] <- "gene"

## make a dataframe to run in CSI with 0.5 tpm threshold
# first for WT select the correct data
control.tpm.data <- control.tpm.data[,1:30]
head(control.tpm.data)

control.TF.0.5tpm <- control.tpm.data[ rownames(control.tpm.data) %in% TFs_to_use0.5tpm$gene,]
head(control.TF.0.5tpm)
dim(control.TF.0.5tpm)
colnames(control.TF.0.5tpm)

# re-organise
control.TF.0.5tpm <- control.TF.0.5tpm[,c(1,4,7,10,13,16,19,22,25,28,2,5,8,11,14,17,20,23,26,29,3,6,9,12,15,18,21,24,27,30)]
head(control.TF.0.5tpm)
colnames(control.TF.0.5tpm)

control.TF.0.5tpm <- rbind(c(rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),3)),
                           control.TF.0.5tpm )
colnames(control.TF.0.5tpm) <- rep(c("WT1","WT2","WT3"),each=10)
head(control.TF.0.5tpm)
dim(control.TF.0.5tpm)

write.csv(file="control_tpm_FLB_0.5tpm_input_CSI_only_WT_DE_TF.csv",control.TF.0.5tpm)


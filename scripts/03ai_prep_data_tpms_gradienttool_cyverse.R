# Script to prepare data for gradient tool on cyverse to look at genes differentially expressed across time in RNA-seq timecourse
# Philippa Borrill 08.02.2016


## read in data and only select genes with > 0.5 tpm ##

out_dir_control <- "Y:\\PB_AFLF\\RefSeqv1.1_masked\\control_timecourse\\gradient_tool\\"


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

# calculate the average max tpm
colnames(tpmData_FLB[,31:40])
tpmData_FLB$maxtpm <- apply(tpmData_FLB[,31:40],1,max)
head(tpmData_FLB)

# select only rows with a maxtpm >0.5
tpm_max_FLB <- tpmData_FLB[which(tpmData_FLB$maxtpm>0.5),]
head(tpm_max_FLB)

## remove LC genes
head(tpm_max_FLB)
tpm_max_FLB$Row.names <- rownames(tpm_max_FLB)
tpm_max_FLB <- tpm_max_FLB[!grepl("LC",tpm_max_FLB$Row.names), ]
head(tpm_max_FLB)
dim(tpm_max_FLB)

# remove tpm and maxtpm columns
tpm_max_FLB <- tpm_max_FLB[,1:30]
head(tpm_max_FLB)
dim(tpm_max_FLB)

# re-organise
tpm_FLB_ordered <- tpm_max_FLB[,c(1,4,7,10,13,16,19,22,25,28,2,5,8,11,14,17,20,23,26,29,3,6,9,12,15,18,21,24,27,30)]
head(tpm_FLB_ordered)
colnames(tpm_FLB_ordered) <- c(rep(c("3", "7", "10", "13", "15", "17", "19", "21", "23", "26"),3))
head(tpm_FLB_ordered)

write.csv(file=paste0(out_dir_control,"control_tpm_FLB_0.5tpm_input_gradient_tool.csv"),tpm_FLB_ordered)


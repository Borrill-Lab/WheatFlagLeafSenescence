# Aim is to run combine samples to gene expression level from transcript level

# Philippa Borrill
#06-03-2018

#Steps will be:

#1: Summarise counts per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

#2: Summarise tpm per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

##### #1: Summarise counts per gene ########## 

setwd("Y:\\PB_AFLF\\RefSeqv1.1_masked\\")

#source("https://bioconductor.org/biocLite.R")
#biocLite("tximportData")
#install.packages("readr")
library(tximportData)
library(readr)

# read in list of studies to use
study_list <- c("control_timecourse")
study_list

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.table("Y:\\PB_AFLF\\RefSeqv1.1\\transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)

for (i in study_list) {
  study <- i
setwd(paste0("Y:\\PB_AFLF\\RefSeqv1.1_masked\\",study))

#study <- "control_timecourse"

# make vector pointing to the kallisto results files   ########
samples <- read.table(paste0(study,"_samples.txt"), header=F)
samples
colnames(samples) <- c("sample")

files <- file.path("Y:\\PB_AFLF\\RefSeqv1.1_masked", study, "kallisto_results",samples$sample, "abundance.tsv", fsep ="\\")
files
names(files) <- paste0(samples$sample)
head(files)
all(file.exists(files))

library(tximport)
# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
names(txi)

# move into directory where I will save this analysis
setwd("Y:\\PB_AFLF\\RefSeqv1.1_masked\\expression_per_gene")

# to see counts summarised per gene
#head(txi$counts)
#colnames(txi$counts)

# save counts summarised per gene
write.table(txi$counts, file=paste0(study,"_count.tsv"),sep = "\t")

# to see tpm summarised per gene
#head(txi$abundance)
#colnames(txi$abundance)

# save tpm summarised per gene
write.table(txi$abundance, file=paste0(study,"_tpm.tsv"),sep = "\t")

# see lengths summarised per gene
head(txi$length)

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file=paste0(study,"_gene_lengths.csv"))

}

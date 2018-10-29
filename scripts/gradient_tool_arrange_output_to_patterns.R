# Aim is to look at gradient tool output and group genes into clusters which have the same expression pattern
# Philippa Borrill 
# 29.10.2018

## first for normalised

setwd("Y:\\PB_AFLF\\RefSeqv1.1_masked\\control_timecourse\\gradient_tool\\control_gradient_tool_tpm_normalised\\")
out_dir <- "Y:\\PB_AFLF\\RefSeqv1.1_masked\\control_timecourse\\gradient_tool\\control_gradient_tool_tpm_normalised\\GO_enrichment_zscore2_clusters\\"

data <- read.csv(file="out.csv") # this is the file produced by Cyverse
head(data)
dim(data)

# what I want to do is re-arrange the data so I have 1 row per gene, each DAA (X) as a column, with the tscore (i.e. zscore) in the box

library("tidyr")

data_sel <- data[,c(1,2,7)]
head(data_sel)

data_spread <- spread(data_sel, X, tscore)
head(data_spread)

data_spread[data_spread$item =="TraesCS1D02G327600",]

# now convert tscore to downreg (-1), upreg (+1) or flat (0)
up_down_data <- data_spread
head(up_down_data)
row.names(up_down_data) <- up_down_data$item
up_down_data <- up_down_data[,-1]
head(up_down_data)

# replace anything between -2 and 2 with 0 # 2 is the zscore threshold I will use
up_down_data[up_down_data <= 2 & up_down_data >= -2] <- 0
head(up_down_data)
tail(up_down_data)
up_down_data[up_down_data < -2 ] <- -1
head(up_down_data)
tail(up_down_data)
up_down_data[up_down_data > 2 ] <- 1
head(up_down_data)
tail(up_down_data)

tail(up_down_data,100)
up_down_data[up_down_data$`10` > 0,]

# now sort according to each column in turn
sorted_up_down_data <- up_down_data[
  order(up_down_data[,1],
        up_down_data[,2],
        up_down_data[,3],
        up_down_data[,4],
        up_down_data[,5],
        up_down_data[,6],
        up_down_data[,7],
        up_down_data[,8],
        up_down_data[,9],
        up_down_data[,10]),]

head(sorted_up_down_data)
tail(sorted_up_down_data)
dim(sorted_up_down_data)

# use unite to make a new column of the "pattern"
sorted_up_down_data <- (unite(sorted_up_down_data,pattern, remove=F, `3`:`26`, sep="_"))
head(sorted_up_down_data)
tail(sorted_up_down_data)

dim(unique(sorted_up_down_data[,1:10]))
length(unique(sorted_up_down_data$pattern))
# there are 213 different patterns
write.csv(sorted_up_down_data,file="sorted_up_down_data_HC_only.csv")

# how many genes are in each pattern (not much point in keeping patterns with very few genes - should group them together- but how to decide that?)

library("dplyr")

number_per_pattern <- sorted_up_down_data %>%
  group_by(pattern) %>%
  summarize(n= length(`3`))

write.csv(file=paste0(out_dir,"number_per_pattern.csv"), number_per_pattern)

head(number_per_pattern)
nrow(number_per_pattern) # number of patterns
nrow(as.data.frame(number_per_pattern)[number_per_pattern$n > 10,]) # number of patterns with >10 genes
nrow(as.data.frame(number_per_pattern)[number_per_pattern$n > 20,]) # number of patterns with >10 genes

library(lumi)
library(edgeR)
#rawdata <- read.csv("miRNA_raw_counts_HEK_Pellet.csv", header= TRUE, row.names=1)
#rawdata<- read.csv("miRNA_raw_counts_Correct_Name_ReOrdered.csv", header= TRUE, row.names= 1)
rawdata<- read.delim("miRNA_raw_counts_CORRECT_NAMES.txt", header= TRUE, row.names= 1)
rawdata<- rawdata[c(43:54)]
rawdata<-rawdata[,order(names(rawdata))]
#Below is only PC3 data
# Filtering and Normalization
# compute the effective library size by using TMM normalization
keep <- rowSums(rawdata>10) >=2
y <- rawdata[keep,]
write.table(y, "PC3exoFiltered.txt", sep="\t")

##To Filter and subset data 
# Data input
# load the edgeR library
library(edgeR)
library(lumi)
# read in the raw data
rawdata <- read.delim("miRNA_raw_counts_CORRECT_NAMES.txt")
# create a DGEList object. In this case, raw data is in columns 2 to 55, and unique identifiers are in 1 (miRNA name)
y <- DGEList(counts=Cav1PN4[,1:6], genes=rownames(rawdata))
# Filtering and Normalization
# compute the effective library size by using TMM normalization
keep <- rowSums(cpm(y)>10) >=2
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)

#make QC plots
d <-cpm(y, normalized.lib.sizes=TRUE)
d<-t(d)
dist<-dist(d)
hc<-hclust(dist)
plot(hc)
# Export the graph to a PDF
pdf(file="dendrogramPC3Cav1Sub.pdf", paper="a4r")
plot(hc)
dev.off()

# An MDS plot shows the biological coefficient of variation between the samples. The two 
# dimensions are the biggest and second biggest sources of variation within the data. Again, this plot 
# seems to suggest that the biggest differentiators of these samples are not related to the properties of 
# the genes tested.
d <-t(d)
plotMDS(d, col=c(rep(1,15), rep(2, 15), rep(1, 12), rep(2, 12))) 
legend("topright", legend = c("Exosome", "Pellet"), col = 1:2, pch = 15) 
plotMDS(d, col=c(rep(1,30), rep(2,24))) 
legend("topright", legend = c("HEK", "PC3"), col = 1:2, pch = 15) 
# Export the graph to a PDF
pdf(file="MDSCavin1withSubPC3exo.pdf", paper="a4r")
plotMDS(d, col=c(rep(1,15), rep(2, 15), rep(1, 12), rep(2, 12))) 
legend("topright", legend = c("Exosome", "Pellet"), col = 1:2, pch = 15) 

plotMDS(d, col=c(rep(1,30), rep(2,24))) 
legend("topright", legend = c("HEK", "PC3"), col = 1:2, pch = 15)
dev.off()


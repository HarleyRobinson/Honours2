#Here, I'll look at the top expressed genes that are 
#also identified through prediction. 
data=read.csv("PC3PelletC1DownSigAKA363.csv", header=TRUE)
keep<- data[data$log2FoldChange>=0, ]
keep<- keep[order(keep$log2FoldChange, decreasing= FALSE),]
write.csv(keep, "OrderedmiR363proteins.csv")

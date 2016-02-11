library(DESeq2)
#Always filter and normalise the data first (FilteringandSubsetting.r)
rawdata<- read.csv(file="miRNA_analysis_kerry.csv", header= TRUE, row.names=1)
Cav2PN4= rawdata[c(13:15, 1:3)]
Cav1PN4= rawdata[c(1:3, 5, 6, 7)]
CaveolinPN4= rawdata[c(1:3, 4:6)]
Cav3PN4= rawdata[c(7:9, 1:3)]

#Analysis for Cavin
exp_designCav = data.frame(row.names = colnames(Cav1PN4),
                            condition =c("Test", "Test", "Test", "Cont", "Cont", "Cont"),
                            libType = c("paired-end", "paired-end", "paired-end", "paired-end", "paired-end", "paired-end"))
condition= relevel(exp_designCav$condition, "Test")
levels(exp_designCav)
head(exp_designCav)
ddsCav=DESeqDataSetFromMatrix(countData=Cav1PN4, colData=exp_designCav, design=~condition)
ddsCav=DESeq(ddsCav)
resCav= results(ddsCav, alpha=0.05)
write.csv(resCav, "HEKEXOCaveVCav3.csv")

rlogoddsdds=rlogTransformation(ddsCav1)
png("plotPCA_Test.png")
plotPCA(rlogoddsdds)
text(rownames(rlogoddsdds), pos= 3)
dev.off()

#graph below gives best mirna given high BHpval annnnd pvalue. 
png("plotMA_CaveVcav3HEK_minpval.png")
DESeq2::plotMA(resCav, ylim=c(-3, 3))
low<- resCav[order(resCav$pvalue), ]
topGene<- rownames(resCav)[which(low$pvalue <=0.05)]
with(resCav[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()
exo<- read.csv("EXOcavin1VGFP.csv", header=TRUE, row.names= 1)
cell <- read.csv("CELLcavin1VGFP.csv", header= TRUE, row.names= 1)
exoMirs<- rownames(exo)
cell<- cell[exoMirs, ]
FCFC<- as.data.frame(cell$log2FoldChange-exo$log2FoldChange)
rownames(FCFC)<- exoMirs
FCFC2<- as.data.frame(exo$log2FoldChange-cell$log2FoldChange)
rownames(FCFC2)<- exoMirs
write.csv(FCFC, "FCFCminusAll.csv")
FCFC<- na.omit(FCFC)
fc<- FCFC$`cell$log2FoldChange - exo$log2FoldChange`
range(fc)
breaks= seq(-1.2, 1.2, by=0.2)
fc.cut=cut(fc, breaks, right=FALSE)
fc.freq= table(fc.cut)
downers<- exo$log2FoldChange<0
MirSeqs<- read.csv("MirSeqNameVsSeqSep.csv", header=FALSE)
seq<- 0
d= NULL
i=1
for (i in(1:nrow(exoMirs))) {
  matched= grep(exoMirs[i, 1], MirSeqs$V1)
  print(exoMirs[i, 1])
  print(MirSeqs[matched, 2])
  #d= rbind(d, data.frame(exoMirs[i], MirSeqs[matched, 2]))
  i=i+1
  }

write.csv(d, "NamesOfUsedMirs.csv")
names<- read.csv("NamesOfUsedMirs.csv", header= FALSE, row.names=1)
read<- read.csv("FCFCminus2.csv", header= TRUE)
exoMirs<- read.csv("MirsUpLookingFornames.csv", header= FALSE)

i=1
for (i in range(1:100)) {
  mir<- sample(mir)
  i=i+1
}
print(mir)
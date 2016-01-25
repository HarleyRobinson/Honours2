library(topGO)
library(annotate)
library(org.Hs.eg.db)
predicted<- read.delim("PCcav2DownGenes.txt", header= FALSE)
Data<- read.csv("deg_cavin1_gfp.csv", header= TRUE)
PredictedGenes<- predicted$V2
DataGenes<- Data$symbol
Overlap<- unique(PredictedGenes[DataGenes%in% PredictedGenes])
CodingGenes<- Data[Data$biotype=="protein_coding",]
SigCG<- CodingGenes[CodingGenes$padj<=0.1, ]
UsefulGenes<- SigCG[SigCG$symbol%in%Overlap, ]
write.csv(UsefulGenes, "PC3PelletC2DownSig.csv")
Filtered<- UsefulGenes[UsefulGenes$log2FoldChange>=0, ]
names<- get(rownames(predictions), data= 'org.Hs.eg')
GOdata= new('topGOdata', ontology= 'BP', allGenes= Filtered$symbol, annot= annFUN.GO2genes, GO2genes = allGoGenes, geneSel= selection, nodeSize=10)

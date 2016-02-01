library(topGO)
library(annotate)
library(org.Hs.eg.db)

predicted<- read.delim("SharedmiR20b5p.txt", header= FALSE)
Data<- read.csv("deg_cavin1_gfp.csv", header= TRUE)
PredictedGenes<- predicted$V2
DataGenes<- Data$symbol
Overlap<- unique(PredictedGenes[PredictedGenes%in%DataGenes])
CodingGenes<- Data[Data$biotype=="protein_coding",]
SigCG<- CodingGenes[CodingGenes$padj<=0.1, ]
rownames(SigCG)<- make.names(SigCG[, 2], unique=TRUE)
UsefulGenes<- Overlap[Overlap%in%SigCG$symbol]
rownames(Data)<- make.names(Data[, 2], unique=TRUE)
Genes<- Data[as.character(UsefulGenes), ]
write.csv(Genes, "OrderedmiR20proteins.csv")
Filtered<- UsefulGenes[UsefulGenes$log2FoldChange>=0, ]
#names<- get(rownames(predictions), data= 'org.Hs.eg')
#GOdata= new('topGOdata', ontology= 'BP', allGenes= Filtered$symbol, annot= annFUN.GO2genes, GO2genes = allGoGenes, geneSel= selection, nodeSize=10)

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") 
results <- getBM(attributes = c("go_id"), filters = "refseq_mrna", values = c("9606.ENSP00000404029"), mart = mart) 
nipps<- 0
for (name in plist) {
  name= as.character(name)
  nipps[name] <- getBM(attributes = c("go_id"), filter="uniprot_swissprot", values = c(name), mart = mart)
}
name<- "Q9NVA2"
results <- getBM(attributes = c("go_id"), filter="uniprot_swissprot", values = c("Q53H82"), mart = mart) 
GORNAbind<- 0
i= 0
for (list in nipps) {
  GORNAbind[i]<- "GO:0003723"%in%list
  i=i+1
}
#GO:0003723 for RNA bind
#GO:0045121 for membrane raft
#GO:0016020 for membrane 
ProtBind<-which(GORNAbind ==1, arr.ind=TRUE)
Candidate<- plist[ProtBind]

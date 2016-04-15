FusBind<- read.csv("FUSinter2.csv", header=FALSE)
DRMpresent<- read.csv("RaftSigParsed.csv", header=TRUE) 
DRMdata<- read.csv("lipidraftSig.csv", header=TRUE)
FinalGraph= cbind(DRMdata, DRMpresent)
write.csv(FinalGraph, "LRsigGeneName.csv")
keep<- intersect(DRMpresent$Gene.Name, FusBind$V8)
keep2<- intersect(DRMpresent$Gene.Name, FusBind$V9)
rownames(FinalGraph)= FinalGraph$Gene.Name
FusBindDRMFound<- FinalGraph[keep2,]
write.csv(FusBindDRMFound, "FusBindingDRMfoundproteins.csv")
Kbind<- read.csv("HNRNPKinter.csv", header=FALSE)
Keep<- intersect(Kbind$V9, DRMpresent$Gene.Name)
KbindDRMfound<- FinalGraph[Keep, ]
write.csv(KbindDRMfound, "KbindingDRMfoundProteins.csv")

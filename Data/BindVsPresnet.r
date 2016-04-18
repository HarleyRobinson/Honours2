FusBind<- read.csv("FUSinter2.csv", header=FALSE)
DRMpresent<- read.csv("RaftinsigParsed.csv", header=FALSE) 
DRMdata<- read.csv("DRMinsig.csv", header=TRUE)
FinalGraph= cbind(DRMdata, DRMpresent)
write.csv(FinalGraph, "LRinsigGeneName.csv")
keep<- intersect(DRMpresent$V1, FusBind$V8)
keep2<- intersect(DRMpresent$V1, FusBind$V9)
rownames(FinalGraph)= FinalGraph$V1
FusBindDRMFound<- FinalGraph[keep2,]
write.csv(FusBindDRMFound, "FusBindingDRMfoundproteinsINSIG.csv")
Kbind<- read.csv("HNRNPKinter.csv", header=FALSE)
Keep<- intersect(Kbind$V9, DRMpresent$V1)
KbindDRMfound<- FinalGraph[Keep, ]
write.csv(KbindDRMfound, "KbindingDRMfoundProteinsINSIG.csv")


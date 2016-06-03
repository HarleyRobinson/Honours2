library(topGO)
#using files generated from ComparingPredicted...
data=read.csv("PartnersofInterestingMirs.csv", header=FALSE)
#Sorting for only genes that are expected to be changed in expression
#For Upreg miR, protein should be negative and vice versa. 
keep<- data[1:819, ]
#write.csv(keep, "PC3PelletmiR20bfiltered.csv")
#next part needs annotate. 

all<- factor(as.integer(data[, 1]%in%keep[, 1]))
names(all)<- data[, 1]
selection= function(x) TRUE
GOdata<- new("topGOdata", ontology="BP", allGenes= all, geneSel=function(p) p==1, description= "miR targets", annot=annFUN.org, mapping="org.Hs.eg.db", ID= "symbol")
results.ks= runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks
allRes= GenTable(GOdata, KS= results.ks, orderBy= 'KS', topNodes=25)
allRes[,c('GO.ID', 'Term', 'KS')]
write.csv(allRes[,c('GO.ID', 'Term', 'KS')], "InterestingMirsenrichment.csv")

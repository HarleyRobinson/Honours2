

pred6b= getPredictedTargets(mir6, sources= c("miranda", "targetscan"), species ='hsa', method = 'max', min_src=2)
pred5b= getPredictedTargets(mir5, sources= c("miranda", "targetscan"), species ='hsa', method = 'max', min_src=2)
pred4b= getPredictedTargets(mir4, sources= c("miranda", "targetscan"), species ='hsa', method = 'max', min_src=2)
pred3b= getPredictedTargets(mir3, sources= c("miranda", "targetscan"), species ='hsa', method = 'max', min_src=2)
pred20bb= getPredictedTargets(mir2, sources= c("miranda", "targetscan"), species ='hsa', method = 'max', min_src=2)
predictionsb = getPredictedTargets(mir, sources= c("miranda", "targetscan"), species ='hsa', method = 'max', min_src=2)
totalb<- rbind(predictionsb, pred20bb, pred3b, pred4b, pred5b, pred6b)

rankedGenes<- totalb[, 'rank_product']
selection= function(x) TRUE
allGoGenes= annFUN.org(whichOnto= 'BP', feasibleGenes= NULL, mapping= 'org.Hs.eg.db', ID= 'entrez')
GOdata= new('topGOdata', ontology= 'BP', allGenes= rankedGenes, annot= annFUN.GO2genes, GO2genes = allGoGenes, geneSel= selection, nodeSize=10)

results.ks= runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks
allRes= GenTable(GOdata, KS= results.ks, orderBy= 'KS', topNodes=20)
allRes[,c('GO.ID', 'Term', 'KS')]
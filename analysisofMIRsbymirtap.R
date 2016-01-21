library(miRNAtap)
library(org.Hs.eg.db)
library(annotate)
library(topGO)
mir= 'miR-363'
predictions = getPredictedTargets(mir, species ='hsa', method = 'max', min_src=2)
rownames(predictions)<- getSYMBOL(rownames(predictions), data= 'org.Hs.eg')
predictions<- cbind(predictions, c(rep(mir, nrow(predictions))))


mir2= 'miR-20b'
pred20b= getPredictedTargets(mir2, species ='hsa', method = 'max', min_src=2)
rownames(pred20b)<- getSYMBOL(rownames(pred20b), data= 'org.Hs.eg')
pred20b<- cbind(pred20b, c(rep(mir2, nrow(pred20b))))


mir3= "miR-106a"
pred3= getPredictedTargets(mir3, species ='hsa', method = 'max', min_src=2)
rownames(pred3)<- getSYMBOL(rownames(pred3), data= 'org.Hs.eg')
pred3<- cbind(pred3, c(rep(mir3, nrow(pred3))))

mir4= 'miR-17'
pred4= getPredictedTargets(mir4, species ='hsa', method = 'max', min_src=2)
rownames(pred4)<- getSYMBOL(rownames(pred4), data= 'org.Hs.eg')
pred4<- cbind(pred4, c(rep(mir4, nrow(pred4))))

mir5= 'miR-20a'
pred5= getPredictedTargets(mir5, species ='hsa', method = 'max', min_src=2)
rownames(pred5)<- getSYMBOL(rownames(pred5), data= 'org.Hs.eg')
pred5<- cbind(pred5, c(rep(mir5, nrow(pred5))))

mir6= 'miR-125b'
pred6= getPredictedTargets(mir6, species ='hsa', method = 'max', min_src=2)
rownames(pred6)<- getSYMBOL(rownames(pred6), data= 'org.Hs.eg')
pred6<- cbind(pred6, c(rep(mir6, nrow(pred6))))

#mir7= 'miR-203a'
#pred7= getPredictedTargets(mir7, sources= c("miranda", "targetscan"), species ='hsa', method = 'max', min_src=2)
#rownames(pred7)<- getSYMBOL(rownames(pred7), data= 'org.Hs.eg')
#pred7<- cbind(pred7, c(rep(mir7, nrow(pred7))))


total<- rbind(predictions, pred20b, pred3, pred4, pred5, pred6)
write.csv(total, 'mir-363predictions.csv', row.names=TRUE) 

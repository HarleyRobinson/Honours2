rd<- read.csv("miRNA_analysis_kerry.csv", header=TRUE, row.names= 1)
rawdata<- rd[c(1:3, 5, 6, 8)]
rdPellet<- read.delim("miRNA_raw_counts_CORRECT_NAMES.txt", header=TRUE, row.names=1)
Pellet<- rdPellet[c(31:42)]
Pellet<- Pellet[c(2, 6, 10, 1, 5, 9)]
name<- intersect(rownames(Pellet), rownames(rawdata))
data<- Pellet[cbind(name), ]
dataPellet<- data[order(rownames(data)), ]
rawdata<- rawdata[cbind(name), ]
rawdata<- rawdata[order(rownames(data)), ]
fullset<- cbind(rawdata, dataPellet)
rf<- as.data.frame(rowMeans(rawdata[4:6]))
rg<- as.data.frame(rowMeans(dataPellet[4:6]))
Fullmeans<- as.data.frame(c(rf, rg))
rownames(Fullmeans)<- rownames(rf)
plot(rf$`rowMeans(rawdata[4:6])`, rg$`rowMeans(dataPellet[4:6])`, main="GFP Correlation", xlab="Exo", ylab="Cell", pch=19)
text(rf$`rowMeans(rawdata[4:6])`, rg$`rowMeans(dataPellet[4:6])`, label=rownames(rf), cex= 0.7)
#making normalised data
subsetCell<- as.data.frame(rg[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"),])
rownames(subsetCell)<- c("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p")
subsetExo<- as.data.frame(rf[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"),])
rownames(subsetExo)<- c("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p")
p<- plot(subsetExo$`rf[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, subsetCell$`rg[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, main="ex", xlab="Exo", ylab="Cell", pch=19)
text(subsetExo$`rf[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, subsetCell$`rg[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, labels=rownames(subsetExo), cex=0.7, pos=3)
p + geom_label_repel()
reg1<- lm(rg$`rowMeans(dataPellet[4:6])`~rf$`rowMeans(rawdata[4:6])`-1)
abline(reg1)
textplot(rf[imp, ], rg[imp, ], imp)
library(ggrepel)
imp<- c("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p")
text(rf[imp,], rg[imp,], labels=imp)
#####Attempting to do FC correlation plots
Cell<-read.csv("CELLcavin1VGFP.csv", header= TRUE, row.names= 1)
Exo<- read.csv("EXOcavin1VGFP.csv", header= TRUE, row.names= 1)
name<- intersect(rownames(Cell), rownames(Exo))
data<- Cell[cbind(name), ]
CellD<- data[order(rownames(data)), ]
ExoD<- Exo[name,]
ExoD<- ExoD[order(rownames(ExoD)), ]
plot(CellD$log2FoldChange, ExoD$log2FoldChange, main= "MicroRNA expression FC between GFP and Cavin-1", xlab="CELL FC", ylab= "EXO FC", pch=19, xlim=c(-3.1, 1.5), ylim=c(-3.1, 1.5))
reg1<- lm(CellD$log2FoldChange~ExoD$log2FoldChange)
#summary(lm(CellD$log2FoldChange~ExoD$log2FoldChange))
#abline(0, 0.9)
abline(0.4, 0.9)
abline(-0.4,0.9)
text(CellD[imp, "log2FoldChange"], ExoD[imp, "log2FoldChange"], imp, ylim=c(-3.5, 1.5), xlim= c(-3.5, 1.5), cex = 1.2)

#tring to get the labels to repel
install.packages(c("wordcloud","tm"),repos="http://cran.r-project.org")
library(wordcloud)
library(tm)

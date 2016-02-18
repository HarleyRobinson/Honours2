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
rf<- as.data.frame(rowMeans(rawdata[1:3]))
rg<- as.data.frame(rowMeans(dataPellet[1:3]))
plot(rf$`rowMeans(rawdata[1:3])`, rg$`rowMeans(dataPellet[1:3])`, main="eg", xlab="Exo", ylab="Cell", pch=19)
text(rf$`rowMeans(rawdata[1:3])`, rg$`rowMeans(dataPellet[1:3])`, label=rownames(rf), cex= 0.7)
#making normalised data
subsetCell<- as.data.frame(rg[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"),])
rownames(subsetCell)<- c("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p")
subsetExo<- as.data.frame(rf[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"),])
rownames(subsetExo)<- c("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p")
plot(subsetExo$`rf[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, subsetCell$`rg[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, main="ex", xlab="Exo", ylab="Cell", pch=19)
text(subsetExo$`rf[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, subsetCell$`rg[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`, labels=rownames(subsetExo), cex=0.7, pos=3)
reg1<- lm(subsetCell$`rg[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`~subsetExo$`rf[cbind("hsa-miR-363-3p", "hsa-miR-20b-5p", "hsa-miR-19a-3p", "hsa-miR-148a-3p", "hsa-miR-574-3p", "hsa-miR-146a-5p", "hsa-miR-30a-5p", "hsa-miR-10b-5p", "hsa-miR-200a-3p"), ]`)
abline(reg1)

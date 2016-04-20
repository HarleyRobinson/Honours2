
pwm<- read.csv("IdealMotifAlldownpesudo.csv", header= TRUE, row.names=1)
seq<- read.csv("seqs.csv", header=FALSE, row.names=1)
n=1
windowscore<-0
for (n in 1:(length(seq)-4)) {
  window<- seq[25, n:(n+3)]
  print(rownames(seq[25,]))
  i=1
  d=0
  for (column in window) {
    d[i]<-pwm[as.character(column), i]
    i=i+1
  }
  windowscore[n]<- d[1]*d[2]*d[3]*d[4]
  n=n+1
}
windowscore
max(windowscore)
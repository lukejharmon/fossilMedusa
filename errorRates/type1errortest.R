source("errorRates/newbirthdeathtree.R")
library(ape)
library(fossilMEDUSA)

numberOfBreaks<-numeric(1000)
for(i in 1:1000) {
  
tt<-birthdeathTree(b=1, d=0.5, taxaStop=100, return.all.extinct=F)

cc<-collapseTreeAndFossils(tt, endClades=10, timeFraction=0)
cc

#layout(matrix(1:3, ncol=1))
#plot(tt$tree, show.tip.label=F)
#plot(cc$meduTree)
#plotDtt(cc$fullDtt, fullTime=tt$t)

rr<-makeMedusaRichness(cc)

res<-runFossilMEDUSA(cc$meduTree, rr)
ss<-summarize.MEDUSA(res)
numberOfBreaks[i]<-ss
}

table(numberOfBreaks)

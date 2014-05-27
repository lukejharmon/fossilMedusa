source("errorRates/newbirthdeathtree.R")
library(ape)

tt<-birthdeathTree(b=1, d=0.5, taxaStop=100)

cc<-collapseTreeAndFossils(tt, endClades=5, timeFraction=0)
cc

layout(matrix(1:3, ncol=1))
plot(tt$tree, show.tip.label=F)
plot(cc$meduTree)
plotDtt(cc$fossilDtt)

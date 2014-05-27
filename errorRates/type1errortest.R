source("errorRates/newbirthdeathtree.R")

tt<-birthdeathTree(b=1, d=0, taxaStop=100)
plot(tt$tree)

cc<-collapseTreeAndFossils(tt, endClades=10)
cc

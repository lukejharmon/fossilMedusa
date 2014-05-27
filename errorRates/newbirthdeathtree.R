`birthdeathTree` <-
function (b, d, timeStop=0, taxaStop=0, seed=0, print.seed=FALSE, return.all.extinct=TRUE, keepAppearances=TRUE){

	if(seed==0) seed=set.seed.clock(print=print.seed);

	if(timeStop==0 & taxaStop==0)
	stop("Must have stopping criterion\n");
	

	
	while(1) {

		edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
		edge.length <- rep(NA, 2)
		stem.depth <- numeric(2)
		alive<-rep(TRUE, 2) # marker for live lineages
		t <- 0 #time at any point in the tree
		next.node<-4
		if(keepAppearances) {
			app<-cbind(c(2, 3), c(0, 0), c(NA, NA))
		}


		############
		repeat{
			if(taxaStop) if(sum(alive)>=taxaStop) break;
			if(sum(alive)==0) break;
			dt<-rexp(1, sum(alive)*(b+d));
	 	 	t<-t+dt;
	  		if(timeStop) if(t>=timeStop) {
				t<-timeStop;
				break;
	  		}
	  		r<-runif(1)
         	if(r<=b/(b+d)) {#this creates a bifucation in the tree
	        	random_lineage <- round(runif(1, min=1, max=sum(alive)))
				e<-matrix(edge[alive,], ncol=2)
				parent<-e[random_lineage,2]
				alive[alive][random_lineage]<-FALSE
				edge<-rbind(edge, c(parent, next.node), c(parent, next.node+1))
				if(keepAppearances) {
					app<-rbind(app, c(next.node+1, t, NA))
					app[app[,1]==parent,1]<-next.node
				}
				next.node<-next.node+2
				alive<-c(alive, TRUE, TRUE)
				stem.depth<-c(stem.depth, t, t)
				x<-which(edge[,2]==parent)
				edge.length[x]<-t-stem.depth[x]
				edge.length<-c(edge.length, NA, NA)
			
            }

         	else {#This terminates one of the current lineages on the tree
            	random_lineage <- round(runif(1, min=1, max=sum(alive)))
		    	edge.length[alive][random_lineage]<-t-stem.depth[alive][random_lineage];
          	    if(keepAppearances) {
          	    	killed<-edge[alive,2][random_lineage]
					app[app[,1]==killed,3]<-t
				}
				alive[alive][random_lineage]<-FALSE

            }
      	}
      
		if(return.all.extinct==T | sum(alive)>1) break;
	}
	
	edge.length[alive]<-t-stem.depth[alive]
	

	tips<-!(edge[,2] %in% edge[,1])
	tt<-edge[tips,2]
	
	it<-as.numeric(levels(as.factor(edge[,1])))
	
	tCount<-1
	iCount<-1
	
	newEdge<-matrix(nrow=dim(edge)[1], ncol=dim(edge)[2])
	newApp<-matrix(nrow=dim(app)[1], ncol=dim(app)[2])

	for(i in 1:max(edge)) {
		if(i %in% tt) {
			rr<-tCount
			tCount<-tCount+1
		} else {
			rr<-iCount+sum(tips)
			iCount<-iCount+1		
		}
		ok<-edge==i
		newEdge[ok]<-rr
		newApp[which(app[,1]==i), 1]<-rr


	}
			
	edge<-newEdge
	app[,1]<-newApp[,1]



	tip.label<-edge[tips,2]
		
	
    	mode(tip.label) <- "character"
    	obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label, Nnode=length(tip.label)-1)
    	class(obj) <- "phylo"


    
    if(keepAppearances) {
    	return(list(tree=obj, app=app, t=t))	
    } else {
    	return(obj)
    }

}	
	
`set.seed.clock` <-
function(print=F){
	date = date()
 	seed1 = as.numeric(strsplit(substring(date,12,19),":")[[1]])%*%c(1,100,10000)
 	seed <- runif(1, min=0, max=50) * seed1
 	set.seed(seed)
 	if(print) cat("Seed = ", seed, "\n");
 	seed[1,1]
}

plotIntervals<-function(bb) {
	pres<-bb$t

	plot(bb$app[,2], 1:dim(bb$app)[1], xlim=c(0, pres), pch=19, cex=0.3)


	for(i in 1:dim(bb$app)[1]) {
		lines(bb$app[i,2:3], c(i,i))
		if(is.na(bb$app[i,3]))
			lines(c(bb$app[i,2],pres), c(i,i), col="red")
	}
}

`bdVarTree` <-
function (b, d, change, timeStop=0, taxaStop=0, seed=0, print.seed=FALSE, return.all.extinct=TRUE, keepAppearances=FALSE, rateCats=c(0.5, 1, 2)){

	if(seed==0) seed=set.seed.clock(print=print.seed);

	if(timeStop==0 & taxaStop==0)
		stop("Must have stopping criterion\n");
	

	
	while(1) {

		edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
		edge.length <- rep(NA, 2)
		stem.depth <- numeric(2)
		sRates<-rep(1, 2)
		eRates<-rep(1, 2)
		alive<-rep(TRUE, 2)
		t <- 0 #time at any point in the tree
		next.node<-4
		if(keepAppearances) {
			app<-cbind(c(2, 3), c(0, 0), c(NA, NA))
		}
		branchChanges<-rbind(c(1, 0, 1, 1, 1), c(1, 0, 2, 1, 1))
		colnames(branchChanges)<-c("Lineage", "Time", "Type", "Beginning", "End")

		############
		repeat{
			if(taxaStop) if(sum(alive)>=taxaStop) break;
			if(sum(alive)==0) break;
			rateVector<-c(sRates[alive]*b, eRates[alive]*d, rep(change, sum(alive)))
			dt<-rexp(1, sum(rateVector));
	 	 	t<-t+dt;
	  		if(timeStop) if(t>=timeStop) {
				t<-timeStop;
				break;
	  		}
	  		probVector<-cumsum(rateVector/sum(rateVector))

	  		r<-runif(1)
	  		event<-min(which(probVector>r))
	  	
	  		
         	if(event<=sum(alive)) {#this creates a bifucation in the tree
	        	random_lineage <- event
				e<-matrix(edge[alive,], ncol=2)
				parent<-e[random_lineage,2]
				
				ps<-sRates[alive][random_lineage]
				pe<-sRates[alive][random_lineage]
				sRates<-c(sRates, ps, ps)
				eRates<-c(eRates, pe, pe)
				
				alive[alive][random_lineage]<-FALSE
				edge<-rbind(edge, c(parent, next.node), c(parent, next.node+1))
				if(keepAppearances) {
					app<-rbind(app, c(next.node+1, t, NA))
					app[app[,1]==parent,1]<-next.node
				}
				next.node<-next.node+2

				
				alive<-c(alive, TRUE, TRUE)
				stem.depth<-c(stem.depth, t, t)
				x<-which(edge[,2]==parent)
				edge.length[x]<-t-stem.depth[x]
				edge.length<-c(edge.length, NA, NA)

			
            } else if(event<=sum(alive)*2) {#This terminates one of the current lineages on the tree
            	random_lineage <- event-sum(alive)
		    	edge.length[alive][random_lineage]<-t-stem.depth[alive][random_lineage];
          	    if(keepAppearances) {
          	    	killed<-edge[alive,2][random_lineage]
					app[app[,1]==killed,3]<-t
				}
				alive[alive][random_lineage]<-FALSE

            } else {
            	random_lineage <- event-2*sum(alive)
            	rr<-runif(1)
            	if(rr<0.5) {
            		et<-1
            		or<-sRates[alive][random_lineage]
            		while(1) {
            			nr<-sample(rateCats)[1]
            			if(nr!=or) break
            			}
            		sRates[alive][random_lineage]<-nr
				} else {
					et<-2
					or<-eRates[alive][random_lineage]
					while(1) {
						
						}
					eRates[alive][random_lineage]<-sample(rateCats)[1]
					nr<-eRates[alive][random_lineage]
				}
				if(or!=nr)
					branchChanges<-rbind(branchChanges, c(edge[alive,2][random_lineage], t, et, or, nr))
            }
      	}
      
		if(return.all.extinct==T | sum(alive)>1) break;
	}
	edge.length[alive]<-t-stem.depth[alive]
	
	
	n<--1;
	for(i in 1:max(edge)) {
		if(any(edge[,1]==i)) {
			edge[which(edge[,1]==i), 1]<-n
			edge[which(edge[,2]==i), 2]<-n
			app[which(app[,1]==i), 1]<-n
			n<-n-1
		}
	}

	mm<-match(app[,1], edge[,2])
	edge[edge[,2]>0,2]<-1:sum(edge>0)


	tips<-!(edge[,2] %in% edge[,1])
	tip.label<-edge[tips,2]
	
	app[,1]<-edge[mm,2]

	#tip.label<-1:sum(edge>0)
	
	
    	mode(tip.label) <- "character"
    	obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label)
    	class(obj) <- "phylo"
    	obj<-old2new.phylo(obj)
    	obj<-read.tree(text=write.tree(obj))
    
    if(keepAppearances) {
    	return(list(tree=obj, app=app, t=t, branchChanges=branchChanges))	
    } else {
    	return(obj)
    }
}

`bdVarTreeSetTimes` <-
function (b=0.1, d=0.05, timeStop=15, taxaStop=0, seed=0, print.seed=FALSE, return.all.extinct=TRUE, keepAppearances=FALSE, timesToChange=c(5, 6), newS=c(0.01, 10), newE=c(0.01, 0.01)) {

	if(seed==0) seed=set.seed.clock(print=print.seed);

	if(timeStop==0 & taxaStop==0)
		stop("Must have stopping criterion\n");
	

	
	while(1) {

		edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
		edge.length <- rep(NA, 2)
		stem.depth <- numeric(2)
		sRates<-rep(1, 2)
		eRates<-rep(1, 2)
		alive<-rep(TRUE, 2)
		t <- 0 #time at any point in the tree
		next.node<-4
		if(keepAppearances) {
			app<-cbind(c(2, 3), c(0, 0), c(NA, NA))
		}
		branchChanges<-rbind(c(1, 0, 1, 1, 1), c(1, 0, 2, 1, 1))
		colnames(branchChanges)<-c("Lineage", "Time", "Type", "Beginning", "End")

		changeCounter<-1

		############
		repeat{
			if(taxaStop) if(sum(alive)>=taxaStop) break;
			if(sum(alive)==0) break;
			rateVector<-c(sRates[alive]*b, eRates[alive]*d)
			dt<-rexp(1, sum(rateVector));
			makeChange<-F
			if(changeCounter<=length(timesToChange)) {
				if((t+dt)>timesToChange[changeCounter]) {
					makeChange=T
					}
			} 
			
			if(makeChange) {
				t<-timesToChange[changeCounter]
			} else t<-t+dt;
			
	  		if(timeStop) if(t>=timeStop) {
				t<-timeStop;
				break;
	  		}
	  		probVector<-cumsum(rateVector/sum(rateVector))

	  		r<-runif(1)
	  		event<-min(which(probVector>r))
	  	
	  		if(makeChange) {
	  			random_lineage <- sample(1:sum(alive))[1]
            	et<-1
            	or<-sRates[alive][random_lineage]
            	sRates[alive][random_lineage]<-newS[changeCounter]
            	nr<-sRates[alive][random_lineage]
            	branchChanges<-rbind(branchChanges, c(edge[alive,2][random_lineage], t, et, or, nr))

				et<-2
				or<-eRates[alive][random_lineage]
				eRates[alive][random_lineage]<-newE[changeCounter]
				nr<-eRates[alive][random_lineage]
            	branchChanges<-rbind(branchChanges, c(edge[alive,2][random_lineage], t, et, or, nr))
            	changeCounter<-changeCounter+1

	  		} else if(event<=sum(alive)) {#this creates a bifucation in the tree
	        	random_lineage <- event
				e<-matrix(edge[alive,], ncol=2)
				parent<-e[random_lineage,2]
				
				ps<-sRates[alive][random_lineage]
				pe<-eRates[alive][random_lineage]
				sRates<-c(sRates, ps, ps)
				eRates<-c(eRates, pe, pe)
				
				alive[alive][random_lineage]<-FALSE
				edge<-rbind(edge, c(parent, next.node), c(parent, next.node+1))
				if(keepAppearances) {
					app<-rbind(app, c(next.node+1, t, NA))
					app[app[,1]==parent,1]<-next.node
				}
				next.node<-next.node+2

				
				alive<-c(alive, TRUE, TRUE)
				stem.depth<-c(stem.depth, t, t)
				x<-which(edge[,2]==parent)
				edge.length[x]<-t-stem.depth[x]
				edge.length<-c(edge.length, NA, NA)

			
            } else if(event<=sum(alive)*2) {#This terminates one of the current lineages on the tree
            	random_lineage <- event-sum(alive)
		    	edge.length[alive][random_lineage]<-t-stem.depth[alive][random_lineage];
          	    if(keepAppearances) {
          	    	killed<-edge[alive,2][random_lineage]
					app[app[,1]==killed,3]<-t
				}
				alive[alive][random_lineage]<-FALSE

            } else {
            	
            }
      	}
      
		if(return.all.extinct==T | sum(alive)>1) break;
	}
	
	edge.length[alive]<-t-stem.depth[alive]
	

	tips<-!(edge[,2] %in% edge[,1])
	tt<-edge[tips,2]
	
	it<-as.numeric(levels(as.factor(edge[,1])))
	
	tCount<-1
	iCount<-1
	
	newEdge<-matrix(nrow=dim(edge)[1], ncol=dim(edge)[2])
	newApp<-matrix(nrow=dim(app)[1], ncol=dim(app)[2])
	newBC<-matrix(nrow=dim(branchChanges)[1], ncol=dim(branchChanges)[2])

	for(i in 1:max(edge)) {
		if(i %in% tt) {
			rr<-tCount
			tCount<-tCount+1
		} else {
			rr<-iCount+sum(tips)
			iCount<-iCount+1		
		}
		ok<-edge==i
		newEdge[ok]<-rr
		newApp[which(app[,1]==i), 1]<-rr
		newBC[which(branchChanges[,1]==i), 1]<-rr


	}
			
	edge<-newEdge
	app[,1]<-newApp[,1]
	branchChanges[,1]<-newBC[,1]

	print(edge)
	print(branchChanges)

	tip.label<-edge[tips,2]
		
	
    	mode(tip.label) <- "character"
    	obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label, Nnode=sum(alive))
    	class(obj) <- "phylo"


    
    if(keepAppearances) {
    	return(list(tree=obj, app=app, t=t, branchChanges=branchChanges))	
    } else {
    	return(obj)
    }
}


collapseTreeForMedusa<-function(phy, timeFraction=0.5, endClades=0) {
  	phenotype<-rep(1, length(phy$tip.label))
	names(phenotype)<-phy$tip.label
	rootDepth<-max(branching.times(phy))
	threshold<-rootDepth*timeFraction

	if(timeFraction==0 & endClades==0)
		return(NULL)

	while(1) {
		bt<-branching.times(phy)
		
		if(timeFraction!=0) {
			ct<-which(bt<threshold)
			if(length(ct)==0) break;
			nodeOrder<-order(bt[ct], decreasing=T)
			collapseNode<-names(bt[ct][nodeOrder])[1]
		} else if(endClades!=0) {
			if(length(bt)<=(endClades-1)) break;
			nodeOrder<-order(bt)
			collapseNode<-names(bt[nodeOrder])[1]

		}
		ns<-node.leaves(phy, collapseNode)
		whichTip<-which(phy$tip.label==ns[1])
		whichCollapse<-which(phy$tip.label %in% ns[-1])
		phenotype[whichTip]<-phenotype[whichTip]+sum(phenotype[whichCollapse])
		names(phenotype)[whichTip]<-paste(ns, collapse=".")

		phenotype<-phenotype[-whichCollapse]

		phy$tip.label[whichTip]<-paste(ns, collapse=".")
		for(j in 2:length(ns))
			phy<-drop.tip(phy, ns[j])
		
	}

	phy$phenotype<-phenotype
	phy
}

summaryChanges<-function(phy, branchChanges, ...) {
	for(i in 1:nrow(branchChanges)){
		nl<-node.leaves(phy, branchChanges[i,1])
		cat("Lineage: ", branchChanges[i,1], "\n")
		cat("Descendants:\n")
		cat(nl, "\n")
		
		}
		
	colors<-rep(1, length(phy$edge.length))
	oo<-as.numeric(as.factor(branchChanges[,1]))+1
	thickness<-rep(1, length(phy$edge.length))
	for(i in 1:nrow(branchChanges)) {
		ee<-match(branchChanges[i,1], phy$edge[,2])	
		if(!is.na(ee)) {
			colors[ee]<-oo[i]
			thickness[ee]<-thickness[ee]+1
		}
	}
	
	
	plot.phylo(phy, edge.color=colors, edge.width=thickness, ...)
		
	}
	
collapseTreeAndFossils<-function(data, timeFraction=0.5, endClades=0) {
	
	fullTree<-data$tree
	allC<-getAllClades(fullTree)
	allCSize<-numeric(length(allC))
	for(i in 1:length(allC)) allCSize[i]<-length(allC[[i]])

	survTree<-prune.extinct.taxa(fullTree)
	allSurvivors<-survTree$tip.label
	
	meduTree<-collapseTreeForMedusa(survTree, timeFraction, endClades)
	
	fullClade<-list()
	Nl<-numeric(length(phy$tip.label))
	Ne<-numeric(length(phy$tip.label))
	
	fossilDtt<-list()
	for(i in 1:length(phy$tip.label)) {
		tip<-meduTree$tip.label[i]
		clade<-strsplit(tip, split="\\.")[[1]]
		notInClade<-!(allSurvivors %in% clade)
		good<-rep(F, length(allC))
		for(j in 1:length(allC)) {
			ok1<-sum(!(clade %in% allC[[j]]))
			ok2<-sum(allSurvivors[notInClade] %in% allC[[j]])
			if(ok1==0 & ok2==0) good[j]<-T
		}
		mm<-max(allCSize[good])
		fullClade[[i]]<-allC[[which(good & allCSize==mm)]]
		Nl[i]<-length(clade)
		Ne[i]<-length(fullClade[[i]])-Nl[i]
		
		fossils<-!(fullClade[[i]] %in% allSurvivors)
		m<-match(fullClade[[i]][fossils], data$app[,1])
		aa<-data$app[m,]
		fossilDtt[[i]]<-getFossilDtt(aa)
	}	
	
	bc<-data$branchChanges
	newbc<-bc
	newbc<-cbind(newbc, rep(0, nrow(bc)))
	colnames(newbc)<-c(colnames(bc), "isInTriangle")
	
	allM<-getAllClades(meduTree)
	for(i in 1:nrow(branchChanges)) {
		t1<-node.leaves(fullTree, bc[i,1])
		ok<-rep(F, length(allM))
		size<-numeric(length(allM))
		for(j in 1:length(allM)) {
			tt<-strsplit(allM[[j]], split="\\.")
			taxa<-tt[[1]]
			if(length(tt)>1) for(k in 2:length(tt)) taxa<-c(taxa, tt[[k]])
			size[j]<-length(taxa)
			if(sum(!(taxa %in% t1))==0) ok[j]<-T
			}
		if(sum(ok)==0) { #The change is inside one of the triangles
		  newbc[i,6]<-1
		  ok[]<-F
		  for(j in 1:length(allM)) {
			tt<-strsplit(allM[[j]], split="\\.")
			taxa<-tt[[i]]
			if(length(tt)>1) for(k in 2:length(tt)) taxa<-c(taxa, tt[[k]])
			size[j]<-length(taxa)
			if(sum(taxa %in% t1)!=0) ok[j]<-T
		  }
		}	
		large<-max(size[ok])
		newCladeNumber<-which(size==max(size[ok]) & ok)
		newbc[i,1]<-newCladeNumber
		
	}
	
}

getFossilDtt<-function(dd) {
	times<-sort(as.numeric(dd[,2:3]))
	nd<-numeric(length(times))
	for(i in 1:length(times)) 
		nd[i]<-sum(dd[,2]<=times[i] & dd[,3]>times[i])
	cbind(times, nd)
	}

getAllClades<-function(phy) 
{
	res<-list()
	count=1
	ntip<-length(phy$tip.label)
	for(i in 1:(phy$Nnode+ntip)) {
		l<-node.leaves(phy, i)
		res[[i]]<-l	
	}
	res
}


prune.extinct.taxa<-function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    phy2 <- phy
    phy <- new2old.phylo(phy)
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(as.numeric(phy$edge[, 2]) == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    depth <- max(xx)
    offset <- depth - xx[names(xx) > 0]
    drops <- phy$tip.label[offset > tol]
    if (length(drops) >= (nb.tip - 1)) 
        return(NULL)
    if (length(drops) == 0) 
        return(phy2)
    res <- drop.tip(phy2, drops)
    res
}


#' runFossilMEDUSA: MEDUSA including fossil diversity data
#'
#' This function runs MEDUSA but includes data on fossil diversities through time. Requires 
#' phylogeny 'phy' and (optional) 'richness' (containing extant and potentially fossil
#' counts); 'fossil.richness' optionally passed in separately. If no richness information
#' is provided then it is assumed tips represent single species with comprehensive
#' sampling. Returns a list of models with 0, 1, 2, ..., 'model.limit' partitions.
#' 'model.limit' can be supplied by user, but is overruled if large enough (relative
#' to tree size) that AIC correction factor becomes undefined (different for pure birth vs.
#' birth-death models). 'phy' is assumed to be ultrametric and rescaled (i.e. calibrated) to time.
#'
#' @param phy An object of class 'phylo'
#' @param richness Richness data for extant species
#' @param fossil.richness Fossil richness data
#' @param est.extinction Estimate extinction rates or not?
#' @param fossil.minimum Treat fossils as minimum diversity estimates? If F, numbers are treated as true estimates
#' @param model.limit Maximum number of split models to consider
#' @param epsilon.value Input fixed epsilon value; if NULL, then epsilon is estimated
#' @param mc Use multicore or not 
#' @param ... other arguments; currently unused  
#' @return List with elements
#'  $models, which contains:
#'     $par: i x 2 matrix (for birth-death; i x 1 for pure-birth); the jth row contains 
#'       the speciation and (optional) extinction rate for the jth rate class
#'     $lnLik.part: vector of length i; the jth element is the partial log
#'       likelihood value due to the jth rate class
#'     $lnLik: = sum(lnLik.part); the overall log-likelihood of this model
#'     $split.at: the i+1 locations of splits.  The first element is the root node (i.e. background rate).
#'  $phy: the tree, which may have been manipulated (renamed, pruned) in the analysis
#'  $z: a matrix listing branch times and richnesses; for summarizing results
#'  $anc: a list of all ancestors; for summarizing results
#'  $model.summary: data.frame listing considered models, likelihoods, AIC, etc.
#'
#' The species richness information is assumed to have columns "taxon" and "n.taxa"; 
#' "taxon" must match with a tip.label in the phylogeny "phy". May also include 'exemplar' column,
#' used for incompletely-sampled clades which requiring pruning.
#' 

## *NOTE* Fossil observation currently limited to once per pendant edge.
##   This will be generalized (or at least extended) eventually...

## TO-DO:
 # 2. set epsilon to some value (for all piecewise models), estimate r (num.par = 2*num.models - 1)
 # 6. solve optimzation problem - pass in vector of values  - WORKING
 # 11. for birth-death models, check if break is at tip; if so, only estimate r (factor into AIC)

runFossilMEDUSA <- function(phy, richness=NULL, fossil.richness=NULL, est.extinction=T, fossil.minimum=F, model.limit=20, epsilon.value=NULL, mc=F, ...)
{
	if (is.null(richness))  # Assume tree represents single species tips and is completely sampled
	{
		richness <- data.frame(taxon=phy$tip.label, n.taxa=1, n.fossils=NA, f.time=NA)
	} else {
## Before determining model.limit, prune tree as necessary (from 'taxon' information in 'richness')
		phyData <- prune.tree.merge.data(phy, richness, fossil.richness, fossil.minimum)
		phy <- phyData$phy
		richness <- phyData$richness
	}
	
## Limit on number of piecewise models fitted; based on tree size, AICc correction factor, and flavour of model fitted (i.e. # parameters estimated)
	model.limit <- get.max.model.limit(richness, model.limit, est.extinction)

## Store pertinent information: branch times, richness, ancestors
	obj <- make.cache.medusa(phy, richness)
	
## Keep track of all nodes, internal and pendant (for keeping track of breakpoints)
	pend.nodes <- seq_len(length(phy$tip.label))   # Ideally, calculate pendant splits just once, keep track through various models
	int.nodes <- (length(phy$tip.label)+1):max(phy$edge)
	all.nodes <- c(pend.nodes, int.nodes)
	
	anc <- obj$anc
	z <- obj$z
	z.orig <- z # Save for summarizing models
	
## Pre-fit pendant edges so these values need not be re(re)calculated; amounts to ~25% of calculations
	tips <- list()
	cat("Optimizing parameters for pendant edges...\n\n")
	tips <- lapply(pend.nodes, medusa.ml.tips, z, anc, est.extinction, epsilon.value, fossil.minimum)
	
## 'fit' holds current results; useful for initializing subsequent models
	fit <- medusa.ml.initial(z, est.extinction, epsilon.value, fossil.minimum)
	models <- list(fit)
	
	
	cat("Step 1 (of ", model.limit, "): best LH = ", models[[1]]$lnLik, "\n", sep="")

	for (i in seq_len(model.limit-1))
	{
		node.list <- all.nodes[-fit$split.at]
		if (mc == T)  # multicore (i.e. multithreaded) processing. Only works in Terminal for some reason.
		{
			res <- mclapply(node.list, medusa.ml.update, z, anc, fit, tips, est.extinction, epsilon.value, fossil.minimum)
		} else {
			res <- lapply(node.list, medusa.ml.update, z, anc, fit, tips, est.extinction, epsilon.value, fossil.minimum)
		}
		best <- which.min(unlist(lapply(res, "[[", "AICc")))  # Better for dealing with models of different dimensions
		models <- c(models, res[best])
		z <- medusa.split(node.list[best], z, anc)$z
		fit <- res[[best]]   # keep track of '$split.at' i.e. nodes already considered
		
		cat("Step ", i+1, " (of ", model.limit, "): best LH = ", models[[i+1]]$lnLik, "; break at node ", models[[i+1]]$split.at[i+1],"\n", sep="")
	}
	
	model.summary <- calculate.model.fit.summary(models, phy)
	cat("\n", "Model fit summary:", "\n\n", sep="")
	print(model.summary)
	
#	results <- list(models=models, phy=phy, z=z.orig, anc=anc)
	results <- list(z=z.orig, anc=anc, models=models, phy=phy, model.summary=model.summary)
	
	return(results)
}



# Function to prune tree using 'richness' information, assumed to have minimally two columns, "taxon" and "n.taxa"
  # Perhaps relax on these column names, may cause too many problems
# May also include 'exemplar' column; in that case, rename relevant tip.label before pruning.
# In addition, may have columns 'n.fossils' and 'f.time'; must match names exactly or ambiguous error generated.
 # Taxa without fossil counts should have NA in n.fossil and f.time columns.
# If 'fossil.richness' passed in separately, merge with 'richness'
prune.tree.merge.data <- function(phy, richness, fossil.richness=NULL, fossil.minimum=F)
{
# Rename exemplar taxa with taxon name in richness file
	if (!is.null(richness$exemplar))
	{
# Change relevant tip.labels in phy; individual 'exemplar' may be NA, use original tip.label.
# Ordering in richness file should NOT be assumed to match order of tip.labels
		i.na <- is.na(richness$exemplar)
		phy$tip.label[match(richness$exemplar[!i.na], phy$tip.label)] <- as.character(richness$taxon[!i.na])
	}
	
# Check names against 'richness' (already fixed above). May use 'exemplar', but *must* use 'taxon'
# Update 'n.fossils' and 'f.time' accordingly
# Again, assume ordering is random, but that taxa are consistent across richness files
# Annoying number of possible combinations of file format. Don't be so nice! Force them into one or few.
	if (!is.null(fossil.richness)) # Merge 'fossil.richness' with 'richness'
	{
		if (!is.null(fossil.richness$exemplar) && !is.null(richness$exemplar))
		{
			i.fossil <- match(fossil.richness$exemplar, richness$exemplar)
		} else if (!is.null(fossil.richness$exemplar) && !is.null(richness$taxon))
		{
			i.fossil <- match(fossil.richness$exemplar, richness$taxon)
# This would be odd, but allow for it; need to change tip labels to 'exemplar' from 'fossil.richness'
	# VERY low priority
		} else {
			i.fossil <- match(fossil.richness$taxon, richness$taxon)
		}
		
		if (length(i.fossil) != length(fossil.richness[,1])) {stop("Problem matching extant and fossil taxon information. Fucktard.")}
		
		for (i in 1:length(fossil.richness[,1]))
		{
			richness[i.fossil[i],"n.fossils"] <- fossil.richness[i,"n.fossils"]
			richness[i.fossil[i],"f.time"] <- fossil.richness[i,"f.time"]
		}
	}
# Just for consistency with downstream functions:
	if (is.null(richness$n.fossils)) {richness <- cbind(richness, n.fossils=NA, f.time=NA)}
	
# Check if fossil.minimum==T and if any n.fossil==1; delete and warn
	if (fossil.minimum == T)
	{
		if (!is.na(any(richness$n.fossils == 1)))  # is it possible to be uglier?!? make pretty later.
		{
			if(any(richness$n.fossils == 1))
			{
				drop.fossil <- which(richness$n.fossils == 1)
				drop.fossil.taxa <- as.character(richness$taxon[drop.fossil])
				
				richness$n.fossils[drop.fossil] <- NA
				richness$f.time[drop.fossil] <- NA
				
				cat("Warning: dropped fossil observations since cannot use fossil counts == 1 as minimums:\n")
				print(drop.fossil.taxa)
				cat("\n")
			}
		}
	}
	
	if (length(phy$tip.label) != length(richness[,1]))
	{
# Prune tree down to lineages with assigned richness
		temp <- richness[, "n.taxa"]
		names(temp) <- richness[, "taxon"]
		pruned <- treedata(phy, temp)		# geiger function calling ape (namecheck)
		prunedTree <- pruned$phy
# Check the tree
	#	plotNN(prunedTree)					# Node numbers (ape-style) plotted
		phy <- prunedTree
	}
	return(list(phy=phy, richness=richness))
}



## Original default was to fit 20 models (or less if the tree was small).
## Changing to a stop-criterion e.g. when k = n-1 (i.e. when denominator of AICc correction is undefined).
## k <- (3*i-1) # when BOTH birth and death are estimated
  ## This occurs when i = n/3
  ## if est.extinction=F, k <- (2*i-1); max i = n/2
## n <- (2*num.taxa - 1) == (2*length(richness[,1]) - 1) # i.e. total number of nodes in tree (internal + pendant)
  ## if present, consider each fossil count as one data point; i.e. n <- num.nodes + num.fossils
## Eventually use AICc itself as a stopping criterion (not urgent, as it is very fast).
get.max.model.limit <- function(richness, model.limit, est.extinction)
{
	samp.size <- sum(!is.na(richness[,"n.fossils"])) + (2*length(richness[,1]) - 1)
	
	if (est.extinction == T)
	{
		max.model.limit <- as.integer(samp.size/3)
		if (model.limit > max.model.limit) {model.limit <- max.model.limit}
		cat("Limiting consideration to ", model.limit, " piecewise BD models\n\n", sep="")
	} else {
		max.model.limit <- as.integer(samp.size/2)
		if (model.limit > max.model.limit) {model.limit <- max.model.limit}
		cat("Limiting consideration to ", model.limit, " piecewise pure birth models\n\n", sep="")
	}
	return(model.limit)
}



## The make.cache.medusa function is like the first half of the original splitEdgeMatrix().
## It works through and reorders the edges, then works out start and end times of these
## based on the phylogeny's branching times.
##
## Rather than return a matrix, I am returning a list so that a few
## additional bits can easily be returned.  
##
## In addition, every node's ancestors are also calculated.  The element 'anc' is a list.
## $anc[i] contains the indices within $edge, $t.start, etc., of all ancestors of node 'i'
## (in ape node numbering format).
make.cache.medusa <- function(phy, richness)
{
	if (is.null(richness)) {stop("This should not be possible: You, sir, are a fucktard.")}
	
	n.tips <- length(phy$tip.label)
	n.edges <- nrow(phy$edge)
# Alternatively, for a rooted tree:
#	n.edges <- (2* n.tips) - 2
# unless... will polytomies ever be considered?	
	n.int <- n.edges - n.tips
	
## Ape numbers the tips first
	i.int <- seq_len(n.int)
	interior <- phy$edge[,2] %in% phy$edge[,1]
	bt <- branching.times(phy)
	
# Broken up to more easily parse fossil information; consider only internal edges first
	edges.int <- phy$edge[interior,]
	colnames(edges.int) <- c("anc", "dec")
	
	t.0 <- bt[match(edges.int[,1], (n.tips+1):max(edges.int))]
	t.1 <- c(t.0[i.int] - phy$edge.length[interior])
# Why use 'edge.length' here and not simply edges.int[,2]?
# Seems to give exact numbers, but R says == is false for a few entries; rounding error?
#	t.1a <- bt[match(edges.int[,2], (n.tips+1):max(edges.int))]
	
	z.internal <- cbind(edges.int, t.0, t.f=rep(NA, n.int), t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.int), n.f=rep(NA, n.int), n.t=rep(NA, n.int))
	
# Now, pendant edges	
	edges.pendant <- phy$edge[match(seq_len(n.tips), phy$edge[,2]),]
	colnames(edges.pendant) <- c("anc", "dec")

	t.0 <- bt[match(edges.pendant[,1], (n.tips+1):max(edges.pendant))]
	t.1 <- rep(0, n.tips)

# Ensure richnesses are ordered appropriately
	ext.richness <- richness$n.taxa[match(phy$tip.label, richness$taxon)]
	fossil.richness <- richness$n.fossils[match(phy$tip.label, richness$taxon)]
	fossil.times <- richness$f.time[match(phy$tip.label, richness$taxon)]
		
	z.pendant <- cbind(edges.pendant, t.0, t.f=fossil.times, t.1, t.len=t.0 - t.1,
		n.0=rep(1, n.tips), n.f=fossil.richness, n.t=ext.richness)

# Used for identifying ancestral nodes below i.e. tracking breakpoints
	all.edges <- rbind(edges.int, edges.pendant)
	
	z <- rbind(z.internal, z.pendant)
	z <- cbind(z,partition=rep(1, n.edges))
	rownames(z) <- NULL

	list(z=z, anc=lapply(seq_len(max(all.edges)), ancestors.idx, all.edges)) # And, we're good to go...
}



## This generates the indices of all ancestors of a node, using ape's edge matrix.
ancestors <- function(node.list, all.edges)
{
	ans <- node.list
	repeat
	{
		node.list <- all.edges[all.edges[,1] %in% node.list,2]
		if (length(node.list) > 0) {ans <- c(ans, node.list)} else {break}
	}
	unlist(ans)
}

## The function 'ancestors' returns the indices of all ancestors within the edge matrix.
ancestors.idx <- function(node.list, all.edges)
{
	which(all.edges[,1] == node.list | all.edges[,2] %in% ancestors(node.list, all.edges))
}



# Only used for base model
medusa.ml.initial <- function(z, est.extinction, epsilon.value, fossil.minimum)
{
	rootnode <- min(z[is.na(z[,"n.t"]),1])
	obj <- medusa.ml.fit.partition(1, z, est.extinction, epsilon.value, fossil.minimum=fossil.minimum)
	
	aic.fit <- calculate.AIC(fit=obj, z, est.extinction, epsilon.value)
	
	if (est.extinction)
	{
		list(par=matrix(obj$par, nrow=1, dimnames=list(NULL,c("r", "epsilon"))), lnLik.part=obj$lnLik, 
		   lnLik=obj$lnLik, split.at=rootnode, AIC= aic.fit[1], AICc=aic.fit[2], num.par=aic.fit[3])
	} else {
		list(par=matrix(obj$par, nrow=1,dimnames=list(NULL,"r")), lnLik.part=obj$lnLik, lnLik=obj$lnLik,
		   split.at=rootnode, AIC= aic.fit[1], AICc=aic.fit[2], num.par=aic.fit[3])
	}
}


# Pre-fit values for pendant edges; DON'T recalculate later; should account for ~25% of all calculations
medusa.ml.tips <- function(node, z, anc, est.extinction, epsilon.value, fossil.minimum)
{
	obj <- medusa.split(node, z, anc)
	z <- obj$z
# Partition '2' represents the pendant edge
	fitted <- medusa.ml.fit.partition(2, z, est.extinction, epsilon.value, fossil.minimum=fossil.minimum)
	
	return(fitted)
}



# 'fit' contains parameter value(s) from previous model, used to initialize subsequent model
# Pass in pre-fitted values for pendant edges; DON'T recalculate 
medusa.ml.update <- function(node, z, anc, fit, tips, est.extinction, epsilon.value, fossil.minimum)
{
	obj <- medusa.split(node, z, anc)
	z <- obj$z
	aff <- obj$affected
	if (length(aff) > 2) stop("This should not happen, fucktard.")
	
	op <- fit$par
	sp <- op[aff[1],] # Use previously fit parameter values from clade that is currently being split
# expecting: medusa.ml.fit.partition(partition, z, est.extinction, epsilon.value, sp=c(0.1, 0.05), fossil.minimum)
	fit1 <- medusa.ml.fit.partition(aff[1], z, est.extinction, epsilon.value, sp, fossil.minimum)
	fit2 <- 0
	
	if (node <= (length(z[,1])/2 + 1))  # Number of tips
	{
		fit2 <- tips[[node]]
		} else {
		fit2 <- medusa.ml.fit.partition(aff[2], z, est.extinction, epsilon.value, sp, fossil.minimum)
	}
	
	op[aff[1],] <- fit1$par # Replace parameters with new values for diminished clade
	fit$par <- rbind(op, fit2$par)
	fit$lnLik.part[aff] <- c(fit1$lnLik, fit2$lnLik) # Replace parameters with new values for diminished clade
	fit$split.at <- c(fit$split.at, node)
	fit$lnLik <- sum(fit$lnLik.part)
	
	aic.fit <- calculate.AIC(fit, z, est.extinction, epsilon.value)
	
	fit$AIC <- aic.fit[1]
	fit$AICc <- aic.fit[2]
	fit$num.par <- aic.fit[3]
	
	return(fit)
}



## Split the edge matrix 'z' by adding a partition rooted at node 'node'.
##   Note: in original MEDUSA parlance, this is cutAtStem=T.
## The list 'anc' is a list of ancestors (see make.cache.medusa, above).
##   It returns a list with elements:
##   z: new medusa matrix, with the new partition added
##   affected: indices of the partitions affected by the split (n == 2).
# Add cutAtStem = F option. Allow both within a single analysis?!?
medusa.split <- function(node, z, anc)
{
	part <- z[,"partition"]
	base <- min(part[z[,1] == node | z[,2] == node])
	tag <- max(part) + 1

	i <- anc[[node]]
	idx <- i[part[i] == base]
	z[idx,"partition"] <- tag

	list(z=z, affected=c(unique(part[idx]), tag))
}



# Where optimization *should* to happen...
# sp = initializing values for r & epsilon
# Default values should never be used (except for first model), as the values from the previous model are passed in.
#medusa chokes on really big timetreees because r is close to one. 
###The thing with your tree is that the pendant edges are HUGE (the smallest is 278 MY), such that, for any large-ish value of r we might be considering, B is ~1. log(1 - B) is thus -Inf and R dies. Original MEDUSA is dying too. This should be a pervasive problem with old trees (thanks for bringing it to my attention). The obvious fix that comes to mind is rescaling the tree so that the "unit" of time is, say, 10 MY rather than 1 MY. I'll play with the tree and see what I get. I have been able to get it to run by specifying smaller starting values for r (the maximum value that doesn't crash is 0.04). To change this, go to the following function and change sp (starting parameter values) from c(0.1, 0.05) to c(0.04, 0.05).

### default medusa.ml.fit.partition <- function(partition, z, est.extinction, epsilon.value, sp=c(0.1, 0.05), fossil.minimum)

#bigtrees option (may cause optim to become stuck locally)
medusa.ml.fit.partition <- function(partition, z, est.extinction, epsilon.value, sp=c(0.04, 0.05), fossil.minimum)
{
# Construct likelihood function:
	lik <- make.lik.medusa.part(z[z[,"partition"] == partition,,drop=FALSE], est.extinction, epsilon.value, fossil.minimum)
	
# *** Fix optimization problem ***
	if (est.extinction) # Give a range of values to "navigate the banana"
	{
		fit <- optim(fn=lik, par=sp, method="N", control=list(fnscale=-1))
		list(par=fit$par, lnLik=fit$value)
	} else {
		fit <- optimize(fn=lik, interval=c(0, 1))  # check on a valid range
		list(par=fit$minimum, lnLik=fit$objective)
	}
}



## make.lik.medusa.part: generate a likelihood function for a single partition.
##
## In the pendant calculations, the variables 'A' and 'B' are the A and B terms
## in Foote et al. Science: 283 1310-1314, p. 1313, defined as:
## 
##   A: probability of extinction of one lineage over time 't'
##   B: A * lambda / mu OR A / epsilon
##
## Where there is a single lineage at time 0 (a==1), the calculation is
##   log(1 - A) + log(1 - B) + (n-1)*log(B)
## but this is conditioned on survival by dividing by (1-A)
## (subtracting log(1-A) on a log scale) which cancels to give:
##   log(1 - B) + (n-1)*log(B)
## when pendant richness==1 the calculation is even simpler:
##   log(1 - B)
##
## When a > 1 (more than one starting species; will only happen coming off of a fossil observation),
## the calculations follow Foote et al. 1999, put into the function 'foote.fossil.pendant'.
## 'est.extinction' determines whether yule model is assumed
## 'epsilon.value', if != NULL, fixes epsilon for all partitions, estimates r; not yet implemented
## 'fossil.minimum', if == F, treat fossils as exact counts
## Only one 'partition' is passed in at a time from a proposed split.
make.lik.medusa.part <- function(partition, est.extinction, epsilon.value, fossil.minimum)
{
# Handle internal, fossil, and pendant edges separately
	is.int <- is.na(partition[,"n.t"])
	is.fossil <- !is.na(partition[,"n.f"]) # no 'n.f'? no fossil.
	is.pend <- !is.int == !is.fossil
	
	n.int <- sum(is.int)
	n.fossil <- sum(is.fossil)
	n.pend <- sum(is.pend)
	
	if (n.int + n.fossil + n.pend != length(partition[,1])) stop("You messed up, yo. Fucktard.")
	
## Internal and pendant calculations differ; split'em up
	int  <- partition[is.int,,drop=FALSE]
	foss <- partition[is.fossil,,drop=FALSE]
	pend <- partition[is.pend,,drop=FALSE]
	
	sum.int.t.len <- sum(int[,"t.len"])  # Simply sum all internal edges
	int.t.0 <- int[,"t.0"]
	
# 'n.0' = Foote's 'a', initial diversity; 'nt' = Foote's 'n', final diversity
	pend.n.0 <- pend[,"n.0"] # Foote's 'a': initial diversity
	pend.n.t <- pend[,"n.t"] # Foote's 'n': final diversity
	pend.t.len <- pend[,"t.len"]
	
	foss.n.0 <- foss[,"n.0"]   # Will always be 1
	foss.n.f <- foss[,"n.f"]
	foss.n.t <- foss[,"n.t"]
	t.len.int.to.foss <- foss[,"t.0"] - foss[,"t.f"]
	t.len.foss.to.pend <- foss[,"t.f"] - foss[,"t.1"]
	
# User may pass in epsilon; don't change it, just estimate r
	f <- function(pars)
	{
		if (est.extinction == T)
		{
			r <- pars[1]
			epsilon <- pars[2]
	#		bd <- get.b.d(r, epsilon)   # Not necessary
	#		b <- bd$b
	#		d <- bd$d
			
# Previous version where b & d were estimated:
	#		b <- pars[1]
	#		d <- pars[2]
	#		r <- abs(b - d)
	#		epsilon <- d / b
			
			if (r < 0 | epsilon <= 0 | epsilon >= 1) {return(-Inf)}
		} else {
## Implement pure birth; not currently working. Not immediately urgent.
			r <- pars[1]
			d <- 0
			b <- r
			epsilon <- 0
			
			if (r < 0) {return(-Inf)}
		}
		
		if (n.int == 0) {l.int <- 0} else {
## Likelihood of internal edges from Rabosky et al. (2007) equation (2.3):
			l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - (epsilon * exp(-r * int.t.0))))
# How to do this when pure birth? Is it cool as it is? Ask Luke...
# Under Yule, should be: l.int <- n.int * log(r) - r * sum.int.t.len

		}
		
		if (n.pend == 0) {l.pend <- 0} else {
## Calculations are from the following:
## Rabosky et al. 2007. Proc. Roy. Soc. 274: 2915-2923.
## Foote et al. 1999. Science. 283: 1310-1314
## Raup. 1985. Paleobiology 11: 42-52 [Foote et al. correct the equation [A18] where a > 1]
## Bailey. 1964. The Elements Of Stochastic Processes, With Applications To The Natural Sciences
## Kendall. 1948. Ann. Math. Stat. 19: 1–15.
##
## A = probability of extinction of one lineage over time 't'
## B = A * (lambda/mu)
##
## When there is a single lineage at time 0 (a = 1), the calculation is
##   log(1 - A) + log(1 - B) + (n - 1)*log(B)
## but this is conditioned on survival by dividing by (1-A)
## (subtracting log(1 - A) on a log scale) which cancels to give:
##   log(1 - B) + (n - 1)*log(B)
##      - for n.t == 1, reduces further to log(1-B)
##
## A = mu*(exp((lambda - mu)*t) - 1)) / (lambda*exp((lambda - mu)*t) - mu)
##  let r = (lambda - mu); ert = exp((lambda - mu)*t)
## A = mu*(ert - 1)/(lambda*ert - mu)
##
## B = A * (lambda/mu)
##   = [mu*(ert - 1)/(lambda*ert - mu)] * (lambda/mu)
##   = (lambda*(ert - 1))/(lambda*ert - mu)
##   = (lambda*(ert - 1))/(lambda(ert - mu/lambda))
##   = (ert - 1) / (ert - epsilon)

			if (all(pend.n.0 == 1))   # All pendant nodes begin with richness '1'; calculations simple.
			{
				i.pend.n.t.1 <- which(pend.n.t == 1)   # calculations even simpler: log(1-B)
				i.pend.n.t.n1 <- which(pend.n.t != 1)
				
				ert <- exp(r * pend.t.len)
				B <- (ert - 1) / (ert - epsilon) # Equivalently: B <- (bert - b) / (bert - d)
				
		# Under Yule should be: B <- (ert - 1) / ert
				
	# separate out calculations; likely not any faster (slower?)
#				l.pend.1 <- sum(log(1 - B[i.pend.n.t.1]))
#				l.pend.n1 <- sum(log(1 - B[i.pend.n.t.n1]) + (pend.n.t[i.pend.n.t.n1] - 1)*log(B[i.pend.n.t.n1]))
#				l.pend <- l.pend.1 + l.pend.n1
				
				l.pend <- sum(log(1 - B[i.pend.n.t.1])) + 
				   sum(log(1 - B[i.pend.n.t.n1]) + (pend.n.t[i.pend.n.t.n1] - 1)*log(B[i.pend.n.t.n1]))
				
#				l.pend <- sum(log(1 - B) + (pend.n.t - 1)*log(B)) # equivalent single formula version
			} else stop("Fucktard. How the hell can a pendant edge have n.0 > 1 ?!?")
		}
		
# Fossils!
		if (n.fossil == 0) {l.fossil <- 0} else {

# What is desired: P(n,t,N,T) = P(n,t,1)*P(N,T-t,n)
# - each term (internal node to fossil, and fossil to pendant node) is conditioned on survival (dividing by [1-(P(0,t,a)]^a)]

			if (fossil.minimum) # loop counts from 1 -> n-1 fossils; *** NEED TO SPEED THIS SHIT UP! ASK FITZJOHN ***
			{
## First consider path from internal node to fossil, conditioning on survival (i.e. divide by (1-A)^a, although 'a' is always 1 here):
##  sum(over x from 1 to n-1) [(1-B)*B^(x-1)]
				
				ert <- exp(r * t.len.int.to.foss)
				B.1 <- (ert - 1) / (ert - epsilon) # Equivalently: B.1 <- (bert - b) / (bert - d)
			# Under Yule should be: B <- (ert - 1) / ert
				
				sum.prob.int.fossil <- numeric(length(foss.n.f))
				for (i in 1:n.fossil)
				{
					tmp <- 0
					for (j in 1:foss.n.f[i])
					{
						tmp <- sum((B.1[i]^(j-1)), tmp)
					}
					sum.prob.int.fossil[i] <- (1 - B.1[i]) * tmp # (1-B) brought outside sum, as it is present in every term
				}
				
# Now consider path from fossil to pendant node; looping bits in 'foote.fossil.pendant.minimum'
				ert <- exp(r * t.len.foss.to.pend)
				A <- (ert - 1) / ((ert/epsilon) - 1) # Equivalently: A <- (d*(ert - 1)) / ((b*ert) - d)
				B.2 <- A / epsilon                   # Equivalently: B.2 <- A * (b/d)
				
			# Using raw probabilities:
				sum.prob.fossil.pend <- foote.fossil.pendant.minimum(foss.n.f, foss.n.t, A, B.2)
				l.fossil <- sum(log(1 - (sum.prob.int.fossil * sum.prob.fossil.pend)))
				
			# Using log-scale values:
#				sum.prob.fossil.pend <- foote.fossil.pendant.combined(foss.n.f, foss.n.t, A, B.2, fossil.minimum)
#				l.fossil <- sum(log(1 - (sum.prob.int.fossil * exp(sum.prob.fossil.pend))))
				
			} else { # No looping; pretty fast
				ert <- exp(r * t.len.int.to.foss) # internal node to fossil
				B.1 <- (ert - 1) / (ert - epsilon) # Equivalently: B.1 <- (bert - b) / (bert - d)
		# Under Yule should be: B <- (ert - 1) / ert
				
				ert <- exp(r * t.len.foss.to.pend)     # fossil to pendant node
				A <- (ert - 1) / ((ert/epsilon) - 1) # Equivalently: A <- (d*(ert - 1)) / ((b*ert) - d)
				B.2 <- A / epsilon                   # Equivalently: B.2 <- A * (b/d)
				
				l.fossil <- sum(log(1 - B.1) + (foss.n.f - 1)*log(B.1)) + 
				   foote.fossil.pendant.exact(foss.n.f, foss.n.t, A, B.2)
			}
		}
		l.int + l.pend + l.fossil
	}
}



## Calculates the probability of having *at least* 'a' fossils and ending with exactly 'n' extant species.
## Richness may increase or decrease.
## Should this be computed on a log-scale?
foote.fossil.pendant.minimum <- function(a, n, A, B)
{
#	lA <- log(A)
	A.exp.a <- A^a					# P(0,t,a) i.e. does not survive
#	lB <- log(B)
#	l1mA <- log(1 - A)
#	l1mA.exp.a <- log(1 - A.exp.a)	# P(lineage survives), used to condition upon, logscale
	r1mA.exp.a <- (1 - A.exp.a)		# P(lineage survives), used to condition upon, rawscale
	
#	l1mB <- log(1 - B)
#	l1mA1mB <- l1mA + l1mB
	
	sum.prob.fossil.pend <- numeric(length(a))
	
	for (i in 1:length(a)) # Loop over fossils
	{
		tmp <- numeric(a[i]-1)
		for (k in 1:(a[i]-1)) # For fossil i loop over from 1 to n[i]-1
		{
			ai <- k
			ni <- n[i]
			min.a.n <- pmin.int(ai, ni)
#			min.a.n <- pmin(ai, ni) # *Insanely* slow! pmin.int saves 70s on a 107s job! WTF?!?
	   		j <- seq_len(min.a.n)
			
	# *** Should this be put on a log scale to preserve precision? ***
	  # Tests with reasonable data yield the same answer either way
			foo <- choose(ai, j) * choose(ni - 1, j - 1) * A[i]^(ai - j) * ((1-A[i])*(1-B[i]))^j * B[i]^(ni - j)
			tmp[k] <- sum(foo)
	# Alternate log form:
	  # foo <- lchoose(ai, j) + lchoose(ni - 1, j - 1) + (ai - j)*lA[i] + j*l1mA1mB[i] + (ni - j)*lB[i]
	    }
	#  	tmp[i] <- sum(tmp)/(1-A[i])  # old conditioning; not appropriate for a > 1
			
# Conditioning on survival factors out of the sum:
		sum.prob.fossil.pend[i] <- sum(tmp)/(r1mA.exp.a[i])
	# Alternate log form:
	  #	foo <- (logspace_sum(tmp))
	  #	sum.prob.fossil.pend[i] <- foo - l1mA.exp.a[i]
	}

	return(sum.prob.fossil.pend)
}


## Calculates the probability of having *exactly* 'a' fossils and ending with *exactly* 'n' extant species.
## Use only when 'a', starting richness, > 1 (i.e. if coming off a fossil count).
## For 'off-fossil' calculations, richness may increase or decrease.
foote.fossil.pendant.exact <- function(a, n, A, B)
{
	lA <- log(A)
	A.exp.a <- A^a					# P(0,t,a) i.e. does not survive
	lB <- log(B)
	l1mA <- log(1 - A)
	l1mA.exp.a <- log(1 - A.exp.a)	# P(lineage survives), used to condition upon, logscale
	r1mA.exp.a <- (1 - A.exp.a)		# P(lineage survives), used to condition upon, rawscale
	
	l1mB <- log(1 - B)
	l1mA1mB <- l1mA + l1mB
	
#	min.a.n <- pmin(a, n) # *Insanely* slow!
	min.a.n <- pmin.int(a, n)
	
	l.pend.fossil <- numeric(length(min.a.n))
	
	for (i in 1:length(a))
	{
		ji <- seq_len(min.a.n[i])
		ai <- a[i]
		ni <- n[i]
		
#		tmp <- lchoose(ai, ji) + lchoose(ni - 1, ji - 1) + (ai - ji)*lA[i] + ji*l1mA1mB[i] + (ni - ji)*lB[i]
		
# Non-logspace version much faster
		foo <- choose(ai, ji) * choose(ni - 1, ji - 1) * A[i]^(ai - ji) * ((1-A[i])*(1-B[i]))^ji * B[i]^(ni - ji)
		tmp <- sum(foo)
		
## Conditioning on survival factors out of the sum:
		l.pend.fossil[i] <- sum(tmp)/(r1mA.exp.a[i])
		
#		l.pend.fossil[i] <- logspace_sum(tmp) - l1mA.exp.a[i] # logspace_sum: calculate log(x+y) from log(x) and log(y)
	#	l.pend.fossil[i] <- logspace_sum(tmp) - l1mA[i]  # old conditioning; not appropriate for a > 1
	}

#	sum(l.pend.fossil)
	sum(log(l.pend.fossil))
}



# 'fit' contains '$par' and '$lnlik'
calculate.AIC <- function(fit, z, est.extinction, epsilon.value)
{
## Sample size taken (for now) as the total num.nodes in the tree (internal + pendant) plus num.fossil observations
  # num.nodes = (2*length(phy$tip.label) - 1) == (2*length(richness[,1]) - 1) == length(z[,1]) + 1
	n <- (length(z[,1]) + 1) + sum(!is.na(z[,"n.f"]))
	
## Number of parameters estimated depends on values of est.extinction and epsilon.value
  # Includes both formal parameters AND number of breaks. Note: first model does not involve a break.
## Models where all parameters are estimated:
  # 2 parameters for base model (no breakpoint) + 3 parameters (r, eps, breakpoint) for each subsequent model
## Models where only one parameter is estimated:
  # 1 parameter for base model (no breakpoint) + 2 parameters (r, breakpoint) for each subsequent model
	
	if (length(fit$par) < 3) # i.e. base model
	{
		num.models <- 1
	} else {
		num.models <- length(fit$par[,1])
	}
	
	if (est.extinction == F || !is.null(epsilon.value)) # only one parameter estimated (plus breakpoint)
	{
		k <- 1 + (2 * (num.models - 1))
	} else {
		k <- 2 + (3 * (num.models - 1))
	}
	
	lnLik <- fit$lnLik
	
	AIC <- (-2 * lnLik) + (2*k)
	AICc <- AIC + 2*k*(k+1)/(n-k-1)
	
	aic.fit <- c(AIC, AICc, k)
	return(aic.fit)
}


## Prints out a table of likelihoods, parameters, AIC scores, and AIC weights (deltaAICs are also available, if desired)
## Add functionality to print out coloured tree as well
calculate.model.fit.summary <- function (models, phy, plot.fig=T)
{
	tmp <- matrix(nrow=(length(models)), ncol=6)
	colnames(tmp) <- c("N.Models", "Break.Node", "Ln.Lik", "N.Param", "AIC", "AICc")
	
	w.AIC <- numeric(length(models))
	w.AICc <- numeric(length(models))
	
	for (i in 1:length(tmp[,1]))
	{
		tmp[i,] <- c(i, as.integer(models[[i]]$split.at[i]), models[[i]]$lnLik, models[[i]]$num.par, models[[i]]$AIC, models[[i]]$AICc)
	}
	
	all.res <- as.data.frame(tmp)
	all.res[1,2] <- NA # root node for base model
	
	w.AIC <- calculate.akaike.weights(all.res$AIC)
	w.AICc <- calculate.akaike.weights(all.res$AICc)
	
#	all.res <- cbind(all.res[,c(1:5)], w.AIC=w.AIC$w, AIC.best=w.AIC$best, AICc=all.res$AICc, w.AICc=w.AICc$w, AICc.best=w.AICc$best)
	all.res <- cbind(all.res[,c(1:5)], w.AIC=w.AIC$w, AICc=all.res$AICc, w.AICc=w.AICc$w)
	
#	options(digits=5)
	
	if (plot.fig){plot.AIC.models(all.res)}
	return(all.res)
}

calculate.akaike.weights <- function (AIC)
{
	best <- min(AIC)
	deltaAIC <- AIC-best
	sumDeltaAIC <- sum(exp(-0.5 * deltaAIC))
	w <- (exp(-0.5 * deltaAIC)/sumDeltaAIC)
	
	ind <- character(length(deltaAIC))
	for (i in 1:length(deltaAIC))
	{
		if(deltaAIC[i]==0) {ind[i] <- "*"}
	}
	
	results<-data.frame("AIC"=AIC,"dAIC"=deltaAIC,"w"=w,"best"=ind)
	
	return(results)
}

plot.AIC.parameters <- function (all.res)
{
	ylim <- c(min(all.res[,"AIC"],all.res[,"AICc"]), max(all.res[,"AIC"],all.res[,"AICc"]))
	plot(all.res[,"N.Param"],all.res[,"AICc"], xlab="Number of Parameters", ylab="Model Fit", ylim=ylim, type="l", col="blue")
	points(all.res[,"N.Param"],all.res[,"AICc"], col="blue", pch=21, bg="white")
	points(all.res[,"N.Param"],all.res[,"AIC"], col="black", type="l")
	points(all.res[,"N.Param"],all.res[,"AIC"], col="black", pch=21, bg="white")
	
	legend("topleft", c("AICc","AIC"), pch=21, pt.bg="white", lty=1, col=c("blue", "black"), inset = .05, cex=0.75, bty="n") # 'bottomright' also works
}

plot.AIC.models <- function (all.res)
{
	ylim <- c(min(all.res[,"AIC"],all.res[,"AICc"]), max(all.res[,"AIC"],all.res[,"AICc"]))
	plot(all.res[,"N.Models"],all.res[,"AICc"], xlab="Number of Piecewise Models", ylab="Model Fit", ylim=ylim, type="l", col="blue")
	points(all.res[,"N.Models"],all.res[,"AICc"], col="blue", pch=21, bg="white")
	points(all.res[,"N.Models"],all.res[,"AIC"], col="black", type="l")
	points(all.res[,"N.Models"],all.res[,"AIC"], col="black", pch=21, bg="white")
	
	legend("topleft", c("AICc","AIC"), pch=21, pt.bg="white", lty=1, col=c("blue", "black"), inset = .05, cex=0.75, bty="n") # 'bottomright' also works
}

# Works just as well as the one above
plot.AIC.alt <- function (all.res)
{
  ylim <- c(min(all.res[,"AIC"],all.res[,"AICc"]), max(all.res[,"AIC"],all.res[,"AICc"]))
	plot(all.res[,"N.Param"],all.res[,"AICc"], xlab="Number of Parameters", ylab="Model Fit", ylim=ylim, type="l", col="blue")
	points(all.res[,"N.Param"],all.res[,"AICc"], col="blue", pch=21, bg="white")
	par("new"=TRUE)
	plot(all.res[,"N.Param"],all.res[,"AIC"], xlab="", ylab="", xaxt="n", yaxt="n", ylim=ylim, type="l", col="black")
	points(all.res[,"N.Param"],all.res[,"AIC"], col="black", pch=21, bg="white")
	legend("topleft", c("AICc","AIC"), pch=21, pt.bg="white", lty=1, col=c("blue", "black"), inset = .05, cex=0.75, bty="n") # 'bottomright' also works
}

# Get b and d values from r (b-d) and epsilson (d/b)
get.b.d <- function(r, epsilon)
{
	b <- r/(1-epsilon)
	d <- b-r   # Alternatively: d <- eps*r/(1-eps)
	return(list(b=b, d=d))
}

logspace_sum <- function(logx)
{
	r <- logx[1]
	if (length(logx)>1)
	{
		for (i in 2:length(logx)) {r <- logspace_add(r, logx[i])}
	}
	return(r)	
}

# Calculate log(x+y) from log(x) and log(y)
logspace_add <- function(logx, logy)
{
	if (logx == -Inf)
	{
		return(logy)
	} else {
		max(logx, logy) + log1p(exp (-abs (logx - logy)))
	}
}



# Print out results from a given piecewise model. May be:
#   1. optimal model under AICc: the default
#   2. optimal model under AIC: through setting "aic=T"
#   3. AICc model using threshold improvement (like original MEDUSA): through setting "threshold" to some value
#   4. AIC model using threshold improvement (like original MEDUSA): through setting "aic=T" and "threshold" to some value
#   5. user-selected model: through setting "model" to some integer value (from summary table)
# In any case, also print out base model for comparison
# Need to parse 'results' (list of lists). Contains elements:
#  $models, which contains:
#     $par: i x 2 matrix (for birth-death; i x 1 for pure-birth); the jth row contains 
#       the speciation and (optional) extinction rate for the jth rate class
#     $lnLik.part: vector of length i; the jth element is the partial log
#       likelihood value due to the jth rate class
#     $lnLik: = sum(lnLik.part); the overall log-likelihood of this model
#     $split.at: the i+1 locations of splits.  The first element is the root node (i.e. background rate).
#  $phy: the tree, which may have been manipulated (renamed, pruned) in the analysis
#  $z: a matrix listing branch times and richnesses; for summarizing results
#  $anc: a list of all ancestors; for summarizing results
#  $model.summary: data.frame listing considered models, likelihoods, number of paramters, AICc, AICc, Akaike weights.

summarize.MEDUSA <- function(results, model=NULL, threshold=NULL, aic=F, plotTree=T, time=T, cex=0.5)
{
# Desirables:
#  1. table listing parameter values of selected model
#  2. list parameters of base model
#  3. tree printed with colour-coded edges, node labels to indicate splits
	
# Extract constituent components from results
	fit <- results$models
	phy <- results$phy
	z <- results$z
	anc <- results$anc
	model.summary <- results$model.summary
	
# First, determine which model is desired
	model.id <- 0
	criterion <- character
	if (aic) {criterion <- "AIC"} else {criterion <- "AICc"}
	if (!is.null(model))
	{
		model.id <- model
	} else {
    
    if(is.null(threshold)) threshold<-.threshold.medusa(phy)
    
		if (threshold == 0)
		{
			model.id <- which.min(unlist(lapply(fit, "[[", criterion)))
		} else {   # Find best model using threshold criterion
      model.id <- 1
			while (1)
			{
				if ((model.id + 1) > length(fit)) break;
				if ((unlist(fit[[model.id]][criterion]) - unlist(fit[[model.id+1]][criterion])) < threshold) break;
				model.id <- model.id + 1
			}
		}
	}
		
	break.pts <- fit[[model.id]]$split.at
	opt.model <- as.data.frame(cbind(Split.node=break.pts, fit[[model.id]]$par))
	base.model <- as.data.frame(fit[[1]]$par)
	
	cat("\nEstimated parameter values for model ", model.id, ":\n\n", sep="")
	print.data.frame(opt.model, digits=5)
	opt.weight <- 0
	base.weight <- 0
	if (criterion == "AICc")
	{
		opt.weight <- model.summary$w.AICc[model.id]
		base.weight <- model.summary$w.AICc[1]
	} else {
		opt.weight <- model.summary$w.AIC[model.id]
		base.weight <- model.summary$w.AIC[1]
	}
	cat("\n", criterion, " weight for this model: ", opt.weight, "\n\n", sep="")
	
	if (model.id != 1)
	{
		cat("\nFor comparison, estimated values for the base model are:\n\n")
		print.data.frame(base.model, digits=5, row.names=F)
		cat("\n", criterion, " weight for the base model: ", base.weight, "\n\n", sep="")
	}
	
# Get desired tree-model conformation
	for (i in 1:length(break.pts))
	{
		tmp <- medusa.split(break.pts[i], z, anc)
		z <- tmp$z
	}
	
# Plot tree with purdy colours and labelled nodes (to better map between tree and table)
	if (plotTree)
	{
		margin <- F
		mm <- match(phy$edge[,2], z[,2])
		if (time) {margin=T}
		plot(phy, edge.color=z[mm,10], no.margin=!margin, cex=0.75)
		if (time) axisPhylo(cex.axis=0.75)
		
		for (i in  1:length(break.pts))
		{
			nodelabels(i, node= break.pts[i], frame = "c", font = 1, cex=0.5)
		}
	}
  return(length(break.pts))
}

















# *** NOT used ***

# Experiments in logspace:
#	for (i in 1:length(a)) # Loop over fossils
#	{
#		tmp <- numeric(a[i]-1)
#		for (k in 1:(a[i]-1)) # For fossil i loop over from 1 to n[i]-1
#		{
#			ai <- k
#			ni <- n[i]
#			min.a.n <- pmin(ai, ni)
#	   		j <- seq_len(min.a.n)
#			
#	# *** Should this be put on a log scale to preserve precision? ***
#	  # Tests with reasonable data yield the same answer either way
#		#	foo <- choose(ai, j) * choose(ni - 1, j - 1) * A[i]^(ai-j) * ((1-A[i])*(1-B[i]))^j * B[i]^(ni - j)
#			foo <- lchoose(ai, j) + lchoose(ni - 1, j - 1) + (ai-j)*lA[i] + j*l1mA1mB[i] + (ni-j)*lB[i]
#			tmp[k] <- logspace_sum(foo)
#	# Alternate log form:
#	  # foo <- lchoose(ai, j) + lchoose(ni - 1, j - 1) + (ai-j)*lA[i] + j*l1mA1mB[i] + (ni-j)*lB[i]
#	    }
#	#  	tmp[i] <- sum(foo)/(1-A[i])  # old conditioning; not appropriate for a > 1
#			
#   # Condition on survival outside of the sum, as it factors out of all the subexpressions
#		foo <- (logspace_sum(tmp))
#		sum.prob.fossil.pend[i] <- foo - l1mA.exp.a[i]
#	}


## Print out tree with ape-style node-numbering
## Really only useful for troubleshooting
plotNN <- function (phy, time=T, margin=T, cex=0.5) 
{
	phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge)
	plot.phylo(phy, show.node.label=TRUE, no.margin=!margin, cex=cex) #, label.offset=0.5)
	if (time && !margin) cat("Cannot plot time axis without a margin.\n")
	else if (time && margin) axisPhylo(cex.axis=0.75)
}

# Function to prune tree using richness information; simpler version of 'prune.tree.merge.data' above
# 'richness' is assumed to have minimally two columns, "taxon" and "n.taxa"
# May also include 'exemplar' column; in that case, rename relevant tip.label before pruning
pruneTree <- function(phy, richness)
{
# Rename exemplar taxa with taxon name in richness file
	if (!is.null(richness$exemplar))
	{
# Change relevant tip.labels in phy; individual 'exemplar' may be NA, use original tip.label.
# Ordering in richness file should NOT be assumed to match order of tip.labels
		i.na <- is.na(richness$exemplar)
		phy$tip.label[match(richness$exemplar[!i.na], phy$tip.label)] <- as.character(richness$taxon[!i.na])
	}
	
	if (length(phy$tip.label) != length(richness[,1]))
	{
# Prune tree down to lineages with assigned richness
		temp <- richness[, "n.taxa"]
		names(temp) <- richness[, "taxon"]
		pruned <- treedata(phy, temp)		# geiger function calling ape (namecheck)
		prunedTree <- pruned$phy
# Check the tree
	#	plotNN(prunedTree)					# Node numbers (ape-style) plotted
		phy <- prunedTree
	}
	return(phy)
}

.threshold.medusa<-function (phy) 
{
    if ("multiPhylo" %in% class(phy)) {
        phy <- phy[[1]]
    }
    N <- Ntip(phy)*2
    if(N < 4) stop("Medusa only works with trees that have 4 or more tip clades")
    
    a <- -35.9410523803326
    b <- 6.7372587299747
    c <- -0.100615083407549
    Offset <- 27.5166786643334
    y <- a * (N - b)^c + Offset
    if (y < 0) {
        y <- 0
    }
    cat("Appropriate AICc threshold for tree of ", N, " tips is: ", 
        y, ".\n\n", sep = "")
    return(y)
}


source("fossilMEDUSA_working.R")
load("Whale.Rdata")
load("Mustelids.Rdata")


# Run data treating fossils as exact counts, minimum counts, or not at all

whales.f.exact <- runFossilMEDUSA(phy=whale.tree, richness=whale.richness, model.limit=8, fossil.minimum=F)
whales.f.min <- runFossilMEDUSA(phy=whale.tree, richness=whale.richness, model.limit=8, fossil.minimum=T)
whales.extant <- runFossilMEDUSA(phy=whale.tree, richness=whale.ext.richness, model.limit=8, est.extinction=T)


mustelids.f.exact <- runFossilMEDUSA(phy=mustelid.tree, richness=mustelid.richness, fossil.minimum=F)
mustelids.f.min <- runFossilMEDUSA(phy=mustelid.tree, richness=mustelid.richness, fossil.minimum=T)
mustelids.extant <- runFossilMEDUSA(phy=mustelid.tree, richness=mustelid.ext.richness)

###testing with fungal data
fungus.tree<-read.tree("fung.phy")
fungus.richness<-read.csv("fungus_richness.csv")
fungus.exact <- runFossilMEDUSA(phy=fungus.tree, richness=fungus.richness)
summarize.MEDUSA (fungus.exact)

#testing with green plant data
plants.tree<-read.tree("plants.tree")
plants.richness<-read.csv('GreenPlants_spprichness2.csv')
plants.exact<-runFossilMEDUSA(phy=plants.tree, richness=plants.richness)

summarize.MEDUSA (plants.exact)


###testing with alternative fungal richness coding
fungus.tree<-read.tree("fung.phy")
fungus.richness<-read.csv("fungus_richness2.csv")
fungus.exact <- runFossilMEDUSA(phy=fungus.tree, richness=fungus.richness)
summarize.MEDUSA (fungus.exact)


mustelids.f.min <- runFossilMEDUSA(phy=mustelid.tree, richness=mustelid.richness, fossil.minimum=T)
mustelids.extant <- runFossilMEDUSA(phy=mustelid.tree, richness=mustelid.ext.richness)
# Results from fossilMEDUSA is a list of lists. Annoying, yes; but saves on re-parsing data for summary later.
# Luke suggests making this into a class with specific print functions. At present, contains elements:
#  $models, which contains:
#     $par: i x 2 matrix (for birth-death; i x 1 for pure-birth); the jth row contains 
#       the speciation and (optional) extinction rate for the jth rate class
#     $lnLik.part: vector of length i; the jth element is the partial log likelihood value due to the jth rate class
#     $lnLik: = sum(lnLik.part); the overall log-likelihood of this model
#     $split.at: the i+1 locations of splits.  The first element is the root node (i.e. background rate).
#  $phy: the tree, which may have been manipulated (renamed, pruned) in the analysis
#  $z: a matrix listing branch times and richnesses; for summarizing results
#  $anc: a list of all ancestors; for summarizing results
#  $model.summary: data.frame listing considered models, likelihoods, number of paramters, AICc, AICc, Akaike weights.
#
# e.g. to see fit of all candidate models for a given analysis (say, whales with fossils treated as minimum counts), type:
#  whales.f.min$model.summary
#
# The summary function has the following interface:
# summarize.MEDUSA <- function(results, model=NULL, threshold=0, aic=F, plotTree=T, time=T, cex=0.5)
# Print out results from a given piecewise model. May be:
#   1. optimal model under AICc: the default
#   2. optimal model under AIC: through setting "aic=T"
#   3. AICc model using threshold improvement (like original MEDUSA): through setting "threshold" to some value
#   4. AIC model using threshold improvement (like original MEDUSA): through setting "aic=T" and "threshold" to some value
#   5. user-selected model: through setting "model" to some integer value (from summary table)
#
# For example, if you want to plot an arbitrary model for whales with fossils treated as minimum counts:
#  summarize.MEDUSA (whales.f.min, model=7)

# Sumamrize using default AICc model:

summarize.MEDUSA (whales.f.exact)
summarize.MEDUSA (whales.f.min)
summarize.MEDUSA (whales.extant)

summarize.MEDUSA (mustelids.f.exact)
summarize.MEDUSA (mustelids.f.min)
summarize.MEDUSA (mustelids.extant)


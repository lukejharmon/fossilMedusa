#source("fossilMEDUSA_working.R")
#load("Whale.Rdata")
#load("Mustelids.Rdata")
#make sure your richness file looks like this:
whale.richness
#there are five columns named taxon, exemplar, n.taxa, n.fossils, f.time
# taxon is the name of the exact tip in the tree
#exemplar is the name of the group that tip represents
#n.taxa is the number of tip species in that group
#n.fossils is the number of fossils at some time interval
#f.time is the time of that observation (0 is the present time)

#So...if you do not have any fossil observations for a group like caimins, just keep all of the caimins in the richness files (name each and give a richness of 1). If you DO have oberservations for a group, just pick one member of the clade to represent the richness, then give the fossil richness, n.fossils, at the time you report it, f.time. the tip richness, n.taxa, will be the richness for the present-day.

# Run data treating fossils as exact counts, minimum counts, or not at all

whales.f.exact <- runFossilMEDUSA(phy=whale.tree, richness=whale.richness, model.limit=8, fossil.minimum=F)

#minimum counts, 
whales.f.min <- runFossilMEDUSA(phy=whale.tree, richness=whale.richness, model.limit=8, fossil.minimum=T)
#or not at all (no fossil data used to calculate richness)

whales.extant <- runFossilMEDUSA(phy=whale.tree, richness=whale.ext.richness, model.limit=8, est.extinction=T)

#summarize results
summarize.MEDUSA(whales.f.min)

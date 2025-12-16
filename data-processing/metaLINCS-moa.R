library(METALINCS)

## METALINCS----------------------------------------------------
# First we compute the connectivity enrichment    
res <- computeConnectivityEnrichment(de, nprune=0)
names(res)

# Now compute the MoA enrichment
moa <- computeMoaEnrichment(res) 
names(moa)

# Plot the drugs connectivity using plotDrugConnectivity()
plotDrugConnectivity(res, contr=1)

# Plot the mechanism of action using plotMOA()
plotMOA(moa, contr=1, type="drugClass", ntop=10)
plotMOA(moa, contr=1, type="targetGene", ntop=10)

# Plot the drugs activity map using plotActivationMap()
plotActivationMap(res, nterms = 60, nfc=20, rot=FALSE)

final_integrated@misc$res <- res
final_integrated@misc$moa <- moa

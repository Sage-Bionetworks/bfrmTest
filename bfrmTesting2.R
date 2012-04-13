# bfrmTesting2.R

require(synapseClient)
require(bfrm)
require(affy)
require(parallel)
require(hgu133a.db)

## READ IN BETA CATENIN DATA
egfrEnt <- downloadEntity(138511)
egfrExpr <- ReadAffy(filenames=file.path(egfrEnt$cacheDir, egfrEnt$files))

## NORMALIZE AND SUMMARIZE USING RMA
egfrExprNorm <- rma(egfrExpr, normalize = T, background = F)

egfrMat <- exprs(egfrExprNorm)

egfrSamp <- sapply(strsplit(egfrEnt$files, "_", fixed = T), "[[", 3)
egfrVec <- as.numeric(grepl("wtegf", egfrEnt$files, ignore.case = T))

## RUN BFRM IN THE NONEVOLUTIONARY SPARSE ANOVA MODE
foo <- bfrm(egfrMat, design = egfrVec)
mPPib <- foo@results$mPostPib
topProbeLogical <- mPPib[ , 2] >= 0.99
topProbeInd <- grep("TRUE", topProbeLogical)
topProbeIDs <- rownames(egfrMat)[topProbeInd]
topEGFRSyms <- as.character(mget(topProbeIDs, hgu133aSYMBOL, 
                                 ifnotfound = NA))

## RUN BFRM IN EVOLUTIONARY MODE
# In Asynchronous mode using {parallel}
# asyncJob <- parallel(fooFactor <- evolve(egfrMat, 
#                     init = as.numeric(topProbeInd),
#                     maxVarIter = 30,
#                     maxVars = length(topProbeInd))
#                      )



# In conventional mode
fooFactor <- evolve(egfrMat, 
                    init = as.numeric(topProbeInd),
                    maxVarIter = 30,
                    maxFacs = 5,
                    maxVars = length(topProbeInd))



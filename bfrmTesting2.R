# bfrmTesting2.R

require(synapseClient)
require(bfrm)
require(affy)
require(parallel)
require(hgu133a.db)
require(corpcor)

## READ IN BETA CATENIN DATA
egfrEnt <- downloadEntity(138511)
egfrExpr <- ReadAffy(filenames = file.path(egfrEnt$cacheDir, egfrEnt$files))

## NORMALIZE AND SUMMARIZE USING RMA
egfrExprNorm <- rma(egfrExpr, normalize = T, background = F)

egfrMat <- exprs(egfrExprNorm)

egfrSamp <- sapply(strsplit(egfrEnt$files, "_", fixed = T), "[[", 3)
egfrVec <- as.numeric(grepl("wtegf", egfrEnt$files, ignore.case = T))

## BUILD THE H MATRIX
# Find the AFFX control probesets
affxInd <- grep('AFFX', rownames(egfrMat), ignore.case = F)
controlMat <- egfrMat[affxInd, ]
controlSVD <- fast.svd(controlMat)
hRightNine <- controlSVD$v[ , 1:9]
# hMatrix <- cbind(rep(1, 19), egfrVec, hRightFive)

## RUN BFRM IN THE NONEVOLUTIONARY SPARSE ANOVA MODE
foo <- bfrm(egfrMat, design = egfrVec, control = t(hRightNine))
mPPib <- foo@results$mPostPib
topProbeLogical <- mPPib[ , 2] >= 0.99
topProbeInd <- grep("TRUE", topProbeLogical)
topProbeIDs <- rownames(egfrMat)[topProbeInd]
topEGFRSyms <- as.character(mget(topProbeIDs, hgu133aSYMBOL, 
                                 ifnotfound = NA))

## BRING IN THE GGHEAT CODE ENTITY
ggheatEnt <- loadEntity(274063)
attach(ggheatEnt)

## VISUALIZING THE EGFR DATA SUBSETTED TO THE ANOVA SELECTED INDICES
anovaHeatPlot <- ggheat2(egfrMat[topProbeInd, ], clustering = 'both')


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
                    control = t(hRightNine),
                    priorpsia = 2,
                    priorpsib = 0.005,
                    varThreshold = 0.85,
                    facThreshold = 0.95,
                    maxVarIter = 30,
                    minFacVars = 10,
                    maxFacVars = length(topProbeInd),
                    maxFacs = 50,
                    maxVars = length(topProbeInd)
                    )

# OK, this seems to work. We therefore will need to change the default values 
# for priorpsia and priorpsib

# Now to evaluate the 'projection()' function:
projFoo <- projection(fooFactor, egfrMat)
par(mfrow = c(2,1))
image(fooFactor@results$mF)
image(t(projFoo)) # The projection result is transposed




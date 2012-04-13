# bfrmTesting.R

require(affy)
require(bfrm)
require(synapseClient)
require(Biobase)
require(hgu133a.db)

bcExpr <- ReadAffy(filenames = as.character(list.files("celFiles/", full.names = T)))

## NORMALIZE AND SUMMARIZE USING RMA
bcExprNorm <- rma(bcExpr, normalize = T, background = F)

bcMat <- exprs(bcExprNorm)

bcSamp <- sapply(strsplit(list.files("celFiles"), "_", fixed=T), "[[", 3)
bcBcat <- as.numeric(grepl("bcat", list.files("celFiles")))

## RUN BFRM IN THE SPARSE ANOVA MODE
foo <- bfrm(bcMat, design = bcBcat)
mPPib <- foo@results$mPostPib
topProbeInd <- mPPib[ , 2] >= 0.50
topProbeInd <- grep("TRUE", topProbeLogical)
topProbeIDs <- rownames(bcMat)[topProbeInd]
topBCatSyms <- as.character(mget(topProbeIDs, hgu133aSYMBOL, 
                                ifnotfound = NA))



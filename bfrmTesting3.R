# bfrmTesting3.R

require(synapseClient)
require(Biobase)
require(car)
require(corpcor)
require(ggplot2)
require(snm)
require(bfrm)

## LOAD IN THE EGFR RMA-NORMALIZED DATA
egfrEnt <- loadEntity(299165)

## PULL OUT THE RELEVANT DATA OBJECTS
egfrEset <- egfrEnt$objects$rmaEset
egfrMeta <- egfrEnt$objects$treatment
egfrExpress <- exprs(egfrEset)

## FOR DIAGNOSTICS, GENERATE A PRINCOMP PLOT
svdExp <- fast.svd(egfrExpress)
svdDF <- as.data.frame(svdExp$v[ , 1:4])
colnames(svdDF) <- paste("PrinComp", 1:4, sep = "")
rownames(svdDF) <- colnames(egfrExpress)

# A PRINCOMP PLOT UNADJUSTED
pcPlot1 <- ggplot(svdDF, aes(PrinComp1, PrinComp2)) +
  geom_point() +
  opts(title = "PC Plot of Ectopic EGFR Data") +
  ylab("PrinComp2") +
  xlab("PrinComp1") +
  opts(plot.title = theme_text(size = 14))

# USE SNM TO REMOVE THE TREATMENT EFFECT AND VISUALIZE FOR UNKNOWN ARTEFACT
adjustMM = model.matrix(~ factor(egfrMeta))
tmpFit <- snm(egfrExpress, adj.var = model.matrix(~ factor(egfrMeta)), 
              int.var = data.frame(factor(1:ncol(egfrExpress))),
              rm.adj = T)
svdTmp <- fast.svd(tmpFit$norm.dat - rowMeans(tmpFit$norm.dat), tol=0)
tmpSvdDF <- as.data.frame(svdTmp$v[ , 1:2])
colnames(tmpSvdDF) <- paste("PrinCompTxAdj", 1:2, sep = "")
rownames(tmpSvdDF) <- colnames(egfrExpress)

# A PRINCOMP PLOT OF ADJUSTED
pcPlot2 <- ggplot(tmpSvdDF, aes(PrinCompTxAdj1, PrinCompTxAdj2)) +
  geom_point(aes(colour = factor(egfrMeta))) +
  opts(title = "PC Plot of Ectopic EGFR Data after Removing Treatment") +
  ylab("PrinCompTxAdj2") +
  xlab("PrinCompTxAdj1") +
  opts(plot.title = theme_text(size = 14))
pcPlot2


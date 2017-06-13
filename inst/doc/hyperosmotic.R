## ----options, echo=FALSE----------------------------------------------------------------
options(width = 90)

## ----canprot, message=FALSE-------------------------------------------------------------
library(canprot)
data(canprot)

## ----datasets---------------------------------------------------------------------------
datasets <- pdat_osmotic()

## ----comptab, message=FALSE-------------------------------------------------------------
comptab <- lapply_canprot(datasets, function(dataset) {
  pdat <- get_pdat(dataset, "pdat_osmotic")
  ZC_nH2O(pdat, plot.it=FALSE)
}, varlist="pdat_osmotic")

## ----xsummary, results="asis"-----------------------------------------------------------
library(xtable)
xsummary(comptab)

## ----diffplot, fig.width=6, fig.height=6, fig.align="center"----------------------------
col <- rep("black", length(datasets))
col[grepl("=ASC", datasets)] <- "orange"
diffplot(comptab, col=col)


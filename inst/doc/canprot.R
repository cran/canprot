## ----setup, include=FALSE---------------------------------------------
library(canprot)
oldopt <- options(width = 72)

## ----vignettes, echo = FALSE------------------------------------------
functions <- grep("pdat_", ls("package:canprot"), value = TRUE)
vignettes <- gsub("pdat_", "", functions)
vignettes

## ----reset, include=FALSE-----------------------------------------------------
options(oldopt)


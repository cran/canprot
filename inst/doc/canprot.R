## ----setup, include=FALSE-----------------------------------------------------
library(canprot)

## ----vignettes, echo = FALSE--------------------------------------------------
functions <- grep("pdat_", ls("package:canprot"), value = TRUE)
vignettes <- gsub("pdat_", "", functions)
vignettes


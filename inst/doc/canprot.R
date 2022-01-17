## ----setup, include=FALSE---------------------------------------------
library(canprot)
oldopt <- options(width = 72)

## ----HTML, include=FALSE----------------------------------------------
ZC <- "<i>Z</i><sub>C</sub>"
nH2O <- "<i>n</i><sub>H<sub>2</sub>O</sub>"
H2O <- "H<sub>2</sub>O"
O2 <- "O<sub>2</sub>"

## ----AG---------------------------------------------------------------
AG <- data.frame(Ala = 1, Gly = 1)

## ----AG_metrics-------------------------------------------------------
ZCAA(AG)
H2OAA(AG)

## ----AG_reaction, message = FALSE-------------------------------------
CHNOSZ::basis("QEC")
CHNOSZ::subcrt("alanylglycine", 1)$reaction

## ----LYSC_CHICK-------------------------------------------------------
AA <- CHNOSZ::pinfo(CHNOSZ::pinfo("LYSC_CHICK"))
CHNOSZ::protein.length(AA)
AA

## ----LYSC_metrics-----------------------------------------------------
ZCAA(AA)
H2OAA(AA)

## ----LYSC_reaction, message = FALSE-----------------------------------
CHNOSZ::subcrt("LYSC_CHICK", 1)$reaction

## ----protcomp---------------------------------------------------------
(pc <- protcomp("P24298"))
ZCAA(pc$aa)

## ----protcomp_archaea-------------------------------------------------
aa_file <- system.file("extdata/aa/archaea/JSP+19_aa.csv.xz", package = "canprot")
pc <- protcomp("D4GP79", aa_file = aa_file)
ZCAA(pc$aa)

## ----GRAVY_pI---------------------------------------------------------
proteins <- c("LYSC_CHICK", "RNAS1_BOVIN", "AMYA_PYRFU")
AA <- CHNOSZ::pinfo(CHNOSZ::pinfo(proteins))
pI(AA)
GRAVY(AA)

## ----pdat_3D_datasets-------------------------------------------------
pdat_3D()

## ----pdat_3D----------------------------------------------------------
pdat <- pdat_3D("DKM+20")
str(pdat)

## ----get_comptab------------------------------------------------------
(get_comptab(pdat))

## ----get_comptabs, message = FALSE, results = "hide"------------------
datasets <- pdat_3D()
pdats <- lapply(datasets, pdat_3D)
comptabs <- lapply(pdats, get_comptab)

## ----diffplot, fig.width=5, fig.height=5, fig.align = "center"--------
diffplot(comptabs)
title("3D cell culture vs monolayers")

## ----vignettes, echo = FALSE------------------------------------------
functions <- grep("pdat_", ls("package:canprot"), value = TRUE)
vignettes <- gsub("pdat_", "", functions)
vignettes

## ----reset, include=FALSE-----------------------------------------------------
options(oldopt)


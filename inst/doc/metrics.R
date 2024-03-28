## ----setup, include=FALSE---------------------------------------------
library(canprot)
library(CHNOSZ)
oldopt <- options(width = 72)

## ----HTML, include=FALSE----------------------------------------------
Zc <- "<i>Z</i><sub>C</sub>"
nC <- "<i>n</i><sub>C</sub>"
nH2O <- "<i>n</i><sub>H<sub>2</sub>O</sub>"
H2O <- "H<sub>2</sub>O"
O2 <- "O<sub>2</sub>"

## ----AG_reaction, message = FALSE-------------------------------------
CHNOSZ::basis("QEC")
CHNOSZ::subcrt("alanylglycine", 1)$reaction

## ----AG_metrics-------------------------------------------------------
AG <- data.frame(Ala = 1, Gly = 1)
nH2O(AG, terminal_H2O = 1)

## ----LYSC_CHICK_metrics-----------------------------------------------
AA <- CHNOSZ::pinfo(CHNOSZ::pinfo("LYSC_CHICK"))
plength(AA)
Zc(AA)
nH2O(AA)

## ----LYSC_reaction, message = FALSE, echo = 1:2-----------------------
(reaction <- CHNOSZ::subcrt("LYSC_CHICK", 1)$reaction) # print the reaction
(H2Ocoeff <- with(reaction, coeff[name == "water"])) # print the coefficient on H2O
stopifnot( -(H2Ocoeff + 1) / plength(AA) == nH2O(AA))

## ----CHNOSZ_proteins--------------------------------------------------
iprotein <- CHNOSZ::pinfo(c("LYSC_CHICK", "RNAS1_BOVIN", "AMYA_PYRFU", "CSG_HALJP"))
AAcomp <- CHNOSZ::pinfo(iprotein)

## ----GRAVY------------------------------------------------------------
G_calc <- GRAVY(AAcomp)
# https://web.expasy.org/cgi-bin/protparam/protparam1?P00698@19-147@
# https://web.expasy.org/cgi-bin/protparam/protparam1?P61823@27-150@
# https://web.expasy.org/cgi-bin/protparam/protparam1?P49067@2-649@
G_ref <- c(-0.472, -0.663, -0.325)
stopifnot(all.equal(round(G_calc[1:3], 3), G_ref, check.attributes = FALSE))

## ----pI---------------------------------------------------------------
pI_calc <- pI(AAcomp)
# Reference values calculated with ProtParam
# LYSC_CHICK: residues 19-147 (sequence v1)
# RNAS1_BOVIN: residues 27-150 (sequence v1)
# AMYA_PYRFU: residues 2-649 (sequence v2)
# CSG_HALJP: residues 35-862 (sequence v1)
pI_ref <- c(9.32, 8.64, 5.46, 3.37)
stopifnot(all.equal(as.numeric(pI_calc), pI_ref))

## ----MW---------------------------------------------------------------
# Per-residue molecular weight multiplied by number of residues
MWcalc <- MW(AAcomp) * plength(AAcomp)
# Add terminal groups
MWcalc <- MWcalc + 18.01528
# Reference values for molecular weights of proteins
MWref <- c(14313.14, 13690.29, 76178.25)
stopifnot(all.equal(round(MWcalc[1:3], 2), MWref, check.attributes = FALSE))

## ----reset, include=FALSE-----------------------------------------------------
options(oldopt)


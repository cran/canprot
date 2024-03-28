## ----setup, include=FALSE-----------------------------------------------------
library(canprot)
library(CHNOSZ)
oldopt <- options(width = 80)

## ----HTML, include=FALSE------------------------------------------------------
Zc <- "<i>Z</i><sub>C</sub>"
nC <- "<i>n</i><sub>C</sub>"
nH2O <- "<i>n</i><sub>H<sub>2</sub>O</sub>"
H2O <- "H<sub>2</sub>O"
O2 <- "O<sub>2</sub>"

## ----read_fasta_KHAB17--------------------------------------------------------
fasta_file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
aa <- read_fasta(fasta_file)

## ----aa_KHAB17----------------------------------------------------------------
aa

## ----Zc_KHAB17, out.width = "75%", fig.align = "center", fig.width = 6, fig.height = 4----
xlab <- "Ancestral sequences (older to younger)"
plot(Zc(aa), type = "b", xaxt = "n", xlab = xlab, ylab = cplab$Zc)
names <- gsub(".*_", "", aa$protein)
axis(1, at = 1:6, names)
abline(v = 3.5, lty = 2, col = 2)
axis(3, at = 3.5, "GOE (proposed)")

## ----human_aa_ALAT1-----------------------------------------------------------
(aa <- human_aa("P24298"))
Zc(aa)

## ----DKM20_UniProt------------------------------------------------------------
up <-   c("Q92743", "P43490", "P52895", "P98160", "P23142", "P17301",
"U3KQK0", "Q15582", "Q9HCJ1", "P36222", "P27701", "Q08380", "P08572",
"P00734", "P22413", "O43657", "P35625", "O75348", "P02649", "P13861",
"P10620", "Q9H3N1", "A8K878", "P13611", "P07305", "E7ESP4", "Q9Y625",
"Q5ZPR3", "P62266", "Q96AQ6", "Q8N357", "Q13217", "Q9Y230", "Q9Y639",
"Q86W92", "C9JF17", "Q96PK6", "O95671", "P01033", "Q13501", "P69905",
"Q9Y5X1", "P50281", "Q9UBG0", "O60831", "P02751", "O43854", "P61803",
"J3KN66", "P42765", "P36543", "P15121", "Q16563", "Q12884", "P27695",
"P12110", "P07686", "Q92598", "Q02818", "Q07954", "O60493", "P40939",
"Q9Y3I0", "P51149", "P46776", "P46778", "P62805")

down <- c("J3KN67", "Q9Y490", "J3KNQ4", "E7EVA0", "Q01082", "J3KQ32",
"P54136", "Q9Y696", "Q01995", "Q15404", "P62714", "Q09666", "P07814",
"E7EQR4", "P46821", "O75369", "P02452", "P08123", "P54577", "P01023",
"Q6ZN40", "P42224", "B4DUT8", "Q13443", "Q9HCE1", "Q6DKJ4", "P50552",
"P35222", "P20908", "Q15417", "O75822", "P17812", "P05997", "P04080",
"O43294", "P08243", "P02458")

## ----DKM20_plot, out.width = "75%", fig.align = "center", fig.width = 7, fig.height = 5----
aa_down <- human_aa(down)
aa_up <- human_aa(up)
bp_names <- paste0(c("Down (", "Up ("), c(nrow(aa_down), nrow(aa_up)), c(")", ")"))

par(mfrow = c(1, 2))

Zclist <- list(Zc(aa_down), Zc(aa_up))
names(Zclist) <- bp_names
boxplot(Zclist, ylab = cplab$Zc, col = c(4, 2))
names(Zclist) <- c("x", "y")
p <- do.call(wilcox.test, Zclist)$p.value
legend("bottomleft", paste("p =", round(p, 3)), bty = "n")
title("Cabon oxidation state", font.main = 1)

nH2Olist <- list(nH2O(aa_down), nH2O(aa_up))
names(nH2Olist) <- bp_names
boxplot(nH2Olist, ylab = cplab$nH2O, col = c(4, 2))
names(nH2Olist) <- c("x", "y")
p <- do.call(wilcox.test, nH2Olist)$p.value
legend("bottomleft", paste("p =", round(p, 3)), bty = "n")
title("Stoichiometric hydration state", font.main = 1)

## ----reset, include=FALSE-----------------------------------------------------
options(oldopt)


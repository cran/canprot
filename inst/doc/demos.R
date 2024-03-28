## ----setup, include=FALSE-----------------------------------------------------
library(canprot)
library(CHNOSZ)
oldopt <- options(width = 80)

## ----HTML, include=FALSE------------------------------------------------------
Zc <- "<i>Z</i><sub>C</sub>"

## ----echo = FALSE-------------------------------------------------------------
knitr::read_chunk("../demo/thermophiles.R")

## ----thermophiles_demo_body, out.width = "100%", fig.align = "center", fig.width = 8, fig.height = 6, echo = FALSE, message = FALSE, dpi = 150----
# Change this to TRUE to make the image for the README
# (with only the Nitrososphaeria plot)
do.pdf <- FALSE
if(do.pdf) pdf("thermophiles.pdf", width = 4.2, height = 4.2)
if(!do.pdf) par(mfrow = c(1, 2))

par(mar = c(4, 4.5, 2, 1))

# Define metrics to calculate
metrics <- c("Zc", "S0g")
xlab <- cplab[[metrics[1]]]
ylab <- cplab[[metrics[2]]]

# Define common y-axis limit
ylim <- c(0.3000, 0.3045)

# Define colors
col2.8 <- adjustcolor(2, alpha.f = 0.8)
col4.8 <- adjustcolor(4, alpha.f = 0.8)
col2.5 <- adjustcolor(2, alpha.f = 0.5)

if(!do.pdf) {

## Plot 1: methanogen genomes

# Read amino acid composition
aafile <- system.file("extdata/aa/methanogen_aa.csv", package = "canprot")
aa <- read.csv(aafile)
# Set fill color for thermophiles
ithermo <- aa$ref > 50
bg <- ifelse(ithermo, col2.8, col4.8)
# Set point symbol for Class I and II
pch <- ifelse(aa$abbrv == "Class I", 21, 22)
# Calculate and plot metrics
metvals <- calc_metrics(aa, metrics)
plot(metvals, pch = pch, bg = bg, xlab = xlab, ylab = ylab, ylim = ylim)
# Add convex hull around thermophiles
add_hull(metvals[ithermo, ], col = col2.5, border = NA)
# Add legend and title
text(-0.225, 0.303, "Thermophiles", col = 2)
text(-0.18, 0.3003, "Mesophiles", col = 4)
legend("bottomleft", c("Class I", "Class II"), pch = c(21, 22), cex = 0.9)
title("Methanogen genomes", font.main = 1)

}

## Plot 2: Nitrososphaeria MAGs

# Read genome data
datfile <- system.file("extdata/aa/nitrososphaeria_MAGs.csv", package = "canprot")
dat <- read.csv(datfile)
# Read amino acid compositions
aafile <- system.file("extdata/aa/nitrososphaeria_aa.csv", package = "canprot")
aa <- read.csv(aafile)
# Match accessions
idat <- match(aa$organism, dat$Accession)
dat <- dat[idat, ]
# Filter by family
ifam <- !dat$Family %in% c("Nitrosopumilaceae", "Nitrososphaeraceae")
dat <- dat[ifam, ]
aa <- aa[ifam, ]
# Assign point symbol and color
pch <- sapply(dat$Respiration.type, switch, Anaerobic = 21, Aerobic = 22, Microaerobic = 23)
bg <- sapply(dat$Habitat.type, switch, Thermal = col2.8, Nonthermal = col4.8)
# Calculate and plot metrics
metvals <- calc_metrics(aa, metrics)
# Make plot
plot(metvals, pch = pch, bg = bg, xlab = xlab, ylab = ylab, ylim = ylim)
# Add convex hull around MAGs from thermal habitats
ithermal <- dat$Habitat.type == "Thermal"
add_hull(metvals[ithermal, ], col = col2.5, border = NA)
# Add legend and title
text(-0.215, 0.303, "Thermal habitats", col = 2)
text(-0.175, 0.3013, "Nonthermal habitats", col = 4)
legend("bottomright", c("Anaerobic", "Aerobic", "Microaerobic"), pch = c(21, 22, 23), cex = 0.9)
title("Nitrososphaeria MAGs", font.main = 1)

if(do.pdf) dev.off()

## ----echo = FALSE-------------------------------------------------------------
knitr::read_chunk("../demo/locations.R")

## ----locations_demo_setup, echo = FALSE---------------------------------------
# Change to TRUE to make a PDF image
do.pdf <- FALSE
if(do.pdf) pdf("locations.pdf", width = 8, height = 4.5)

## ----locations_demo_body, out.width = "100%", fig.align = "center", fig.width = 8, fig.height = 4.5, dpi = 150----
# Read SI table
file <- system.file("extdata/protein/TAW+17_Table_S6_Validated.csv", package = "canprot")
dat <- read.csv(file)
# Keep only proteins with validated location
dat <- dat[dat$Reliability == "Validated", ]
# Keep only proteins with one annotated location
dat <- dat[rowSums(dat[, 4:32]) == 1, ]

# Get the amino acid compositions
aa <- human_aa(dat$Uniprot)
# Put the location into the amino acid data frame
aa$location <- dat$IF.main.protein.location

# Use top locations (and their colors) from Fig. 2B of Thul et al., 2017
locations <- c("Cytosol","Mitochondria","Nucleoplasm","Nucleus","Vesicles","Plasma membrane")
col <- c("#194964", "#2e6786", "#8a2729", "#b2333d", "#e0ce1d", "#e4d71c")
# Keep the proteins in these locations
aa <- aa[aa$location %in% locations, ]
## Keep only proteins with length between 100 and 2000
#aa <- aa[plength(aa) >= 100 & plength(aa) <= 2000, ]

# Get amino acid composition for proteins in each location
# (Loop over groups by piping location names into lapply)
aalist <- lapply(locations, function(location) aa[aa$location == location, ] )

# Setup plot
par(mfrow = c(1, 2))
titles <- c(Zc = "Carbon oxidation state", pI = "Isoelectric point")
# Calculate Zc and pI
for(metric in c("Zc", "pI")) {
  datlist <- lapply(aalist, metric)
  bp <- boxplot(datlist, ylab = cplab[[metric]], col = col, show.names = FALSE)
  add_cld(datlist, bp)
  # Make rotated labels
  x <- (1:6) + 0.1
  y <- par()$usr[3] - 1.5 * strheight("A")
  text(x, y, locations, srt = 25, adj = 1, xpd = NA)
  axis(1, labels = FALSE)
  title(titles[metric], font.main = 1)
}

## ----echo = FALSE-------------------------------------------------------------
knitr::read_chunk("../demo/redoxins.R")

## ----redoxins_demo_body, out.width = "75%", fig.align = "center", fig.width = 6, fig.height = 5, echo = FALSE, message = FALSE, dpi = 100----
# Data file with protein IDs, sequence start/stop positions, and midpoint potentials
data_file <- system.file("extdata/fasta/redoxin.csv", package = "canprot")
dat <- read.csv(data_file)
# Drop PDI (a human protein) 20240304
dat <- dat[dat$protein != "PDI", ]

# Read header lines
fasta_file <- system.file("extdata/fasta/redoxin.fasta", package = "canprot")
headers <- read_fasta(fasta_file, type = "header")
# Locate the sequences in the FASTA file
iseqs <- sapply(dat$ID, grep, x = headers)
# Loop over proteins
aalist <- lapply(1:nrow(dat), function(i) {
  # Read the amino acid composition of this protein
  read_fasta(fasta_file, iseq = iseqs[i], start = dat$start[i], stop = dat$stop[i])
})
aa <- do.call(rbind, aalist)

# Make ferredoxin-thioredoxin reductase dimer (variable chain/catalytic chain)
iFTR <- grep("FTR", dat$protein)
aa[iFTR[1], 6:24] <- colSums(aa[iFTR, 6:24])
aa <- aa[-iFTR[2], ]
dat$protein[iFTR[1]] <- paste(dat$protein[iFTR], collapse = ":")
dat <- dat[-iFTR[2], ]

# Calculate Zc
Zc_values <- Zc(aa)

# Point symbols for E. coli and spinach
pch <- rep(19, length(Zc_values))
pch[dat$organism=="spinach"] <- 0

# Start plot
par(las = 1)
plot(dat$E0, Zc_values, pch = pch, xlim = c(-450, -100), ylim = c(-0.28, -0.04),
  xlab = expression(list(italic(E)*degree*"'", mV)),
  ylab = expression(italic(Z)[C]))
# Add dashed lines
lines(dat$E0[dat$organism == "ecoli"], Zc_values[dat$organism == "ecoli"], lty = 2)
lines(dat$E0[dat$organism == "spinach"], Zc_values[dat$organism == "spinach"], lty = 2)

# Add labels
pos <- rep(1, length(Zc_values))
pos[dat$organism == "ecoli"] <- 4
pos[dat$organism == "spinach"] <- 2
dx <- dy <- numeric(length(Zc_values))
dx[dat$protein == "DsbA"] <- -10
dy[dat$protein == "DsbA"] <- -0.012
dx[dat$protein == "DsbC"] <- -20
dy[dat$protein == "DsbC"] <- -0.012
text(dat$E0 + dx, Zc_values + dy, dat$protein, pos = pos)

# Add legend
legend("bottomright", pch = c(19, 0), legend = c(expression(italic("E. coli")), "spinach"))

## ----reset, include=FALSE-----------------------------------------------------
options(oldopt)


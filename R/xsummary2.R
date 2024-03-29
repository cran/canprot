# canH2O/R/xsummary2.R
# make table for new vignettes
# 20191206

xsummary2 <- function(comptab1, comptab2) {
  ct1 <- do.call(rbind, comptab1)
  ct2 <- do.call(rbind, comptab2)
  # get all data
  # include medians or means for up and down groups for making summary .csv files 20200125
  out <- data.frame(
    dataset = ct1$dataset,
    description = ct2$description,
    n1 = ct1$n1,
    n2 = ct1$n2,

    ZC.down = ct1$ZC.median1,
    ZC.up = ct1$ZC.median2,
    ZC.diff = ct1$ZC.diff,

    nH2O.down = ct1$nH2O.median1,
    nH2O.up = ct1$nH2O.median2,
    nH2O.diff = ct1$nH2O.diff,

    nAA.down = ct2$nAA.median1,
    nAA.up = ct2$nAA.median2,
    nAA.diff = ct2$nAA.diff,

    MW.down = ct2$MW.median1,
    MW.up = ct2$MW.median2,
    MW.diff = ct2$MW.diff,

    stringsAsFactors = FALSE
  )

  # prepare table
  x <- out[, c("dataset", "description", "n1", "n2", "ZC.diff", "nH2O.diff",
               "nAA.diff", "MW.diff")]
  # get the publication key from the dataset name
  publication <- sapply(strsplit(x$dataset, "_"), "[", 1)
  # format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # italicize species names (first words in description, surrounded by underscore)
  x$description <- gsub("^_", "<i>", x$description)
  x$description <- gsub("_", "</i>", x$description)
  # combine the publication and description
  x$description <- paste0(publication, " (", x$description, ")")
  # datasets have letter labels
  x$dataset <- c(letters, LETTERS)[1:nrow(x)]
  # to save column width, change "dataset" to "set"
  colnames(x)[1] <- "set"
  # multiply values of ZC and nH2O by 1000
  x[, 5:6] <- x[, 5:6] * 1000
  # multiply values of MW by 100
  x[, 8] <- x[, 8] * 100
  # round values
  x[, 5:8] <- round(x[, 5:8])

  # put markers around negative values
  for(icol in 5:8) {
    ineg <- x[, icol] < 0
    ineg[is.na(ineg)] <- FALSE
    x[ineg, icol] <- paste("**", x[ineg, icol], "**")
  }

  # create xtable
  x <- xtable::xtable(x, align=c("c", "l", "l", rep("r", ncol(x) - 2)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))

  # make the marked negative values bold
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # change "NaN" to "NA"
  x <- gsub("NaN", "NA", x, fixed=TRUE)

  # add headers that span multiple columns
  span_empty6 <- "<td align=\"center\" colspan=\"6\"></td>"
  span_empty2 <- "<td align=\"center\" colspan=\"2\"></td>"
  border <- paste("<table border=1> <tr>", span_empty6, span_empty2, "</tr>")
  x <- gsub("<table border=1>", border, x, fixed=TRUE)

  # more formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>down</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>up</sub>", x, fixed=TRUE)
  x <- gsub("ZC.diff", "&Delta;<i>Z</i><sub>C</sub>", x, fixed=TRUE)
  x <- gsub("nH2O.diff", "&Delta;<i>n</i><sub>H<sub>2</sub>O</sub>", x, fixed=TRUE)
  x <- gsub("nAA.diff", "&Delta;<i>n</i><sub>AA</sub>", x, fixed=TRUE)
  x <- gsub("MW.diff", "&Delta;MW", x, fixed=TRUE)

  # done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(out))
}

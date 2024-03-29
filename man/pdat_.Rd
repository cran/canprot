\encoding{UTF-8}
\name{pdat_}
\alias{pdat_}
\alias{pdat_breast}
\alias{pdat_colorectal}
\alias{pdat_lung}
\alias{pdat_liver}
\alias{pdat_pancreatic}
\alias{pdat_prostate}
\alias{pdat_hypoxia}
\alias{pdat_secreted}
\alias{pdat_3D}
\alias{pdat_glucose}
\alias{pdat_osmotic_bact}
\alias{pdat_osmotic_euk}
\alias{pdat_osmotic_halo}
\alias{.pdat_multi}
\alias{.pdat_osmotic}
\title{Get Protein Expression Data}
\description{
  Get differentially expressed proteins and amino acid compositions.
}

\usage{
  pdat_breast(dataset = 2020)
  pdat_colorectal(dataset = 2020)
  pdat_liver(dataset = 2020)
  pdat_lung(dataset = 2020)
  pdat_pancreatic(dataset = 2020)
  pdat_prostate(dataset = 2020)
  pdat_hypoxia(dataset = 2020)
  pdat_secreted(dataset = 2020)
  pdat_3D(dataset = 2020)
  pdat_glucose(dataset = 2020)
  pdat_osmotic_bact(dataset = 2020)
  pdat_osmotic_euk(dataset = 2020)
  pdat_osmotic_halo(dataset = 2020)
  .pdat_multi(dataset = 2020)
  .pdat_osmotic(dataset = 2017)
}

\arguments{
  \item{dataset}{character, dataset name}
}

\details{
The \code{pdat_} functions assemble lists of up- and down-regulated proteins and retrieve their amino acid compositions using \code{\link{protcomp}}.
The result can be used with \code{\link{get_comptab}} to make a table of chemical metrics that can then be plotted with \code{\link{diffplot}}.

If \code{dataset} is \samp{2020} (the default) or \samp{2017}, the function returns the names of all datasets in the compilation for the respective year.

Each dataset name starts with a reference key indicating the study (i.e. paper or other publication) where the data were reported.
The reference keys are made by combining the first characters of the authors' family names with the 2-digit year of publication.

If a study has more than one dataset, the reference key is followed by an underscore and an identifier for the particular dataset.
This identifier is saved in the variable named \code{stage} in the functions, but can be any descriptive text.

To retrieve the data, provide a single dataset name in the \code{dataset} argument.
Protein expression data is read from the CSV files stored in \code{extdata/expression/}, under the subdirectory corresponding to the name of the \code{pdat_} function.
Some of the functions also read amino acid compositions (e.g. for non-human proteins) from the files in \code{extdata/aa/}.

Descriptions for each function:

\itemize{
\item \code{pdat_colorectal}, \code{pdat_pancreatic}, \code{pdat_breast}, \code{pdat_lung}, \code{pdat_prostate}, and \code{pdat_liver} retrieve data for protein expression in different cancer types.

\item \code{pdat_hypoxia} gets data for cellular extracts in hypoxia and \code{pdat_secreted} gets data for secreted proteins (e.g. exosomes) in hypoxia.

\item \code{pdat_3D} retrieves data for 3D (e.g. tumor spheroids and aggregates) compared to 2D (monolayer) cell culture.

\item \code{.pdat_osmotic} retrieves data for hyperosmotic stress, for the 2017 compilation only.
In 2020, this compilation was expanded and split into \code{pdat_osmotic_bact} (bacteria), \code{pdat_osmotic_euk} (eukaryotic cells) and \code{pdat_osmotic_halo} (halophilic bacteria and archaea).

\item \code{pdat_glucose} gets data for high-glucose experiments in eukaryotic cells.

\item \code{.pdat_multi} retrieves data for studies that have multiple types of datasets (e.g. both cellular and secreted proteins in hypoxia), and is used internally by the specific functions (e.g. \code{pdat_hypoxia} and \code{pdat_secreted}).
}

}

\value{
A list consisting of:
\describe{
  \item{\code{dataset}}{Name of the dataset}
  \item{\code{description}}{Descriptive text for the dataset, used for making the tables in the vignettes (see \code{\link{mkvig}})}
  \item{\code{pcomp}}{UniProt IDs together with amino acid compositions obtained using \code{\link{protcomp}}}
  \item{\code{up2}}{Logical vector with length equal to the number of proteins; TRUE for up-regulated proteins and FALSE for down-regulated proteins}
}
}

\examples{
# List datasets in the 2017 complilation for colorectal cancer
pdat_colorectal(2017)
# Get proteins and amino acid compositions for one dataset
pdat_colorectal("JKMF10")
}

\concept{Protein data}

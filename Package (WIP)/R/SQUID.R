#' @title Quantify mitochondrial stress signatures in a dataset
#'
#' @description Generates a mitochondrial stress score for the inputted mRNA and/or metabolomics matrices using Stress Quantification Using Integrated Datasets (SQUID).
#'
#' @param mRNA Normalized gene expression counts from RNA sequencing. Count matrix should be arranged with genes as rows and samples as columns. At least one of mRNA or metabolomics must be inputted.
#' @param metabolites Normalized metabolite abundance levels from targeted metabolomics. Abundance matrix should be arranged with metabolites as rows and samples as columns. At least one of mRNA or metabolomics must be inputted.
#' @param metadata Sample metadata table. Sample ID column should be specified using the merge_by column.
#' @param geneType Gene ID form. Accepted values include Ensembl ID ('ENSEMBL'), gene symbol ('SYMBOL'), or Entrez ID ('ENTREZID).
#' @param species Species used for mRNA abundance. Must be 'human' or 'mouse'.
#' @param metLookup Lookup table to match inputted metabolites to metabolites from the mitochondrial stress dataset. If not provided, only the metabolites with matching names to the original mitochondrial stress dataset will be used.
#' @param ncomp Number of components (stress signatures) to generate using sPLS-DA
#' @param log.mRNA Whether mRNA counts have been log2 transformed
#' @param log.metabolites Whether metabolite abundances have been log2 transformed
#' @param seed Random number generator seed value for pseudo-pairing samples from existing mitochondrial stress dataset. Defaults to 1. Seed values can be compared using the compareSeeds function.
#' @param verbose Whether to output progress messages during stress score calculation
#' @param merge_by Column of metadata containing sample names for merging metadata and outputted stress scores
#'
#' @return An object containing some or all of the following:
#'
#' \item{output}{Maximum likelihood estimates for the parameters.}
#' \item{compPlot}{Plot showing component score patterns for components 1 through \var{ncomp}}
#' \item{data.mRNA}{mRNA gene expression counts used for analysis. Likely differs from the inputted mRNA table due to data processing prior to DIABLO analysis. Only returned if mRNA counts are inputted at start.}
#' \item{data.metabolites}{Metabolite abundance levels used for analysis. Likely differs from the inputted mRNA table due to data processing prior to DIABLO analysis. Only returned if metabolomics values are inputted at start.}
#' \item{scores.mRNA}{Stress scores associated with each sample from mRNA analysis. Only returned if mRNA counts are inputted at start.}
#' \item{scores.metabolites}{Stress scores associated with each sample from metabolite analysis. Only returned if metabolite data are inputted at start.}
#' \item{mergedScores.mRNA}{Stress scores associated with each sample from mRNA analysis, merged with inputted metadata. Only returned if mRNA counts and metadata table are inputted at start.}
#' \item{mergedScores.metabolites}{Stress scores associated with each sample from metabolite analysis. Only returned if metabolomics values and metadata table are inputted at start.}
#' \item{loadings.mRNA}{Gene loadings from DIABLO analysis, used to calculate stress score on inputted mRNA data. Only returned if mRNA counts are inputted at start.}
#' \item{loadings.metabolites}{Metabolite loadings from DIABLO analysis, used to calculate stress score on inputted metabolite data. Only returned if metabolite data are inputted at start.}
#'
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
#' @importFrom mixOmics "block.splsda"

SQUID <- function(mRNA = NULL,
                  metabolites = NULL,
                  metadata = NULL,
                  geneType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                  species = c("human", "mouse"),
                  metLookup = NULL,
                  ncomp = 4,
                  log.mRNA = F,
                  log.metabolites = F,
                  seed = 1,
                  verbose = T,
                  merge_by = 1) {
  print("Hello, world!")
}

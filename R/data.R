#' Contaminant protein information from Frankenfield et al. 2022
#'
#' A dataset containing contaminant protein information commonly found in mass
#' spectrometry experiments, compiled from Frankenfield et al. 2022.
#'
#' @format A data frame with 381 rows and 9 variables:
#' \describe{
#'   \item{Uniprot.ID}{UniProt protein identifier}
#'   \item{Entry.name}{UniProt entry name}
#'   \item{Status}{Review status (reviewed/unreviewed)}
#'   \item{Protein.names}{Full protein names}
#'   \item{Gene.names}{Gene names}
#'   \item{Organism}{Source organism}
#'   \item{Length}{Protein length in amino acids}
#'   \item{Source.of.Contamination}{Contamination source category}
#'   \item{Protein.Ids}{Copy of Uniprot.ID for compatibility}
#' }
#' @source Frankenfield et al. 2022 Supplemental Table S1
"Frankenfield_et_al_2022"
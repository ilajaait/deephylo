#' Create a phylogeny under the Constant Rate Birth-Death (CRBD) model
#'
#' @param n.tips number of tips
#' @param lambda speciation rate
#' @param mu extinction rate
#'
#' @return phylogeny (ape format)
createPhyloCRBD <- function(n.tips, lambda, mu) {
  lambda > mu || stop("Speciation rate inferior to extinction rate.")
  diversitree::trees(
    pars = c(lambda, mu), type = "bd", max.taxa = n.tips, n = 1
  )[[1]]
}

#' Compute total number of nodes of a phylogeny
#'
#' @param phylo phylogeny (ape format)
#'
#' @return int
numberNodesTotal <- function(phylo) {
  2 * length(phylo$tip.label) - 1
}

#' Compute number of tips of a phylogeny
#'
#' @param phylo phylogeny (ape format)
#'
#' @return int
numberTips <- function(phylo) {
  length(phylo$tip.label)
}

#' Compute number of internal nodes of a phylogeny
#'
#' @param phylo phylogeny (ape format)
#'
#' @return int
numberNodesInternal <- function(phylo) {
  length(phylo$tip.label) - 1
}

#' Compute phylogeny height
#'
#' @param phylo phylogeny (ape format)
#'
#' @return float
getHeight <- function(phylo) {
  max(castor::get_all_distances_to_root(phylo))
}

#' Is a phylogeny node a tip?
#'
#' @param phylo phylogeny (ape format)
#' @param i node index
#'
#' @return logical
isTip <- function(phylo, i){
  i <= numberTips(phylo)
}

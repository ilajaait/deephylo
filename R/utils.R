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

#' Compute number of edges of a phylogeny
#'
#' @param phylo (ape format)
#'
#' @return int
numberEdges <- function(phylo) {
  length(phylo$edge.length)
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

#' Find root of a phylogeny
#'
#' @param phylo phylogeny (ape format)
#'
#' @return int
getRoot <- function(phylo){
  numberTips(phylo) + 1
}

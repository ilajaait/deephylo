#' Create node dataframe
#'
#' Create a dataframe to fill with node attributes. Rows correspond to nodes,
#' column to attributes. All cells are filled with NAs.
#'
#' @param phylo phylogeny (ape format)
#'
#' @return data.frame
createNodesDf <- function(phylo) {
  n.nodes <- 2 * length(phylo$tip.label) - 1
  n.tips <- length(phylo$tip.label)
  df.nodes <- data.frame(
    index = rep(NA, n.nodes), # index of the node
    istip = rep(NA, n.nodes), # is the node a tip?
    distroot = rep(NA, n.nodes), # distance to the root
    part = rep(NA, n.nodes), # (3*dist) %/% h + 1 in {1;2;3}
    depth = rep(NA, n.nodes), # distance the root (in edges)
    colless = rep(NA, n.nodes), # colless score
    stair = rep(NA, n.nodes) # staircaseness
  )
  df.nodes
}

#### Node data frame ####

#' Create node dataframe
#'
#' Create a dataframe to fill with node attributes. Rows correspond to nodes,
#' column to attributes. All cells are filled with NAs.
#'
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
initializeNodesDf <- function(phylo) {
  n.nodes <- numberNodesTotal(phylo)
  n.tips <- numberTips(phylo)
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

#' Fill the "index" column of node dataframe
#'
#' Fill the "index" column of the nodes dataframe with their respective index.
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
#'
#' @seealso \code{\link{createNodesDf}}
fillNodesIndex <- function(df.nodes, phylo) {
  n.nodes <- numberNodesTotal(phylo)
  df.nodes["index"] <- 1:n.nodes
  df.nodes
}

#' Fill the "istip" column of node dataframe
#'
#' Fill the "istip" column of the node dataframe with corresponding logicals.
#' If the node is a tip `TRUE`, else `FALSE`.
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
#' @seealso \code{\link{createNodesDf}}
fillNodesIsTip <- function(df.nodes, phylo) {
  all(!is.na(df.nodes$index)) || stop("NA in column 'index'.")
  n.tips <- numberTips(phylo)
  df.nodes["istip"] <- isTip(phylo, df.nodes$index)
  df.nodes
}

#' Fill the "dist" column of node dataframe
#'
#' Fill the "dist" column of the node dataframe with corresponding distances
#' to root.
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
#' @seealso \code{\link{createNodesDf}},
#'     \code{\link[castor]{get_all_distances_to_root}}
fillNodesDist <- function(df.nodes, phylo) {
  df.nodes$distroot <- castor::get_all_distances_to_root(phylo)
  df.nodes
}

#' Fill the "part" column of node dataframe
#'
#' Fill the "part" column of the node dataframe with corresponding phylogeny
#' part. The phylogeny is cut in 3 equals regarding time.
#' Let's be D>0 the phylogeny height (maximal distance to the root),
#' then nodes with distance to the root d such that:
#' - d<D/3 belongs to part 1;
#' - D/3<d<2D/3 belongs to part 2;
#' - d>2D/3 belongs to part 3.
#' For more detail, see
#' \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}.
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)

#'
#' @return dataframe
#' @seealso \code{\link{createNodesDf}},
fillNodesPart <- function(df.nodes, phylo) {
  dist <- df.nodes$dist
  height <- getHeight(phylo)
  greater1 <- as.numeric(dist > height / 3)
  greater2 <- as.numeric(dist > as.numeric(2 * height / 3))
  df.nodes$part <- 1 + greater1 + greater2
  df.nodes
}

#' Fill the "colless" column of node dataframe
#'
#' Fill the "colless" column of the node dataframe with corresponding colless
#' score. The colless score is equal to the absolute difference of the number
#' of tips on the left and right side of the node.
#' For more detail, see
#' \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}.
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)

#'
#' @return dataframe
#' @seealso \code{\link{createNodesDf}},
#'     \code{\link{getNodeColless}}
fillNodesColless <- function(df.nodes, phylo) {
  n.nodes <- numberNodesTotal(phylo)
  istip <- df.nodes$istip
  for (i in 1:n.nodes) {
    !istip[i] && (df.nodes$colless[i] <- getNodeColless(phylo, i))
  }
  df.nodes
}

#' Fill the "stair" column of node dataframe
#'
#' Fill the "stair" column of the node dataframe with corresponding
#' stairecaseness. The stairecaseness is equal to the ratio of the number
#' of tips on each side.
#' For more detail, see
#' \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}.
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)

#'
#' @return dataframe
#' @seealso \code{\link{createNodesDf}},
#'     \code{\link{getNodeStair}}
fillNodesStair <- function(df.nodes, phylo) {
  n <- numberNodesTotal(phylo)
  istip <- df.nodes$istip
  for (i in 1:n) {
    !istip[i] && (df.nodes$stair[i] <- getNodeStair(phylo, i))
  }
  df.nodes
}

#' Fill the "depth" column of node dataframe
#'
#' Fill the "depth" column of the node dataframe with corresponding
#' depth. The depth is the number of edges between the node and the root.
#' For instance, for the following simple phylogeny root--->--->B,
#' depth(A) = 1 and depth(B) = 2
#' For more detail, see
#' \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}.
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)

#'
#' @return dataframe
#' @seealso \code{\link{createNodesDf}},
fillNodesDepth <- function(df.nodes, phylo) {
  df.nodes$depth <- castor::get_all_distances_to_root(phylo,
    as_edge_count = TRUE
  )
  df.nodes
}

#' Fill all columns of node dataframe
#'
#' @param df.nodes dataframe created by \code{\link{createNodesDf}}
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
fillNodesAll <- function(df.nodes, phylo) {
  df.nodes <- fillNodesIndex(df.nodes, phylo)
  df.nodes <- fillNodesIsTip(df.nodes, phylo)
  df.nodes <- fillNodesDist(df.nodes, phylo)
  df.nodes <- fillNodesPart(df.nodes, phylo)
  df.nodes <- fillNodesDepth(df.nodes, phylo)
  df.nodes <- fillNodesColless(df.nodes, phylo)
  df.nodes <- fillNodesStair(df.nodes, phylo)
  df.nodes
}

#' Create node dataframe
#'
#' Initialize and fill the node dataframe of the given phylogeny.
#'
#' @param phylo (ape format)
#'
#' @return dataframe
createNodesDf <- function(phylo) {
  phylo |>
    initializeNodesDf() |>
    fillNodesAll(phylo)
}

#' Compute the colless score of a node
#'
#' The colless score is defined as the difference of the number of left and
#' right tips of a node. Thus, this score can only be computed for internal
#' nodes (i.e. with at least one descendant).
#' For more detail, see
#' \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}.
#'
#' @param phylo phylogeny (ape format)
#' @param i node index
#'
#' @return int
getNodeColless <- function(phylo, i) {

  # Tests
  i > numberTips(phylo) || stop("Node is a tip, can't compute colless for tip.")
  i <= numberNodesTotal(phylo) || stop("Index out of range. Too large.")

  # Get node children
  children <- phangorn::Children(phylo, i)
  child1 <- children[1]
  child2 <- children[2]

  # Get their number of descendants
  n.descendant1 <- length(phangorn::Descendants(phylo, child1)[[1]])
  n.descendant2 <- length(phangorn::Descendants(phylo, child2)[[1]])

  abs(n.descendant1 - n.descendant2)
}

#' Compute the stairecaseness of a node
#'
#' Compute the ratio between the minimum and maximum number of tips on each side
#' for a node. Can be computed only for internal node.
#'
#' @param phylo phylogeny (ape format)
#' @param i node index
#'
#' @return float
getNodeStair <- function(phylo, i) {

  # Tests
  i > numberTips(phylo) || stop("Node is a tip, can't compute colless for tip.")
  i <= numberNodesTotal(phylo) || stop("Index out of range. Too large.")

  # Get node children
  children <- phangorn::Children(phylo, i)
  child1 <- children[1]
  child2 <- children[2]

  # Get their number of descendants
  n.descendant1 <- length(phangorn::Descendants(phylo, child1)[[1]])
  n.descendant2 <- length(phangorn::Descendants(phylo, child2)[[1]])
  n.desc <- c(n.descendant1, n.descendant2)

  min(n.desc) / max(n.desc)
}

#### end ####

#### Edge data frame ####

#' Create edge dataframe
#'
#' Create a dataframe to fill with edge attributes. Rows correspond to edges,
#' column to attributes. All cells are filled with NAs.
#'
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
initializeEdgesDf <- function(phylo) {
  n.edges <- numberEdges(phylo)
  df.edges <- data.frame(
    index = rep(NA, n.edges),
    isext = rep(NA, n.edges),
    length = rep(NA, n.edges),
    part = rep(NA, n.edges), # phylogeny part (1, 2 or 3)
    parent = rep(NA, n.edges),
    child = rep(NA, n.edges)
  )
  df.edges
}

#' Fill the "index" column of edge dataframe
#'
#' Fill the "index" column of the edge dataframe with their respective index.
#'
#' @param df.edges dataframe created by \code{\link{createEdgesDf}}
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
#'
#' @seealso \code{\link{createEdgesDf}}
fillEdgesIndex <- function(df.edges, phylo) {
  n.edges <- numberEdges(phylo)
  df.edges$index <- 1:n.edges
  df.edges
}

#' Fill the "isext" column of the edge dataframe
#'
#' An edge is said externed if it is linked to at least one tip.
#' The column is filled with logicals.
#' If the edge is external a `TRUE` is inserted, else it is a `FALSE`.
#'
#' @param df.edges dataframe
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
fillEdgesIsExt <- function(df.edges, phylo) {
  df.edges$isext <- isTip(phylo, df.edges$child)
  df.edges
}

#' Fill the "parent" and "child" columns of edge dataframe
#'
#' An edge link two species: a parent to its child.
#' As time goes, the edge goes from parent to child.
#'
#' @param df.edges dataframe
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
fillEdgesParentChild <- function(df.edges, phylo) {
  df.edges[c("parent", "child")] <- phylo$edge
  df.edges
}

#' Fill the "length" column of edge dataframe
#'
#' Insert edge lengths in the adequate cells of the edge dataframe.
#' Edge length corresponds to the time separating the parent from its child.
#'
#' @param df.edges dataframe
#' @param phylo phylogeny (ape format)
#'
#' @return dataframe
fillEdgesLength <- function(df.edges, phylo) {
  df.edges$length <- phylo$edge.length
  df.edges
}

#' Fill the "part" column of edge dataframe
#'
#' The phylogeny is sliced into 3 equal parts.
#' The part of the edge is determined by the part of the child,
#' e.g. if the edge goes from a parent belonging to part 3 and to a child
#' belonging to part 2, then the edge belongs to part 2.
#'
#' @param df.edges dataframe
#' @param df.nodes dataframe
#'
#' @return dataframe
#' @export
#' @examples
fillEdgesPart <- function(df.edges, df.nodes) {
  df.edges$part <- df.nodes$part[df.edges$child]
  df.edges
}

#' Fill all columns of edge dataframe
#'
#' @param df.edges dataframe
#' @param phylo phylogeny (ape format)
#' @param df.nodes dataframe
#'
#' @return dataframe
fillEdgesAll <- function(df.edges, phylo, df.nodes) {
  df.edges <- fillEdgesIndex(df.edges, phylo)
  df.edges <- fillEdgesLength(df.edges, phylo)
  df.edges <- fillEdgesParentChild(df.edges, phylo)
  df.edges <- fillEdgesPart(df.edges, df.nodes)
  df.edges <- fillEdgesIsExt(df.edges, phylo)
  df.edges
}

#' Create edge dataframe
#'
#' Initialize and fill the edge dataframe of the given phylogeny.
#'
#' @param phylo (ape format)
#' @param df.nodes (dataframe)
#'
#' @return dataframe
createEdgesDf <- function(phylo, df.nodes) {
  phylo |>
    initializeEdgesDf() |>
    fillEdgesAll(phylo, df.nodes)
}

#### end ####

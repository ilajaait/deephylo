#### Exported functions ####

#' Summary statistics of a phylogeny
#'
#' Compute 84 summary statistics on the phylogeny.
#' We split the summary statistics in 3 categories :
#'   * branch length, 25 summary statistics
#'   * topological, 8 summary statistics
#'   * LTT, 51 summary statistics
#' These summary statistics are listed below.
#'
#' Branch length summary statistics:
#'  \tabular{ll}{
#'   Variable \tab Description\cr
#'   height \tab phylogeny height\cr
#'   mean.all \tab mean of all branch lengths\cr
#'   median.all \tab median of all branch lengths\cr
#'   var.all \tab variance of all branch lengths\cr
#'   mean.ext \tab mean of external branch lengths\cr
#'   median.ext \tab median of external branch lengths\cr
#'   var.ext \tab variance of external branch lengths\cr
#'   mean.int1 \tab mean of internal branch lengths of part 1\cr
#'   median.int1 \tab median of internal branch lengths of part 1\cr
#'   var.int1 \tab variance of internal branch lengths of part 1\cr
#'   mean.int2 \tab mean of internal branch lengths of part 2\cr
#'   median.int2 \tab median of internal branch lengths of part 2\cr
#'   var.int2 \tab variance of internal branch lengths of part 2\cr
#'   mean.int3 \tab mean of internal branch lengths of part 3\cr
#'   median.int3 \tab median of internal branch lengths of part 3\cr
#'   var.int3 \tab variance of internal branch lengths of part 3\cr
#'   mean.intext1 \tab mean.int1 / mean.ext\cr
#'   median.int1 \tab median.int1 / median.ext\cr
#'   var.int1 \tab var.int1 / var.ext\cr
#'   mean.intext2 \tab mean.int2 / mean.ext\cr
#'   median.int2 \tab median.int2 / median.ext\cr
#'   var.int2 \tab var.int2 / var.ext\cr
#'   mean.intext3 \tab mean.int3 / mean.ext\cr
#'   median.int3 \tab median.int3 / median.ext\cr
#'   var.int3 \tab var.int3 / var.ext
#' }
#'
#' Topological summary statistics:
#'  \tabular{ll}{
#'   Variable \tab Description\cr
#'   colless \tab absolute difference in left and right tips, summed over int.
#'   nodes\cr
#'   sackin \tab depth summed over tips\cr
#'   widthdepth \tab ratio of maximal width over maximal depth\cr
#'   deltaw \tab maximal consecutive width difference\cr
#'   maxladder \tab maximal ladder size\cr
#'   inladder \tab proportion of nodes in ladder\cr
#'   imbalance \tab proportion of imbalanced nodes (i.e. colless !=0)\cr
#'   stair \tab ratio of tips on each side summed over internal nodes}
#'
#' LTT summary statistics:
#'  \tabular{ll}{
#'   Variable \tab Description\cr
#'   LTT_slope`i` \tab slope of LTT i-th part in semilog scale (i in \[1,10\])\cr
#'   LTT_t`i` \tab i-th time coordinate of binned LTT (i in \[1,20\])\cr
#'   LTT_N`i` \tab i-th lineage coordinate of binned LTT (i in \[1,20\])\cr
#'   n_tips \tab phylogeny number of tips
#'   }
#'
#' @param phylo phylogeny (ape format)
#'
#' @return list
#' @export
#'
#' @examples
#' phylo <- createPhyloCRBD(100, list(lambda = 1, mu = 0.1))
#' sumStats(phylo) # compute its summary statistics
sumStats <- function(phylo){

  # Set up dataframes
  df.nodes <- createNodesDf(phylo)
  df.edges <- createEdgesDf(phylo, df.nodes)

  # Compute summary statistics
  sumstat.bl <- sumStatsBranchLength(phylo, df.edges)
  sumstat.topo <- sumStatsTopo(phylo, df.nodes)
  sumstat.ltt <- sumStatsLtt(phylo)

  c(sumstat.bl, sumstat.topo, sumstat.ltt)
}

#### end ####

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

#### Summary statistics: branch length ####

#' Compute branch lengths summary statistics
#'
#' Compute the 25 summary statistics related to branch lengths listed below.
#' A branch is external if linked to a tip, other the branch is internal.
#' The branches are sorted in 3 parts regarding the time of the child.
#'  \tabular{ll}{
#'   Variable \tab Description\cr
#'   height \tab phylogeny height\cr
#'   mean.all \tab mean of all branch lengths\cr
#'   median.all \tab median of all branch lengths\cr
#'   var.all \tab variance of all branch lengths\cr
#'   mean.ext \tab mean of external branch lengths\cr
#'   median.ext \tab median of external branch lengths\cr
#'   var.ext \tab variance of external branch lengths\cr
#'   mean.int1 \tab mean of internal branch lengths of part 1\cr
#'   median.int1 \tab median of internal branch lengths of part 1\cr
#'   var.int1 \tab variance of internal branch lengths of part 1\cr
#'   mean.int2 \tab mean of internal branch lengths of part 2\cr
#'   median.int2 \tab median of internal branch lengths of part 2\cr
#'   var.int2 \tab variance of internal branch lengths of part 2\cr
#'   mean.int3 \tab mean of internal branch lengths of part 3\cr
#'   median.int3 \tab median of internal branch lengths of part 3\cr
#'   var.int3 \tab variance of internal branch lengths of part 3\cr
#'   mean.intext1 \tab mean.int1 / mean.ext\cr
#'   median.int1 \tab median.int1 / median.ext\cr
#'   var.int1 \tab var.int1 / var.ext\cr
#'   mean.intext2 \tab mean.int2 / mean.ext\cr
#'   median.int2 \tab median.int2 / median.ext\cr
#'   var.int2 \tab var.int2 / var.ext\cr
#'   mean.intext3 \tab mean.int3 / mean.ext\cr
#'   median.int3 \tab median.int3 / median.ext\cr
#'   var.int3 \tab var.int3 / var.ext
#' }
#'
#' @param df.edges dataframe
#' @param phylo phylogeny (ape format)
#'
#' @return list
#'    * `$stat.all`: statistic computed on all branch lengths
#'    (stat is either mean, median or var)
#'    * `$stat.ext`: statistic computed on external branch lengths only
#'    * `$stat.inti`: statistic computed on internal branch lengths only
#'    belonging to part i (i is either 1,2 or 3)
#'    * `$stat.intext1`: `$stat.inti` / `$stat.ext`
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}
sumStatsBranchLength <- function(phylo, df.edges) {

  # Set up
  sumstat <- list(
    mean.all = NA, median.all = NA, var.all = NA,
    mean.ext = NA, median.ext = NA, var.ext = NA,
    mean.int1 = NA, median.int1 = NA, var.int1 = NA,
    mean.int2 = NA, median.int2 = NA, var.int2 = NA,
    mean.int3 = NA, median.int3 = NA, var.int3 = NA,
    mean.intext1 = NA, median.intext1 = NA, var.intext1 = NA,
    mean.intext2 = NA, median.intext2 = NA, var.intext2 = NA,
    mean.intext3 = NA, median.intext3 = NA, var.intext3 = NA, height = NA
  )
  df.edges.ext <- df.edges[df.edges$isext, ] # internal edge
  df.edges.int <- df.edges[!df.edges$isext, ] # external edge
  branch_length.all <- df.edges$length
  branch_length.ext <- df.edges.ext$length
  branch_length.int1 <- df.edges.int[df.edges.int$part == 1, ]$length
  branch_length.int2 <- df.edges.int[df.edges.int$part == 2, ]$length
  branch_length.int3 <- df.edges.int[df.edges.int$part == 3, ]$length
  branch_length.list <- list(
    branch_length.all, branch_length.ext,
    branch_length.int1, branch_length.int2, branch_length.int3
  )
  n <- length(branch_length.list)

  # Fill sumstat for 'int', 'all' and 'int'
  for (i in 1:n) {
    branch_length <- branch_length.list[[i]]
    stats <- getBranchLengthStats(branch_length)
    sumstat[[3 * (i - 1) + 1]] <- stats$mean
    sumstat[[3 * (i - 1) + 2]] <- stats$median
    sumstat[[3 * (i - 1) + 3]] <- stats$var
  }

  # Fill sumstat for 'intext'
  for (i in 1:3) {
    for (stat_name in c("mean", "median", "var")) {
      name.int <- paste(stat_name, ".int", i, sep = "")
      name.intext <- paste(stat_name, ".intext", i, sep = "")
      name.ext <- paste(stat_name, ".ext", sep = "")
      sumstat[[name.intext]] <- sumstat[[name.int]] / sumstat[[name.ext]]
    }
  }

  sumstat$height <- getHeight(phylo)

  sumstat
}

#' Compute branch lengths statistics
#'
#' Compute the mean, median and variance of a list of branch lengths.
#'
#' @param branch_length vector
#'
#' @return list
#'    * `$mean`: branch lengths mean
#'    * `$median`: branch lengths median
#'    * `$var`: branch length variance
getBranchLengthStats <- function(branch_length) {
  length(branch_length) > 0 || return(list(mean = NA, median = NA, var = NA))
  list(
    mean = mean(branch_length),
    median = stats::median(branch_length),
    var = stats::var(branch_length)
  )
}

##### end ####

#### Summary statistics: topology ####

#' Topological summary statistics
#'
#' There are 8 topological summary statistics listed below.
#'  \tabular{ll}{
#'   Variable \tab Description\cr
#'   colless \tab absolute difference in left and right tips, summed over int.
#'   nodes\cr
#'   sackin \tab depth summed over tips\cr
#'   widthdepth \tab ratio of maximal width over maximal depth\cr
#'   deltaw \tab maximal consecutive width difference\cr
#'   maxladder \tab maximal ladder size\cr
#'   inladder \tab proportion of nodes in ladder\cr
#'   imbalance \tab proportion of imbalanced nodes (i.e. colless !=0)\cr
#'   stair \tab ratio of tips on each side summed over internal nodes}
#'
#' @param phylo phylogeny (ape format)
#' @param df.nodes dataframe
#'
#' @return list
#'     the variable listed above can be accessed with `$` + variable name.
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}
sumStatsTopo <- function(phylo, df.nodes) {

  # Set up
  sumstat <- list(
    "colless" = NA, "sackin" = NA, "widthdepth" = NA, "deltaw" = NA,
    "maxladder" = NA, "inladder" = NA, "imbalance" = NA, "stair" = NA
  )

  # Fill
  sumstat$colless <- sum(df.nodes$colless, na.rm = TRUE)
  sumstat$sackin <- getSackin(df.nodes)
  sumstat$widthdepth <- getWidthDepthRatio(df.nodes)
  sumstat$deltaw <- getMaxWidthDiff(df.nodes)
  sumstat$maxladder <- getMaxLadder(phylo)
  sumstat$inladder <- getInLadderNodesProp(phylo)
  sumstat$imbalance <- getImbalancedNodesProp(df.nodes)
  sumstat$stair <- sum(df.nodes$stair, na.rm = TRUE)

  sumstat
}

#' Is the node in a ladder?
#'
#' A node is in a ladder if one of its children is a tip AND the other is not.
#'
#' @param phylo phylogeny (ape format)
#' @param i node index
#'
#' @return logical
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}
isInLadder <- function(phylo, i) {
  !isTip(phylo, i) || return(FALSE)
  children <- phangorn::Children(phylo, i)
  sum(isTip(phylo, children)) == 1
}

#' Traverse a phylogeny and write ladders
#'
#' Traverse the phylogeny in preoder and write node indexes belonging to
#' ladders. Nodes that don't belong to ladder are marked with `-1`.
#' A node is in a ladder if one of its children is a tip AND the other is not.
#'
#' @param phylo phylogeny (ape format)
#' @param i node index
#' @param preorder vector
#'
#' @return vector with ladders
#'
#' @seealso \code{\link{isInLadder}}, \code{\link{getChildLeftRight}}
writeLadders <- function(phylo, i, preorder) {
  value <- -1
  !isInLadder(phylo, i) || (value <- i)
  preorder <- c(preorder, value)
  children <- getChildLeftRight(phylo, i)
  child.left <- children$right
  child.right <- children$left
  if (!is.na(child.left)) {
    preorder <- writeLadders(phylo, child.left, preorder)
  }
  if (!is.na(child.right)) {
    preorder <- writeLadders(phylo, child.right, preorder)
  }
  preorder
}

#' Left and right child of a node
#'
#' The left children is the node the further from the root, the right children
#' is the closest to the root.
#' If node is a tip returns NA.
#'
#' @param phylo phylogeny (ape format)
#' @param i node index
#'
#' @return list
#'    * `$left`: node left child
#'    * `$right`: node right child
getChildLeftRight <- function(phylo, i) {

  # Set up
  child <- list("left" = NA, "right" = NA)
  !isTip(phylo, i) || return(child) # tip has no child, return NA
  dist <- castor::get_all_distances_to_root(phylo)
  child.disordered <- phangorn::Children(phylo, i)

  if (all(isTip(phylo, child.disordered))) { # 2 tips, same distance, random
    child$left <- child.disordered[1]
    child$right <- child.disordered[2]
  } else {
    dist.child <- c(dist[child.disordered[1]], dist[child.disordered[2]])
    child$left <- child.disordered[which(dist.child == max(dist.child))]
    child$right <- child.disordered[which(dist.child == min(dist.child))]
  }

  child
}

#' Extract ladders from vector
#'
#' A node is in a ladder if one of its children is a tip AND the other is not.
#' Ladder of size one are deleted.
#'
#' @param ladder vector output of \code{\link{writeLadders}}
#'
#' @return list
extractLadders <- function(ladder) {
  n <- length(ladder)
  ladder.list <- list()
  l <- c()
  for (i in 1:n) {
    x <- ladder[i]
    if (x != -1) {
      l <- c(l, x)
    } else if (x == -1 & length(l) == 1) {
      l <- c()
    } else if (x == -1 & length(l) >= 2) {
      ladder.list <- append(ladder.list, list(l))
      l <- c()
    }
  }
  return(ladder.list)
}

#' Get ladders of a phylogeny
#'
#' A node is in a ladder if one of its children is a tip AND the other is not.
#' Ladder should be at least of size two.
#'
#' @param phylo phylogeny (ape format)
#'
#' @return list of ladders
#'
#' @seealso \code{\link{writeLadders}}, \code{\link{extractLadders}}
getLadders <- function(phylo) {
  root <- getRoot(phylo)
  ladders <- writeLadders(phylo, root, c())
  extractLadders(ladders)
}

#' Size of the largest ladder of the phylogeny normalized by the number of tips
#'
#' A node is in a ladder if one of its children is a tip AND the other is not.
#' Ladder should be at least of size two.
#'
#' @param phylo phylogeny (ape format)
#'
#' @return numeric
#' @seealso \code{\link{writeLadders}}, \code{\link{extractLadders}},
#' \code{\link{getLadders}}
getMaxLadder <- function(phylo) {
  ladders <- getLadders(phylo)
  length(ladders) > 0 || return(NA)
  max.ladder <- max(sapply(ladders, length))
  max.ladder / numberTips(phylo)
}

#' Proportion of internal nodes in ladder
#'
#' A node is in a ladder if one of its children is a tip AND the other is not.
#' Ladder should be at least of size two.
#'
#' @param phylo phylogeny (ape format)
#'
#' @return numeric
#'
#' @seealso \code{\link{writeLadders}}, \code{\link{extractLadders}},
#' \code{\link{getLadders}}
getInLadderNodesProp <- function(phylo) {
  ladders <- getLadders(phylo)
  nodes_inladder <- length(unlist(ladders))
  nodes_inladder / numberNodesInternal(phylo)
}

#' Proportion of imbalanced nodes
#'
#' A node is imbalanced if it has a different number of left and right tips.
#' Equivalently, the node has a stairecaseness different from 1.
#'
#' @param df.nodes dataframe
#'
#' @return numeric
getImbalancedNodesProp <- function(df.nodes) {
  n.imbalanced <- sum(df.nodes$stair != 1 & !is.na(df.nodes$stair))
  n.internal <- sum(!df.nodes$istip)
  n.imbalanced / n.internal
}

#' Width over depth of a phylogeny
#'
#' Depth of a node is the number of edges between the node and the root.
#' Width at a given depth is the number of nodes at the considered depth.
#' The width of the phylogeny is the maximum of the widths for all possible
#' depths.
#'
#' @param df.nodes dataframe
#'
#' @return numeric
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}
getWidthDepthRatio <- function(df.nodes) {
  width <- max(table(df.nodes$depth))
  depth <- max(df.nodes$depth)
  width / depth
}

#' Sackin score
#'
#' Sum of the depth over the phylogeny tips.
#'
#' @param df.nodes data.frame
#'
#' @return int
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}
getSackin <- function(df.nodes) {
  sum(df.nodes[df.nodes$istip, ]$depth)
}

#' Maximal consecutive width difference of a phylogeny
#'
#' Width at a given depth is the number of nodes at the considered depth.
#' The width of the phylogeny is the maximum of the widths for all possible
#' depths.
#'
#' @param df.nodes data.frame
#'
#' @return numeric
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}
getMaxWidthDiff <- function(df.nodes) {
  table.width <- table(df.nodes$depth)
  n <- length(table.width)
  max(abs(table.width[1:n - 1] - table(df.nodes$depth)[2:n]))
}

#### end ####

#### Summary statistics: LTT ####

#' LTT summary statistics
#'
#' There are 51 LTT summary statistics listed below.
#'  \tabular{ll}{
#'   Variable \tab Description\cr
#'   LTT_slope`i` \tab slope of LTT i-th part in semilog scale (i in \[1,10\])\cr
#'   LTT_t`i` \tab i-th time coordinate of binned LTT (i in \[1,20\])\cr
#'   LTT_N`i` \tab i-th lineage coordinate of binned LTT (i in \[1,20\])\cr
#'   n_tips \tab phylogeny number of tips
#'   }
#'
#' @param phylo phylogeny (ape format)
#'
#' @return list
#'     the variable listed above can be accessed with `$` + variable name.
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1005416}{Saulnier et al., 2014, PLOS}
sumStatsLtt <- function(phylo){
  ltt.coords <- getLttCoords(phylo)
  ltt.slopes <- getLttSlopes(phylo)
  c(ltt.coords, ltt.slopes, list(n_tips=numberTips(phylo)))
}

#' Sample 20 LTT coordinates
#'
#' Take 20 points uniformly distributed regarding the number of lineages.
#' Save LTT coordinates for these 20 points.
#' Phylogeny needs to be at least to have 20 tips.
#'
#' @param phylo phylogeny (ape format)
#'
#' @return list
#'    * `$LTT_coordX` with `coord` either `t` for time or `N` for number of
#'    lineages and `X` between 1 and 20
#'    * for instance, `$LTT_N15` is the number of lineages for bin 15 out of 20
getLttCoords <- function(phylo){

  # Test
  numberTips(phylo) >= 20 || stop("Phylogeny needs to have at least 20 tips.")

  # Set up
  ltt.coord <- ape::ltt.plot.coords(phylo)
  n_events <- length(ltt.coord[,1])
  bins <- 1:20 * (n_events %/% 20)

  # Get binned coordinates
  ltt.coord.binned <- ltt.coord[bins,]

  # Format output
  ltt.coord.binned <- as.list(ltt.coord.binned)
  names(ltt.coord.binned) <- lttCoordsNames()
  ltt.coord.binned
}

#' Create names for binned LTT coordinate list
#'
#' @return vector
lttCoordsNames <- function() {
  name.vec <- c()
  for (coord in c("t", "N")){
    for (i in 1:20){
      name <- paste("LTT_", coord, i, sep="")
      name.vec <- c(name.vec, name)
    }
  }
  name.vec
}

#' Compute 10 LTT slopes
#'
#' Divide the tree in 3 equals part (trough time) and compute the slope for each
#' part (named slope$p.i). Then do the ratio of:
#' - slope of the 2nd over the slope of the 1st part
#' - slope of the 3rd over the slope of the 2nd part
#'
#' @param phylo phylogeny (ape format)
#'
#' @return list
getLttSlopes <- function(phylo){
  numberTips(phylo) >= 20 || stop("Phylogeny needs to have at least 20 tips.")
  n_split <- 10
  slopes <- list()
  coords <- as.data.frame(ape::ltt.plot.coords(phylo))
  coords["N"] <- log(coords["N"])
  coords.split <- splitLttCoords(coords, n = n_split)
  for (i in 1:n_split){
      split <- coords.split[[i]]
      name <- paste("LTT_slope", i, sep="")
      s <- stats::lm(N ~ time, split)$coef[[2]] # compute slope
      slopes[[name]] <- s # save
  }
  slopes
}

#' Split the LTT coordinates regarding the number of events
#'
#' @param coords ltt coordinates (dataframe)
#' @param n number of splits (int, default n=10)
#'
#' @return list
splitLttCoords <- function(coords, n=10){
  coords.split <- list()
  len <- nrow(coords)
  bounds <- as.integer(seq(1, len, length.out=n+1))
  for (i in 1:n){
    inf <- bounds[i] + (i!=1)
    sup <- bounds[i+1]
    split <- coords[inf:sup,]
    coords.split[[i]] <- split
  }
  coords.split
}

#### end ####

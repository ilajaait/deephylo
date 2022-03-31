test_that("Node dataframe is initialized correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0.5)
  df.nodes <- initializeNodesDf(phylo)
  expect_equal(nrow(df.nodes), 2 * 10 - 1) # expected number of rows
  expect_equal(ncol(df.nodes), 7) # expected number of columns
  expect_true(all(is.na(df.nodes))) # filled with NAs
})

test_that("Node dataframe is created correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0.5)
  df.nodes <- initializeNodesDf(phylo)
  df.nodes <- fillNodesAll(df.nodes, phylo)
  expect_equal(df.nodes, createNodesDf(phylo))
})

test_that("Statistics for node dataframe are computed correctly.", {

  # Colless score
  phylo <- createPhyloCRBD(3, 1, 0.5)
  expect_equal(getNodeColless(phylo, 4), 1) # 2 vs. 1 tip
  expect_equal(getNodeColless(phylo, 5), 0) # 1 vs. 1 tip
  expect_error(getNodeColless(phylo, 3)) # node can't be a tip
  expect_error(getNodeColless(phylo, 6)) # node index out of bounds

  # Stairecaseness
  phylo <- createPhyloCRBD(3, 1, 0)
  expect_equal(getNodeStair(phylo, 4), 0.5)
  expect_equal(getNodeStair(phylo, 5), 1) # 1 vs. 1 tip
  expect_error(getNodeStair(phylo, 3)) # node can't be a tip
  expect_error(getNodeStair(phylo, 6)) # node index out of bounds
})

test_that("Fill node dataframe correctly.", {

  # Set up
  phylo <- createPhyloCRBD(10, 1, 0.5)
  df.nodes <- initializeNodesDf(phylo)

  # Check index
  df.nodes <- fillNodesIndex(df.nodes, phylo)
  expect_equal(df.nodes$index, 1:19)

  # Check 'is tip?'
  df.nodes$index[1] <- NA # introduce an error
  expect_error(fillNodesIsTip(df.nodes, phylo))
  df.nodes$index[1] <- 1 # restore dataframe
  df.nodes <- fillNodesIsTip(df.nodes, phylo)
  expect_equal(df.nodes$istip, as.logical(df.nodes$index <= numberTips(phylo)))

  # Check distance to root
  df.nodes <- fillNodesDist(df.nodes, phylo)
  root <- numberTips(phylo) + 1
  expect_equal(df.nodes$distroot[root], 0) # distance between root to root is null
  expect_equal(df.nodes$distroot, castor::get_all_distances_to_root(phylo))

  # Check node parts
  df.nodes <- fillNodesPart(df.nodes, phylo)
  height <- getHeight(phylo)
  dist <- df.nodes$dist
  part <- 1 + as.numeric(dist > height / 3) + as.numeric(dist > 2 * height / 3)
  expect_equal(df.nodes$part, part)

  # Check colless score
  df.nodes <- fillNodesColless(df.nodes, phylo)
  expect_true(all(is.na(df.nodes$colless[1:10]))) # NA for tips
  expect_true(all(!is.na(df.nodes$colless[11:19]))) # no NA for internal nodes

  # Check stairecaseness score
  df.nodes <- fillNodesStair(df.nodes, phylo)
  expect_true(all(is.na(df.nodes$stair[1:10]))) # NA for tips
  expect_true(all(!is.na(df.nodes$stair[11:19]))) # no NA for internal node

  # Check depth
  df.nodes <- fillNodesDepth(df.nodes, phylo)
  expect_equal(df.nodes$depth[root], 0)
  expect_true(all(df.nodes$depth >= 0))

  # Check fill all
  df.nodes.new <- initializeNodesDf(phylo)
  df.nodes.new <- fillNodesAll(df.nodes.new, phylo)
  expect_true(all(df.nodes == df.nodes.new, na.rm = TRUE))
})

test_that("Edge dataframe is initialized correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0.5)
  df.edges <- initializeEdgesDf(phylo)
  expect_equal(nrow(df.edges), numberEdges(phylo)) # expected number of rows
  expect_equal(ncol(df.edges), 6) # expected number of columns
  expect_true(all(is.na(df.edges))) # filled with NAs
})

test_that("Edge dataframe is created correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0.5)
  df.nodes <- createNodesDf(phylo)
  df.edges <- initializeEdgesDf(phylo)
  df.edges <- fillEdgesAll(df.edges, phylo, df.nodes)
  expect_equal(df.edges, createEdgesDf(phylo, df.nodes))
})

test_that("Fill edge dataframe correctly.", {

  # Set up
  phylo <- createPhyloCRBD(10, 1, 0.5)
  df.edges <- initializeEdgesDf(phylo)

  # Check index
  df.edges <- fillEdgesIndex(df.edges, phylo)
  expect_equal(df.edges$index, 1:numberEdges(phylo))

  # Check parent and child
  df.edges <- fillEdgesParentChild(df.edges, phylo)
  expect_true(all(df.edges[c("parent", "child")] == phylo$edge))

  # Check is edge external?
  df.edges <- fillEdgesIsExt(df.edges, phylo)
  expect_equal(sum(df.edges$isext), numberTips(phylo))
  expect_true(all(!is.na(df.edges$is.ext)))

  # Check edge lengths
  df.edges <- fillEdgesLength(df.edges, phylo)
  expect_equal(df.edges$length, phylo$edge.length)

  # Check edge part
  df.nodes <- createNodesDf(phylo)
  df.edges <- fillEdgesPart(df.edges, df.nodes)
  expect_equal(df.edges$part, df.nodes$part[df.edges$child])

  # Check fill all
  df.edges.new <- initializeEdgesDf(phylo)
  df.edges.new <- fillEdgesAll(df.edges, phylo, df.nodes)
  expect_equal(df.edges, df.edges.new)
})

test_that("Branch lengths summary statistis are computed correctly.", {

  # Check mean, median and variance
  expect_equal(
    getBranchLengthStats(c(1, 1, 1)),
    list(mean = 1, median = 1, var = 0)
  )
  expect_equal(
    getBranchLengthStats(c(1, 2, 3)),
    list(mean = 2, median = 2, var = 1)
  )
  expect_equal(
    getBranchLengthStats(c(1, 2, 6)),
    list(mean = 3, median = 2, var = 7)
  )
  expect_equal(
    getBranchLengthStats(c()),
    list(mean = NA, median = NA, var = NA)
  )

  # Check all summary statistics
  phylo <- createPhyloCRBD(100, 1, 0.5)
  df.nodes <- createNodesDf(phylo)
  df.edges <- createEdgesDf(phylo, df.nodes)
  sumstat.bl <- sumStatsBranchLength(phylo, df.edges)
  expect_equal(sumstat.bl$height, getHeight(phylo))
  expect_equal(
    sumstat.bl$mean.int1 / sumstat.bl$mean.ext,
    sumstat.bl$mean.intext1
  )
  expect_equal(
    sumstat.bl$median.int1 / sumstat.bl$median.ext,
    sumstat.bl$median.intext1
  )
  expect_equal(
    sumstat.bl$var.int1 / sumstat.bl$var.ext,
    sumstat.bl$var.intext1
  )
  expect_equal(
    sumstat.bl$mean.int2 / sumstat.bl$mean.ext,
    sumstat.bl$mean.intext2
  )
  expect_equal(
    sumstat.bl$mean.int3 / sumstat.bl$mean.ext,
    sumstat.bl$mean.intext3
  )
  expect_equal(length(sumstat.bl), 25)
})

test_that("Ladders are identified correctly.", {
  phylo <- createPhyloCRBD(3, 1, 0.5)
  expect_true(isInLadder(phylo, 4)) # root for phylo of 3 tips is in ladder
  expect_false(isInLadder(phylo, 3)) # tip can't be in ladder
  expect_false(isInLadder(phylo, 5)) # root child isn't in ladder (has 2 child)
  expect_equal(writeLadders(phylo, 4, c()), c(4, -1, -1, -1, -1))

  phylo <- createPhyloCRBD(3, 1, 0.5)
  ladders <- writeLadders(phylo, 4, c())
  expect_equal(length(unlist(extractLadders(ladders))), 0)

  phylo <- createPhyloCRBD(500, 1, 0.5)
  ladders <- writeLadders(phylo, 501, c())
  ladders.extracted <- extractLadders(ladders)
  expect_equal(ladders.extracted, getLadders(phylo))
  expect_true(getInLadderNodesProp(phylo) <= 1)
  df.nodes <- createNodesDf(phylo)
  expect_true(getImbalancedNodesProp(df.nodes) <= 1)

  skip_if(length(ladders.extracted) == 0)
  expect_true(isInLadder(phylo, ladders[which(ladders != -1)][1]))
  expect_true(isInLadder(phylo, unlist(extractLadders(ladders))[1]))
  max.ladder <- max(sapply(extractLadders(ladders), length))
  expect_equal(getMaxLadder(phylo), max.ladder / numberTips(phylo))
})

test_that("Left and right children identified correctly.", {
  phylo <- createPhyloCRBD(3, 1, 0.5)
  expect_equal(getChildLeftRight(phylo, 3), list(left = NA, right = NA))
  expect_equal(getChildLeftRight(phylo, 4)$right, 5)
})

test_that("Topological summary statistics are computed correctly.", {
  phylo <- createPhyloCRBD(3, 1, 0.5)
  df.nodes <- createNodesDf(phylo)
  expect_equal(getWidthDepthRatio(df.nodes), 1)
  expect_equal(getMaxWidthDiff(df.nodes), 1)
  expect_equal(getSackin(df.nodes), 5)

  phylo <- createPhyloCRBD(50, 1, 0.5)
  df.nodes <- createNodesDf(phylo)
  sumstat <- sumStatsTopo(phylo, df.nodes)
  expect_equal(length(sumstat), 8) # 8 statistics
  expect_equal(is.na(sumstat$maxladder), sumstat$inladder == 0)
})

test_that("LTT summary statistics are computed correctly.", {

  # LTT coord names
  names <- c(
    "LTT_t1", "LTT_t2", "LTT_t3", "LTT_t4", "LTT_t5", "LTT_t6", "LTT_t7",
    "LTT_t8", "LTT_t9", "LTT_t10", "LTT_t11", "LTT_t12", "LTT_t13",
    "LTT_t14", "LTT_t15", "LTT_t16", "LTT_t17", "LTT_t18", "LTT_t19",
    "LTT_t20", "LTT_N1", "LTT_N2", "LTT_N3", "LTT_N4", "LTT_N5", "LTT_N6",
    "LTT_N7", "LTT_N8", "LTT_N9", "LTT_N10", "LTT_N11", "LTT_N12",
    "LTT_N13", "LTT_N14", "LTT_N15", "LTT_N16", "LTT_N17", "LTT_N18",
    "LTT_N19", "LTT_N20"
  )
  expect_equal(lttCoordsNames(), names)

  # LTT coords
  expect_error(getLttCoords(createPhyloCRBD(5, 1, 0.5))) # phylogeny too small
  phylo <- createPhyloCRBD(50, 1, 0.5)
  coord <- getLttCoords(phylo)
  expect_equal(length(coord), 40) # good size
  expect_equal(as.vector(unlist(coord[names[21:40]])), 1:20 * 2) # check lineage

  # LTT slopes
  coords <- as.data.frame(ape::ltt.plot.coords(phylo))
  coord.split <- splitLttCoords(coords)
  expect_true(all(sapply(coord.split, nrow) == 5))
  expect_equal(merge(coord.split[[1]], coord.split[[2]], all=TRUE)$N, 1:10)
  slope1 <- stats::lm(log(N) ~ time, coord.split[[1]])$coef[[2]]
  expect_equal(getLttSlopes(phylo)$LTT_slope1, slope1)
  expect_error(getLttSlopes(createPhyloCRBD(5,1,0)))

  sumstat <- sumStatsLtt(phylo)
  expect_equal(length(sumstat),51)
  expect_equal(sumstat$n_tips,numberTips(phylo))
})

test_that("All summary statistics are computed correctly.", {
  phylo <- createPhyloCRBD(100, 1, 0.5)
  sumstat <- sumStats(phylo)
  expect_equal(length(sumstat), 84)
})

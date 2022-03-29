test_that("Node dataframe is created correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0)
  df.nodes <- createNodesDf(phylo)
  expect_equal(nrow(df.nodes), 2 * 10 - 1) # expected number of rows
  expect_equal(ncol(df.nodes), 7) # expected number of columns
  expect_true(all(is.na(df.nodes))) # filled with NAs
})

test_that("Summary statistics are computed correctly.", {

  # Colless score
  phylo <- createPhyloCRBD(3, 1, 0)
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
  phylo <- createPhyloCRBD(10, 1, 0)
  df.nodes <- createNodesDf(phylo)

  # Check index
  df.nodes <- fillNodesIndex(df.nodes)
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
  expect_equal(df.nodes$dist[root], 0) # distance between root to root is null
  expect_equal(df.nodes$dist, castor::get_all_distances_to_root(phylo))

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
  df.nodes.new <- createNodesDf(phylo)
  df.nodes.new <- fillNodesAll(df.nodes.new, phylo)
  expect_true(all(df.nodes == df.nodes.new, na.rm = TRUE))
})

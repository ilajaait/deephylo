test_that("Create phylogeny correctly.", {
  phylo <- createPhyloCRBD(15, 1, 0)
  expect_equal(numberTips(phylo), 15)
  expect_error(createPhyloCRBD(10, 0, 1))
})

test_that("Compute phylogeny number of nodes, tips and edges correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0)
  expect_equal(numberTips(phylo), 10)
  expect_equal(numberNodesTotal(phylo), 19)
  expect_equal(numberNodesInternal(phylo), 9)
  expect_equal(numberEdges(phylo), numberNodesTotal(phylo) - 1)
})

test_that("Compute height correctly.", {
  phylo <- createPhyloCRBD(15, 1, 0)
  expect_true(getHeight(phylo) > 0)
})

test_that("Identify tips correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0)
  expect_equal(isTip(phylo,1:19), 1:19 <= 10)
})

test_that("Find root correctly.", {
  phylo <- createPhyloCRBD(10, 1, 0)
  expect_equal(getRoot(phylo), 11)
})

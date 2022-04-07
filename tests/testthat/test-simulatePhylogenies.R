test_that("Create phylogeny correctly.", {

  # CRBD
  phylo <- createPhyloCRBD(15, list(lambda = 1, mu = 0))
  expect_equal(numberTips(phylo), 15)
  expect_error(createPhyloCRBD(10, list(lambda = 0, mu = 1)))

  # BiSSE
  phylo <- createPhyloBiSSE(15, list(lambda0 = 1, q = 0.1))
  expect_equal(numberTips(phylo), 15)
  phylo <- createPhyloBiSSE(10, list(lambda0 = 1, q = 0.001))
  expect_true(all(phylo$node.state == 1) | all(phylo$node.state == 0))
  expect_true(all(phylo$tip.state == 1) | all(phylo$tip.state == 0))
  phylo <- createPhyloBiSSE(100, list(lambda0 = 1, q = 0.5))
  expect_false(all(phylo$node.state == 1) | all(phylo$node.state == 0))
  expect_false(all(phylo$tip.state == 1) | all(phylo$tip.state == 0))
})

test_that("Draw parameters correctly.", {

  # Phylo size
  expect_equal(drawPhyloSize(c(10, 10)), 10)
  size <- drawPhyloSize(c(10, 20))
  expect_true(size >= 10 & size <= 20)

  # Rate CRBD
  expect_equal(list(lambda = 0, mu = 0), drawParamCRBD(c(0, 0)))
  params <- drawParamCRBD(c(0, 1))
  expect_true(params$lambda > params$mu)
  expect_true(params$lambda > 0 & params$mu > 0)
  expect_true(params$lambda < 1 & params$mu < 1)
})


test_that("Simulate phylogeny correctly.", {

  # CRBD
  phylo.list <- simulatePhyloCRBD(10, c(20, 50), c(0.1, 1.0))
  expect_equal(length(phylo.list), 3)
  expect_equal(length(phylo.list$phylo), 10)
  expect_equal(nrow(phylo.list$params), 10)
  expect_equal(ncol(phylo.list$params), 2)
  expect_equal(names(phylo.list$params), c("lambda", "mu"))
  expect_true(all(!is.na(phylo.list$params)))

  # BiSSE
  phylo.list <- simulatePhyloBiSSE(10, c(20, 50), list(lambda0 = c(0.1, 1.0), q = c(0, 0.1)))
  expect_equal(length(phylo.list), 3)
  expect_equal(length(phylo.list$phylo), 10)
  expect_equal(nrow(phylo.list$params), 10)
  expect_equal(ncol(phylo.list$params), 2)
  expect_equal(names(phylo.list$params), c("lambda0", "q"))
  expect_true(all(!is.na(phylo.list$params)))
})

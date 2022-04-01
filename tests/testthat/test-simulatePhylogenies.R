test_that("Create phylogeny correctly.", {

  # CRBD
  phylo <- createPhyloCRBD(15, 1, 0)
  expect_equal(numberTips(phylo), 15)
  expect_error(createPhyloCRBD(10, 0, 1))

  # BiSSE
  phylo <- createPhyloBiSSE(15, 1, 0.1)
  expect_equal(numberTips(phylo), 15)
  phylo <- createPhyloBiSSE(10, 1, 0.0000001)
  expect_true(all(phylo$node.state == 1) | all(phylo$node.state == 0))
  expect_true(all(phylo$tip.state == 1) | all(phylo$tip.state == 0))
  phylo <- createPhyloBiSSE(100, 1, 0.5)
  expect_false(all(phylo$node.state == 1) | all(phylo$node.state == 0))
  expect_false(all(phylo$tip.state == 1) | all(phylo$tip.state == 0))
})


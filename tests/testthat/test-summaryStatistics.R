test_that("Creating node dataframe.", {
  phylo <- diversitree::trees(
    pars = c(1., 0), type = "bd", n = 1, max.taxa = 10
  )[[1]]
  df <- createNodesDf(phylo)
  expect_equal(nrow(df), 2 * 10 - 1) # expected number of rows
  expect_equal(ncol(df), 7) # expected number of columns
  expect_true(all(is.na(df))) # filled with NAs
})

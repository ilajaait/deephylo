set.seed(113)
phyloCRBD <- simulatePhyloCRBD(10, c(100, 1000), c(0.1, 1.0))
usethis::use_data(phyloCRBD, overwrite = TRUE)

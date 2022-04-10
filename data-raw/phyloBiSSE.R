set.seed(113)
phyloBiSSE <- simulatePhyloBiSSE(10, c(100, 1000), list(lambda0=c(0.1, 1.0),q=c(0.001,0.1)))
usethis::use_data(phyloBiSSE, overwrite = TRUE)

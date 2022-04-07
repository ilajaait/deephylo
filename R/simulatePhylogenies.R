#' Create one phylogeny under the CRBD model
#'
#' Speciation and extinction rates are constant through time and lineages.
#'
#' @param n.tips number of tips
#' @param params list
#'    * `$lambda`: speciation rate
#'    * `$mu`: extinction rate
#'
#' @return phylogeny (ape format)
#' @export
#'
#' @examples n.tips <- 10 # 10 extant species
#' lambda <- 1 # speciation rate
#' mu <- 0 # extinction rate
#' createPhyloCRBD(n.tips, list(lambda = lambda, mu = mu))
createPhyloCRBD <- function(n.tips, params) {
  lambda <- params$lambda
  mu <- params$mu
  lambda > mu || stop("Speciation rate inferior to extinction rate.")
  diversitree::trees(
    pars = c(lambda, mu), type = "bd", max.taxa = n.tips, n = 1
  )[[1]]
}

#' Create one phylogeny under the BiSSE model
#'
#' The BiSSE diversification model has 6 parameters:
#'    * \eqn{\lambda_0} speciation rate of state 0
#'    * \eqn{\lambda_1} speciation rate of state 1
#'    * \eqn{\mu_0} extinction rate of state 0
#'    * \eqn{\mu_1} extinction rate of state 1
#'    * \eqn{q_{01}} transition rate from state 0 to 1
#'    * \eqn{q_{10}} transition rate from state 1 to 0
#'
#' We assume 4 constraints:
#'    * \eqn{\lambda_1 = 2\lambda_0}
#'    * \eqn{\mu_0 = 0}
#'    * \eqn{\mu_1 = 0}
#'    * \eqn{q_{10} = q_{01}}
#'
#' @param n.tips number of tips
#' @param params list
#'    * `$lambda0`: speciation rate of state 0
#'    * `$q`: transition rate
#'
#' @return phylogeny (ape format)
#' @export
#'
#' @examples n.tips <- 10 # 10 extant species
#' lambda0 <- 1 # speciation rate
#' q <- 0.1 # transition rate
#' createPhyloBiSSE(n.tips, list(lambda0 = lambda0, q = q))
createPhyloBiSSE <- function(n.tips, params) {
  lambda0 <- params$lambda0
  q <- params$q
  diversitree::tree.bisse(c(lambda0, 2 * lambda0, 0, 0, q, q), max.taxa = n.tips)
}

#' Simulate phylogenies under the CRBD model
#'
#' Generate a list of phylogenies under the Constant Rate Birth Death (CRBD)
#' model. Can check and save summary statistics of each phylogeny by activating
#' `sumstat_check`. Moreover true parameters of each phylogeny are also
#' returned.
#'
#' @param n.phylo number of phylogenies to simulate
#' @param tips.lim phylogeny size range
#' @param lambda.lim speciation rate range
#' @param sumstat_check whether to check that summary statistics can be computed
#'    for each phylogeny and save them (logical)
#' @param verbose whether to print progress (logical)
#'
#' @return list
#'    * `$phylo`: list of simulated phylogenies (ape format)
#'    * `$params`: dataframe of phylogeny parameters
#'    * `$sumstat`: dataframe of phylogeny summary statistics
#'    (empty if `ss_check = FALSE`)
#'
#' @export
#'
#' @examples simulatePhyloCRBD(10, c(50, 100), c(0.1, 1.)) # simulate 10 phylo
simulatePhyloCRBD <- function(n.phylo, tips.lim, lambda.lim,
                              sumstat_check = T, verbose = F) {

  # Tests
  length(tips.lim) == 2 || stop("tips.lim should be of length 2.")
  length(lambda.lim) == 2 || stop("lambda.lim should be of length 2.")

  # Simulate
  out <- simulatePhylo(createPhyloCRBD, drawParamCRBD, n.phylo, tips.lim,
                       lambda.lim, sumstat_check, verbose)
  names(out$params) <- c("lambda", "mu")
  out
}

#' Simulate phylogenies under the BiSSE model
#'
#' Generate a list of phylogenies under the Binary State Speciation and Extinction
#' (BiSSE) model. Can check and save summary statistics of each phylogeny by
#' activating `sumstat_check`. Moreover true parameters of each phylogeny are also
#' returned.
#'
#' @param n.phylo number of phylogenies to simulate
#' @param tips.lim phylogeny size range
#' @param params.lim parameter ranges
#' @param sumstat_check whether to check that summary statistics can be computed
#'    for each phylogeny and save them (logical)
#' @param verbose whether to print progress (logical)
#'
#' @return list
#'    * `$phylo`: list of simulated phylogenies (ape format)
#'    * `$params`: dataframe of phylogeny parameters
#'    * `$sumstat`: dataframe of phylogeny summary statistics
#'    (empty if `ss_check = FALSE`)
#'
#' @export
#'
#' @examples simulatePhyloBiSSE(10, c(50, 100), list(lambda0 = c(0.1, 1.0), q = c(0, 0.1)))
simulatePhyloBiSSE <- function(n.phylo, tips.lim, params.lim,
                               sumstat_check = T, verbose = F) {

  # Tests
  length(tips.lim) == 2 || stop("tips.lim should be of length 2.")
  length(params.lim$lambda0) == 2 || stop("lambda0.lim should be of length 2.")
  length(params.lim$q) == 2 || stop("q.lim should be of length 2.")

  # Simulate
  out <- simulatePhylo(createPhyloBiSSE, drawParamBiSSE, n.phylo,
                       tips.lim, params.lim, sumstat_check, verbose)
  names(out$params) <- c("lambda0", "q")
  out
}

#' Backbone to simulate phylogenies
#'
#' Internal function used to simulated phylogenies. Called by \code{\link{simulatePhyloCRBD}}
#' and \code{\link{simulatePhyloBiSSE}}.
#'
#' @param createPhylo function used to create one phylogeny under a diversification model
#' @param drawParam function to draw parameters of the diversification model
#' @param n.phylo number of phylogenies to generate (int)
#' @param tips.lim phylogeny size boundaries
#' @param params.lim parameter boundaries
#' @param sumstat_check logical
#' @param verbose logical
simulatePhylo <- function(createPhylo, drawParam, n.phylo, tips.lim, params.lim,
                          sumstat_check, verbose) {
  # Set up
  phylo.list <- list() # store simulated phylogenies
  params.df <- data.frame(param1 = rep(NA, n.phylo), param2 = rep(NA, n.phylo))
  sumstat.list <- list()

  # Simulation loop
  while (length(phylo.list) < n.phylo) {
    params <- drawParam(params.lim)
    n.tips <- drawPhyloSize(tips.lim)
    phylo <- createPhylo(n.tips, params)
    if (sumstat_check) {
      sumstat <- sumStats(phylo)
      if (all(!is.na(sumstat))) {
        phylo.list <- append(phylo.list, list(phylo))
        sumstat.list <- append(sumstat.list, list(sumstat))
        params.df[length(phylo.list), ] <- params
      }
    } else {
      phylo.list <- append(phylo.list, list(phylo))
      params.df[length(phylo.list), ] <- params
    }
    if (verbose) {
      svMisc::progress(length(phylo.list), n.phylo,
                       progress.bar = TRUE,
                       init = (length(phylo.list) == 1)
      )
    }
  }

  # Output
  sumstat.df <- as.data.frame(do.call(rbind, sumstat.list))
  list(phylo = phylo.list, params = params.df, sumstat = sumstat.df)
}

#' Draw parameters for the CRBD model
#'
#' To ensure that \eqn{\lambda >> \mu}, parameters are drawn as follow:
#'    * \eqn{\lambda} in (`min(lambda.lim)`, `max(lambda.lim)`)
#'    * \eqn{\mu} in (0, \eqn{0.9 \lambda})
#'
#' @param lambda.lim speciation rate boundaries (vector)
#'
#' @return list of parameters
drawParamCRBD <- function(lambda.lim) {
  lambda <- stats::runif(1, lambda.lim[1], lambda.lim[2])
  mu <- stats::runif(1, 0, 0.9 * lambda)
  list(lambda = lambda, mu = mu)
}

#' Draw parameters for the BiSSE model
#'
#' Draw speciation rate of state 0 (lambda0) and transition rate randomly.
#' Other parameters are constrained.
#'
#' @param params.lim parameters boundaries (list of vector)
#'    * `lambda0`: boundaries of speciation rate of state 0
#'    * `q`: boundaries of transition rate
#'
#' @return list of parameters
drawParamBiSSE <- function(params.lim) {
  all(names(params.lim) == c("lambda0", "q")) || stop("Wrong parameter names.")
  lambda0.lim <- params.lim$lambda0
  q.lim <- params.lim$q
  lambda0 <- stats::runif(1, lambda0.lim[[1]], lambda0.lim[[2]])
  q <- stats::runif(1, q.lim[[1]], q.lim[[2]])
  list(lambda0 = lambda0, q = q)
}

#' Draw size of a phylogeny
#'
#' @param tips.lim size boundaries (vector)
#'
#' @return int
drawPhyloSize <- function(tips.lim) {
  round(stats::runif(1, tips.lim[1], tips.lim[2]))
}

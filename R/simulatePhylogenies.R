#' Create one phylogeny under the Constant Rate Birth-Death (CRBD) model
#'
#' @param n.tips number of tips
#' @param lambda speciation rate
#' @param mu extinction rate
#'
#' @return phylogeny (ape format)
#' @export
#'
#' @examples n.tips <- 10 # 10 extant species
#' lambda <- 1 # speciation rate
#' mu <- 0 # extinction rate
#' createPhyloCRBD(n.tips, lambda, mu)
createPhyloCRBD <- function(n.tips, lambda, mu) {
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
#' @param lambda0 speciation rate of state 0
#' @param q transition rate
#'
#' @return phylogeny (ape format)
#' @export
#'
#' @examples n.tips <- 10 # 10 extant species
#' lambda0 <- 1 # speciation rate
#' q <- 0.1 # transition rate
#' createPhyloBiSSE(n.tips, lambda0, q)
createPhyloBiSSE <- function(n.tips, lambda0, q) {
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

  # Set up
  phylo.list <- list() # store simulated phylogenies
  params.df <- data.frame(lambda = rep(NA, n.phylo), mu = rep(NA, n.phylo))
  sumstat.list <- list()

  # Simulation loop
  while (length(phylo.list) < n.phylo) {
    params <- drawParamCRBD(lambda.lim)
    n.tips <- drawPhyloSize(tips.lim)
    phylo <- createPhyloCRBD(n.tips, params$lambda, params$mu)
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

#' Draw size of a phylogeny
#'
#' @param tips.lim size boundaries (vector)
#'
#' @return int
drawPhyloSize <- function(tips.lim) {
  round(stats::runif(1, tips.lim[1], tips.lim[2]))
}

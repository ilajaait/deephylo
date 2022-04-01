#' Create a phylogeny under the Constant Rate Birth-Death (CRBD) model
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

#' Create a phylogeny under the BiSSE model
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
  diversitree::tree.bisse(c(lambda0,2*lambda0,0,0,q,q), max.taxa = n.tips)
}

#' Algebraic degree of semidefinite programming
#'
#' A list of the values of the algebraic degree of semidefinite programming
#' delta(m,n,r) for \code{n=2,3,...,14}.
#' The values are given as strings to avoid rounding errors.
#'
#' @format A list of 13 character matrices; the rows and columns corresponds to
#' m and r, respectively.
#' 
#' @section See also:
#' \code{\link[symconivol]{alg_deg}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @source The degrees have been computed with Maple using a formula given in
#' \url{https://arxiv.org/abs/1509.05241}, \url{https://projecteuclid.org/euclid.kmj/1478073765}.
#' 
"alg_deg_data"



#' Ratios of index probabilities
#'
#' A list of sample counts of a Bernoulli variable with (unnormalized) success
#' and failure probabilities given by Prob\{ind=r\} and Prob\{ind=r+1\}.
#' 
#' The values of the parameters beta, n, r are specified in the names of the list
#' elements, which are of the form \code{beta=_,n=_,r=_} with _ replaced by
#' values. The list contains all combinations of \code{beta=1,2,4}, \code{n=4,5,...,10},
#' and \code{r} between \code{floor((n+1)/2)} and \code{n-2}. The remaining values
#' can be found through the relation Prob\{ind=r\}=Prob\{ind=n-r\}.
#' 
#' @format A list of two-element integer vectors.
#' 
#' @section See also:
#' \code{\link[symconivol]{constr_eigval}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @source The samples have been found with the HMC sampler Stan with a model
#' that simulates the eigenvalues of a Gaussian orthogonal/unitary/symplectic
#' ensemble with index constrained to \{r,r+1\}. The code below shows how this data
#' can be generated.
#' 
#' @examples 
#' \dontrun{
#' library(rstan)
#' 
#' ind_prob <- list()
#' nmin <- 4
#' nmax <- 10
#' filename <- "tmp.stan"
#' N <- 1e6
#' warmup <- 1e3
#' for (n in nmin:nmax) {
#'     for (beta in c(1,2,4)) {
#'         for (r in floor((n+1)/2):(n-2)) {
#'             index <- str_c("beta=",beta,",n=",n,",r=",r)
#'             M <- constr_eigval(beta, n, r, r+1, filename=filename)
#'             stan_samp <- stan( file = filename, data = M$data,
#'                                chains = 1, warmup = warmup,
#'                                iter = N+warmup, cores = 2, refresh = 1e4 )
#'             ef <- rstan::extract(stan_samp)$ef
#'             indprob[[index]] <- c(length(which(ef<0)),length(which(ef>0)))
#'         }
#'     }
#' }
#' file.remove(filename)
#' }
#' 
"ind_prob"



#' Values of estimated limit curve of dimension normalized curvature measures
#'
#' A table of function values of the estimated limit curve of
#' dimension normalized curvature measures mu.
#' 
#' @format A two-column table of function values of mu; the columns corresponds to
#' \code{rho} and \code{kappa=mu(rho)}.
#' 
#' @section See also:
#' \code{\link[symconivol]{mu}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @source The function values have been computed with an estimated value of C=0.2.
#' The code below shows how this data can be generated.
#' 
#' @examples 
#' \dontrun{
#' library(tidyverse)
#' 
#' L <- leigh(1e5)
#' R <- rate(1e5)
#' rho <- seq(0,1,length.out=1e4)
#' Lkappa <- L$lkup_kappa(rho)
#' Rkappa <- R$lkup_R(rho)
#' kappa <- seq(0,1,length.out=1e4)
#' C <- 0.2
#' rhomin <- sapply(kappa, function(x)
#'     return( max(1-sqrt(1-x), min(sqrt(x),rho[which.min(
#'         (x-Lkappa)^2 / (1-abs(2*rho-1)) + 4*C*Rkappa )] ))))
#' 
#' mu_data <- tibble(kappa=kappa,rho=rhomin)
#' }
#' 
"mu_data"



#' Reconstructed values of index constrained curvature measures
#'
#' A list of reconstructed values of index constrained curvature measures;
#' the constraints being of the form r<=ind(x)<=r+s.
#' 
#' The values have been found by the following steps:
#' \itemize{
#'     \item sample eigenvalues from the index constrained Gaussian
#'           orthogonal/unitary/symplectic ensemble using the HMC sampler Stan,
#'     \item convert the samples samples from the corresponding bivariate
#'           chi-bar-squared distribution,
#'     \item reconstruct the weights of the bivariate chi-bar-squared distribution
#'           by running the EM algorithm for 100 steps.
#' }
#' 
#' @format A list of lists of vectors containing the reconstructed weights;
#' the element names of the outer list are of the form "beta=_,n=", while the
#' the element names of the inner lists are of the form "r=_,s=", the "_" of course
#' replaced by corresponding values.
#' 
#' @section See also:
#' \code{\link[symconivol]{constr_eigval}},
#' \code{\link[symconivol]{constr_eigval_to_bcbsq}},
#' \code{\link[symconivol]{estim_em_cm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @source The values have been computed using the HMC sampler stan and functions
#' from the symconivol package. The code below shows how this data can be generated.
#' 
#' @examples 
#' \dontrun{
#' library(tidyverse)
#' library(rstan)
#' 
#' warmup <- 1e3
#' N <- 1e5
#' filename_model <- "tmp.stan"
#' 
#' nmin <- 3
#' nmax <- 10
#' phi_ind <- list()
#' for (beta in c(1,2,4)) {
#'     for (n in nmin:nmax) {
#'         pat <- pat_bnd(beta,n)
#'         index_out <- str_c("beta=",beta,",n=",n)
#'         weights <- list()
#'         for (r in floor((n+1)/2):(n-1)) {
#'             for (s in 0:(n-1-r)) {
#'                 index_in <- str_c("r=",r,",s=",s)
#'                 np <- r
#'                 nf <- s
#'                 nn <- n-r-s
#'                 M <- constr_eigval(pos=(np>0), free=(nf>0), neg=(nn>0),
#'                               filename=filename_model, overwrite=TRUE)
#'                 stan_samp <- stan( file = filename_model, data = M$data,
#'                                    chains = 1, warmup = warmup, iter = N+warmup,
#'                                    cores = 2, refresh = 1e4 )
#'                 samp <- list()
#'                 if (np>0) samp$ep <- rstan::extract(stan_samp)$ep
#'                 if (nf>0) samp$ef <- rstan::extract(stan_samp)$ef
#'                 if (nn>0) samp$en <- rstan::extract(stan_samp)$en
#' 
#'                 m_samp <- constr_eigval_to_bcbsq(pos=(np>0), free=(nf>0), neg=(nn>0),
#'                                                  samp=samp)
#' 
#'                 em <- estim_em_cm(d=pat$d, low=pat$k_low(r), upp=pat$k_upp(r+s),
#'                                   m_samp=m_samp, N=100)
#'                 weights[[index_in]] <- em[101,]
#'             }
#'         }
#'         phi_ind[[index_out]] <- weights
#'     }
#' }
#' file.remove(filename_model)
#' }
#' 
"phi_ind"



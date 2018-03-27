#' Converting samples of index constrained eigenvalues to samples of bivariate chi-bar-squared distribution
#' 
#' \code{constr_eigval_to_bcbsq} converts samples of the index constrained
#' eigenvalue distribution to samples of the corresponding bivariate
#' chi-bar-squared distribution.
#' 
#' @param samp a list containing (a selection of) the elements \'ep\', \'ef\', \'en\',
#'             which are matrices (of the same length) whose rows give the sampled eigenvalues
#' 
#' @return The output of \code{constr_eigval_to_bcbsq} is a two-column matrices
#'         whose rows correspond to the samples from the bivariate
#'         chi-bar-squared distribution.
#' 
#' @section See also:
#' \code{\link[symconivol]{constr_eigval}}, 
#' \code{\link[symconivol]{prepare_em_cm}}, 
#' \code{\link[symconivol]{estim_em_cm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' library(tidyverse)
#' library(rstan)
#' 
#' filename <- "tmp.stan"
#' M <- constr_eigval( beta=2, n=12, ind_low=8, ind_upp=9, filename=filename )
#' stan_samp <- stan( file = filename, data = M$data,
#'                    chains = 1, warmup = 1e3, iter = 1e5, cores = 2, refresh = 1e4 )
#' file.remove(filename)
#' 
#' m_samp <- constr_eigval_to_bcbsq( rstan::extract(stan_samp) )
#' 
#' @export
#'
constr_eigval_to_bcbsq <- function(samp) {
    neg  <- "en" %in% names(samp)
    free <- "ef" %in% names(samp)
    pos  <- "ep" %in% names(samp)
    
    if ( !pos & !free & !neg ) stop("\n Empty model (list elements have to be named \'ep\', \'ef\', \'en\' for positive, free, negative eigenvalues, respectively.")

    if (pos)        n <- dim(samp$ep)[1]
    else if (free)  n <- dim(samp$ef)[1]
    else if (neg)   n <- dim(samp$en)[1]
    
    possq <- vector("numeric",n)
    negsq <- vector("numeric",n)
    
    if (pos) possq <- possq + rowSums(samp$ep^2)
    if (neg) negsq <- negsq + rowSums(samp$en^2)
    if (free) {
        possq <- possq + rowSums(pmax(samp$ef,0)^2)
        negsq <- negsq + rowSums(pmin(samp$ef,0)^2)
    }
    return(cbind(possq,negsq))
}


#' Evaluate bivariate chi-bar-squared samples for maximum likelihood estimation
#' 
#' \code{prepare_em_cm} takes a two-column matrix whose rows form
#' iid samples from a bivariate chi-bar-squared distribution and
#' prepares the data used in maximum likelihood estimation.
#'
#' This function works pretty much exactly as \code{prepare_em} from the
#' \code{conivol} package, the only difference being that the "boundary
#' cases" \code{k==0,n} do not have to be considered/are ignored.
#' In the general case this is not needed, but for the curvature
#' measures this is a useful feature.
#' 
#' @param d the dimension of the bivariate chi-bar squared distribution.
#' @param low lower bound for \code{k}; has to be \code{>0}
#' @param upp upper bound for \code{k}; has to be \code{<d}
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#' 
#' @return The output of \code{prepare_em_cm} is \code{(low-upp+1)} row matrix whose
#'         \code{k}th row contains the products of the density values of the chi_k^2
#'         and chi_(d-k)^2 distributions evaluated in the sample points;
#'         the row-form of the matrix is more convenient for the computations.
#' 
#' @section See also:
#' \code{\link[conivol]{prepare_em}}, 
#' \code{\link[symconivol]{constr_eigval}}, 
#' \code{\link[symconivol]{constr_eigval_to_bcbsq}}, 
#' \code{\link[symconivol]{estim_em_cm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' CM <- curv_meas_exact(4,3)$A[,2]
#' CM <- CM/sum(CM)
#' 
#' m_samp <- conivol::rbichibarsq(1e5,CM)
#' 
#' str( prepare_em_cm( 15, 1, 9, m_samp ))
#' 
#' @export
#'
prepare_em_cm <- function(d, low, upp, m_samp) {
    # low>0, upp<d
    return( apply( m_samp, 1, function(x){dchisq(x[1],low:upp)*dchisq(x[2],(d-low):(d-upp))} ) )
}


.create_mosek_input_em_cm <- function(const) {
    m <- length(const)
    
    mos_inp <- list()
    mos_inp$sense <- "max"
    
    # setting optimizer
    mos_inp$c <- vector("numeric",m)
    opro <- matrix(list(), nrow=5, ncol=m)
    rownames(opro) <- c("type","j","f","g","h")
    for (i in 1:m) {
        opro[ ,i] <- list("LOG", i, const[i], 1.0, 0.0)
    }
    
    mos_inp$scopt <- list(opro=opro)
    
    # variable constraints
    blx <- rep(0, m)            # v_i >= 0
    bux <- rep(1, m)            # v_i <= 1
    mos_inp$bx <- rbind(blx, bux)
    
    # constraint matrix:
    mos_inp$A <- Matrix::Matrix(rep(1,m),1,m)
    
    # constraint rhs:
    mos_inp$bc <- rbind(1,1)
    
    return( mos_inp )
}


.update_mosek_input_em_cm <- function(mos_inp,const) {
    mos_inp$scopt$opro[3, ] <- const
    return(mos_inp)
}


#' Estimating curvature measures from bivariate chi-bar-squared data using EM algorithm
#' 
#' \code{estim_em_cm} produces EM-type iterates from a two-column
#' matrix whose rows form iid samples from a bivariate chi-bar-squared
#' distribution.
#'
#' The sequence of iterates may or may not converge
#' to the maximum likelihood estimate of the mixing weights of the distribution.
#' Log-concavity of the weights is enforced by projecting the logarithms
#' onto the cone of log-concave sequences; this can be turned off by setting
#' \code{no_of_lcc_projections=0}.
#' This function is adapted from \code{estim_em} from the \code{conivol}
#' package, the difference being that the support of the weights is strictly
#' between the boundary cases. It is simplified in that the initial estimate
#' is always the uniform distribution, and the parity equation,
#' which does not hold for curvature measures, will not be enforced.
#' 
#' @param d the dimension of the bivariate chi-bar-squared distribution.
#' @param low lower bound for \code{k}; has to be \code{>0}
#' @param upp upper bound for \code{k}; has to be \code{<d}
#' @param m_samp two-column matrix whose rows from iid samples from a bivariate
#'               chi-bar-squared distribution.
#' @param N the number of iterates that shall be produced.
#' @param no_of_lcc_projections number of projections on the log-concavity cone
#' @param data output of \code{prepare_em_cm(d, low, upp, m_samp)}; this can be called
#'              outside and passed as input to avoid re-executing this
#'              potentially time-consuming step.
#' 
#' @return The output of \code{estim_em_cm} is a list of an \code{(N+1)}-by-\code{(upp-low+1)}
#'         matrix whose rows constitute EM-type iterates, which may or may not
#'         converge to the maximum likelihood estimate of the mixing weights of
#'         the bivariate chi-bar-squared distribution, and the corresponding values
#'         of the log-likelihood function.
#' 
#' @section See also:
#' \code{\link[conivol]{estim_em}}, 
#' \code{\link[symconivol]{constr_eigval}}, 
#' \code{\link[symconivol]{constr_eigval_to_bcbsq}}, 
#' \code{\link[symconivol]{prepare_em_cm}},
#' \code{\link[symconivol]{indnorm_to_unnorm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' CM <- curv_meas_exact(4,3)$A[,2]
#' CM <- CM/sum(CM)
#' 
#' m_samp <- conivol::rbichibarsq(1e5,CM)
#' 
#' d   <- 15
#' low <- 1
#' upp <- 9
#' est <- estim_em_cm( d, low, upp, m_samp )
#' 
#' plot(1+low:upp, CM[1+low:upp])
#' lines(1+low:upp, CM[1+low:upp], col="red")
#' lines(1+low:upp, est[1,])
#' lines(1+low:upp, est[5,])
#' lines(1+low:upp, est[10,])
#' lines(1+low:upp, est[21,])
#' 
#' @export
#'
estim_em_cm <- function(d, low, upp, m_samp, N=20, no_of_lcc_projections=1, data=NULL) {
    if (!requireNamespace("Rmosek", quietly = TRUE))
        stop( paste0("\n Could not find package 'Rmosek'.",
                     "\n See the help entries for more information.") )
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("\n Could not find package 'Matrix'.")
    #################################
    # set return messages from MOSEK (verbose=0 enforces silent mode)
    opts <- list(verbose=0)
    #################################

    # find the values of the chi-squared densities at the sample points
    if (is.null(data))
        data <- symconivol::prepare_em_cm(d, low, upp, m_samp)
    n <- dim(m_samp)[1]
    
    out_iterates <- matrix(0,N+1,upp-low+1)

    # set the starting point for EM
    v <- rep(1,upp-low+1)/(upp-low+1)
    out_iterates[1, ] <- v

    # prepare Mosek inputs
    mos_inp <- .create_mosek_input_em_cm(rep(0,upp-low+1))
    
    # prepare Mosek inputs for log-concavity enforcing
    A_lcc <- matrix(0,upp-low+1,upp-low-1)
    diag(A_lcc) <- 1
    if (upp-low==2) {
        A_lcc[2,1] <- -2
        A_lcc[3,1] <- 1
    } else {
        diag(A_lcc[2:(upp-low),]) <- -2
        diag(A_lcc[3:(upp-low+1),]) <- 1
    }
    mos_inp_lcc <- conivol:::.create_mosek_input_polyh_pol(A_lcc, rep(0,upp-low+1), 0)
    
    for (i in 1:N) {
        denom <- colSums( data * v )
        const <- rowSums( sweep( 1/n * data * v , MARGIN=2, denom, "/") )
        mos_inp <- .update_mosek_input_em_cm(mos_inp,const)
        mos_out <- Rmosek::mosek(mos_inp, opts)
        v <- mos_out$sol$itr$xx
        if (no_of_lcc_projections>0)
            for (i_lcc in (1:no_of_lcc_projections)) {
                mos_inp_lcc <- conivol:::.update_mosek_input_polyh_pol(mos_inp_lcc, log(v))
                mos_out <- Rmosek::mosek(mos_inp_lcc, opts)
                v <- exp(mos_out$sol$itr$xx[1:(upp-low+1)])
                v <- v/sum(v)
            }
        out_iterates[i+1, ] <- v
    }
    
    return( iterates=out_iterates )
}












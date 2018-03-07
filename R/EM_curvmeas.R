#' Converting samples of index constrained eigenvalues to samples of bivariate-chi-bar-squared distribution
#' 
#' \code{constr_eigval_to_bcbsq}
#' 
#' @param pos parameter description
#' @param free parameter description
#' @param neg parameter description
#' @param samp parameter description
#' @param filename parameter description
#' 
#' @return returns what
#' 
#' @section See also:
#' \code{\link[symconivol]{constr_eigval}}, 
#' \code{\link[symconivol]{prepare_em_cm}}, 
#' \code{\link[symconivol]{estim_em_cm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' # some example code
#' 
#' @export
#'
constr_eigval_to_bcbsq <- function(pos, free, neg, samp=NA, filename=NA) {
    if ( !pos & !free & !neg ) stop("\n Empty model.")
    if ( (is.na(samp) && is.na(filename)) || (!is.na(samp) && !is.na(filename)) )
        stop("\n Must give either sample matrix or sample file.")
    if ( !is.na(filename) && !file.exists(filename) )
        stop("\n File with given filename does not exist.")

    if ( !is.na(filename) ) {
        e <- new.env()
        samp_str <- load(file=filename, envir=e)
        samp <- get(samp_str, envir=e)
    }
    
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
#' \code{prepare_em_cm}
#' 
#' @param d parameter description
#' @param low parameter description
#' @param upp parameter description
#' @param m_samp parameter description
#' 
#' @return returns what
#' 
#' @section See also:
#' \code{\link[symconivol]{constr_eigval}}, 
#' \code{\link[symconivol]{constr_eigval_to_bcbsq}}, 
#' \code{\link[symconivol]{estim_em_cm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' # some example code
#' 
#' @export
#'
prepare_em_cm <- function(d, low, upp, m_samp) {
    # low>0, upp<d
    return( list( d=d, low=low, upp=upp, n=dim(m_samp)[1],
                  dens=apply( m_samp, 1, 
                   function(x){dchisq(x[1],low:upp)*dchisq(x[2],(d-low):(d-upp))}
            ) ) )
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
#' \code{estim_em_cm}
#' 
#' @param d parameter description
#' @param low parameter description
#' @param upp parameter description
#' @param m_samp parameter description
#' @param N parameter description
#' @param no_off_lcc_projections parameter description
#' @param data parameter description
#' 
#' @return returns what
#' 
#' @section See also:
#' \code{\link[symconivol]{constr_eigval}}, 
#' \code{\link[symconivol]{constr_eigval_to_bcbsq}}, 
#' \code{\link[symconivol]{prepare_em_cm}},
#' \code{\link[symconivol]{indnorm_to_unnorm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' # some example code
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
        denom <- colSums( data$dens * v )
        const <- rowSums( sweep( 1/data$n * data$dens * v , MARGIN=2, denom, "/") )
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












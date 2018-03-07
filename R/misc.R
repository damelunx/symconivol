#' Exact curvature measures of symmetric cones
#' 
#' \code{curv_meas_exact} returns the exact values of the curvature measures of
#' symmetric cones, which are known (\code{n==1,2,3}).
#' 
#' The known curvature measures are elements in the ring \code{Q[sqrt(2),1/pi]},
#' so the exact values can be given in terms of integers corresponding to a 
#' common denominator and corresponding integer matrices for the coefficients
#' of the natural expansion in \code{1, sqrt(2), 1/pi, sqrt(2)/pi}. These matrices
#' are returned by this function, along with the denominator and the combined
#' matrix of curvature measures.
#' 
#' @param beta Dyson index specifying the underlying (skew-) field:
#'             \describe{
#'               \item{\code{beta==1}:}{real numbers}
#'               \item{\code{beta==2}:}{complex numbers}
#'               \item{\code{beta==4}:}{quaternion numbers}
#'             }
#' @param n size of matrix.
#' 
#' @return The output of \code{curv_meas_exact} is a list of six elements:
#'         \itemize{
#'           \item \code{A}: the combined matrix of curvature measures; \code{A}
#'               is given in terms of the other parameters by
#'               \code{(A_const + A_sqrt2*sqrt(2) + A_piinv/pi + A_sqrt2_pi*sqrt(2)/pi)/denom}
#'           \item \code{A_const}: integer matrix for the constant term
#'           \item \code{A_sqrt2}: integer matrix for the \code{sqrt(2)} term
#'           \item \code{A_piinv}: integer matrix for the \code{1/pi} term
#'           \item \code{A_sqrt2_pi}: integer matrix for the \code{sqrt(2)/pi} term
#'           \item \code{denom}: common denominator of all terms
#'         }
#' 
#' @section See also:
#' \code{\link[symconivol]{alg_deg}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' # considering the case of 3x3 complex unitary matrices
#' CM <- curv_meas_exact(2,3)
#' 
#' # sum of intrinsic volumes is equal to one
#' sum( CM$A )
#' 
#' # sum of even (and odd) index intrinsic volumes is 1/2
#' sum( CM$A %*% rep_len(c(1,0),dim(CM$A)[2]) )
#' 
#' # A is given by combining the remaining matrices and the denominator
#' norm( CM$A - ( CM$A_const + CM$A_sqrt2*sqrt(2) + CM$A_piinv/pi + CM$A_sqrt2_pi*sqrt(2)/pi )/CM$denom )
#' 
#' @export
#'
curv_meas_exact <- function(beta,n) {
    if (!(beta %in% c(1,2,4))) {
        stop("Dyson index beta must be 1, 2, or 4.")
    } else if (n<=0) {
        stop("n must be positive.")
    } else if (n>3) {
        stop("Exact values unknown for n>3.")
    }
    d <- n+beta*choose(n,2)
    A <- matrix(0,d+1,n+1)
    rownames(A) <- paste0("k=",0:d)
    colnames(A) <- paste0("r=",0:n)
    A_const    <- A
    A_sqrt2    <- A
    A_piinv    <- A
    A_sqrt2_pi <- A
    denom <- 0
    if (n==1) {
        denom <- 2
        A_const[1,1] <- 1
    } else if (n==2) {
        if (beta==1) {
            denom <- 4
            A_const[1,1] <- 2
            A_sqrt2[1,1] <- -1
            A_sqrt2[2,2] <- 1
        } else if (beta==2) {
            denom <- 4
            A_const[1,1] <- 1
            A_const[2,2] <- 1
            A_piinv[1,1] <- -2
            A_piinv[3,2] <- 2
        } else if (beta==4) {
            denom <- 24
            A_const[1,1] <- 6
            A_const[2,2] <- 3
            A_const[4,2] <- 3
            A_piinv[1,1] <- -16
            A_piinv[3,2] <- 16
        }
        
    } else if (n==3) {
        if (beta==1) {
            denom <- 4
            A_const[1,1] <- 1
            A_const[2,2] <- 2
            A_const[4,2] <- -1
            A_sqrt2_pi[1,1] <- -2
            A_sqrt2[2,2] <- -1
            A_sqrt2_pi[3,2] <- 2
            A_sqrt2[4,2] <- 1
        } else if (beta==2) {
            denom <- 16
            A_const[1,1]  <- 2
            A_const[2,2]  <- 3
            A_const[6,2]  <- 3
            A_piinv[1,1]  <- -6
            A_piinv[2,2]  <- -6
            A_piinv[3,2]  <- 8
            A_piinv[4,2]  <- 8
            A_piinv[5,2]  <- 4
            A_piinv[6,2]  <- -8
        } else if (beta==4) {
            denom <- 960
            A_const[1,1]  <- 120
            A_const[2,2]  <- 105
            A_const[4,2]  <- 60
            A_const[6,2]  <- 90
            A_const[8,2]  <- -60
            A_const[10,2] <- 165
            A_piinv[1,1]  <- -376
            A_piinv[2,2]  <- -280
            A_piinv[3,2]  <- 192
            A_piinv[4,2]  <- 160
            A_piinv[5,2]  <- 384
            A_piinv[7,2]  <- 152
            A_piinv[8,2]  <- 256
            A_piinv[9,2]  <- 24
            A_piinv[10,2] <- -512
        }
    }
    A_const <- A_const + A_const[1+d:0,1+n:0]
    A_sqrt2 <- A_sqrt2 + A_sqrt2[1+d:0,1+n:0]
    A_piinv <- A_piinv + A_piinv[1+d:0,1+n:0]
    A_sqrt2_pi <- A_sqrt2_pi + A_sqrt2_pi[1+d:0,1+n:0]
    
    A <- ( A_const + A_sqrt2*sqrt(2) + A_piinv/pi + A_sqrt2_pi*sqrt(2)/pi )/denom
    return( list( A=A, A_const=A_const, A_sqrt2=A_sqrt2, A_piinv=A_piinv,
                  A_sqrt2_pi=A_sqrt2_pi, denom=denom ) )
}


#' Pataki bounds
#' 
#' \code{pat_bnd} provides the Pataki inequalities for given \code{beta} and \code{n}.
#' 
#' @param beta Dyson index specifying the underlying (skew-) field:
#'             \describe{
#'               \item{\code{beta==1}:}{real numbers}
#'               \item{\code{beta==2}:}{complex numbers}
#'               \item{\code{beta==4}:}{quaternion numbers}
#'             }
#' @param n size of matrix
#' 
#' @return The output of \code{pat_bnd} is a list of five elements:
#'         \itemize{
#'           \item \code{d}: ambient dimension of corresponding Euclidean space,
#'               given in terms of \code{beta,n} by \code{n+beta*choose(n,2)}
#'           \item \code{k_low}: function returning the lower bound for index \code{k}, given \code{r}
#'           \item \code{k_upp}: function returning the upper bound for index \code{k}, given \code{r}
#'           \item \code{r_low}: function returning the lower bound for index \code{r}, given \code{k}
#'           \item \code{r_upp}: function returning the upper bound for index \code{r}, given \code{k}
#'         }
#' 
#' @section See also:
#' \code{\link[symconivol]{leigh}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' # considering the case of 3x3 complex unitary matrices
#' beta <- 2
#' n <- 3
#' P <- pat_bnd(beta,n)
#' 
#' # ambient dimension
#' P$d
#' 
#' # bounds for index k, given r
#' P$k_low(0:n)
#' P$k_upp(0:n)
#' 
#' # bounds for index r, given k
#' P$r_low(0:P$d)
#' P$r_upp(0:P$d)
#' 
#' @export
#'
pat_bnd <- function(beta,n) {
    return( list( d = n+beta*choose(n,2) ,
                  k_low = function(r) return(r+beta*choose(r,2)) ,
                  k_upp = function(r) return(r+beta*choose(r,2)+beta*r*(n-r)) ,
                  r_low = function(k) return(
                      floor( (beta/2-1+sqrt(beta^2/4+beta*(2*k-1)+1))/beta ) ) ,
                  r_upp = function(k) return(
                      ceiling( (2*beta*n+2-beta-sqrt( (beta-2*beta*n-2)^2-8*beta*k ))/(2*beta) ) )
    ) )
}


#' Leigh's curve
#' 
#' \code{leigh} produces a table and lookup functions for Leigh's curve
#' (see accompanying vignette for definition).
#' 
#' @param N number of intermediate points; size of resulting table
#' 
#' @return The output of \code{leigh} is a list of three elements:
#'         \itemize{
#'           \item \code{tab}: a table of function values with \code{N} elements
#'           \item \code{lkup_rho}: a corresponding lookup function for rho in terms of kappa
#'           \item \code{lkup_kappa}: a corresponding lookup function for kappa in terms of rho
#'         }
#' 
#' @section See also:
#' \code{\link[symconivol]{pat_bnd}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' L <- leigh()
#' ggplot(L$tab, aes(x=rho,y=kappa)) + geom_line()
#' x <- (0:10)/10
#' matrix(c(x, L$lkup_rho(x) ),11,2)
#' matrix(c(x, L$lkup_kappa(x) ),11,2)
#' 
#' @export
#'
leigh <- function(N=1e3) {
    c_a <- function(a) {
        integrand <- function(x) return(sqrt((1-x)/x*(x^2+x+(a-1)/a^2)))
        return(2/pi/(1-(a-1)/a^2)*integrate(integrand, lower=0, upper=1)$value)
    }
    L_a <- function(a) return(a*sqrt(2/(a^2-a+1)))
    a <- seq(1,2,length.out=N)
    c <- sapply(a, c_a)
    L2 <- sapply(a, L_a)
    bnd_a <- function(a) {
        L <- L_a(a)
        return(list(neg_low=-L/a, neg_upp=-L*(1-1/a), pos_low=0, pos_upp=L))
    }
    rho_a <- function(a,moment=0) {
        L <- L_a(a)
        return( function(x) x^moment/pi*sqrt( (L-x)/x * (x+L/a) * (x+(1-1/a)*L) ) )
    }
    int_neg <- function(a) {
        bnd <- bnd_a(a)
        return(integrate(rho_a(a,2),lower=bnd$neg_low,upper=bnd$neg_upp)$value)
    }
    curve_half <- 2*sapply(a, int_neg)
    tab <- tibble( rho=c(rev(1-c),c), kappa=c(rev(curve_half),1-curve_half) )
    lkup_rho <- function(kappa) {
        return( tab$rho[ sapply( kappa, function(x) return(which.min((tab$kappa-x)^2)) ) ] )
    }
    lkup_kappa <- function(rho) {
        return( tab$kappa[ sapply( rho, function(x) return(which.min((tab$rho-x)^2)) ) ] )
    }
    return( list( table=tab, lkup_rho=lkup_rho, lkup_kappa=lkup_kappa ) )
}


#' Algebraic degree of semidefinite programming
#' 
#' \code{alg_deg} produces a table for the algebraic degree of semidefinite
#' programmin (see the accompanying vignette for reference to literature).
#' 
#' @param n size of matrix
#' 
#' @return The output of \code{alg_deg} is a table of the algebraic degrees
#'         for the given value of \code{n}. To avoid loss of loss of precision,
#'         the algebraic degrees are given as strings. See the below example for
#'         a simple conversion to numeric, including the problems associated with
#'         that step.
#' 
#' @section See also:
#' \code{\link[symconivol]{curv_meas_exact}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' AD_str <- alg_deg(6)
#' # convert to integer matrix
#' AD <- matrix( as.numeric(AD_str), dim(AD_str) )
#' colnames(AD) <- colnames(AD_str)
#' rownames(AD) <- rownames(AD_str)
#' # compare both matrices
#' print(AD_str)
#' print(AD)
#' 
#' # doing the same for larger n
#' AD_str <- alg_deg(14)
#' AD <- matrix( as.numeric(AD_str), dim(AD_str) )
#' # the conversion introduces rounding errors
#' i <- which.max(AD_str)
#' print(AD_str[i])
#' print(AD[i],digits=20)
#' 
#' @export
#' 
alg_deg <- function(n) {
    filename <- str_c("data/alg_deg/n=",n,".txt")
    if ( !file.exists(filename) )
        stop(str_c("Data for n=",n," not available."))

    # read file into character vector (first line is just describing text)
    tmp <- readLines(filename)[-1]

    d <- n + choose(n,2)
    T <- matrix("0",d-1,n-1)
    r <- 1
    m <- choose(n,2)
    for (i in 1:length(tmp)) {
        T[m,r] <- tmp[i]
        T[d-m,n-r] <- tmp[i]
        if (m==(choose(n+1,2)-choose(r+1,2))) {
            r <- r+1
            m <- choose(n-r+1,2)
        } else {
            m <- m+1
        }
    }
    colnames(T) <- str_c("r=",1:(n-1))
    rownames(T) <- str_c("m=",1:(d-1))
    return(T)
}






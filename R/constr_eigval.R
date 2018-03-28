#' Stan model creation for sampling index constrained eigenvalues from the Gaussian orthogonal/unitary/symplectic ensembles
#' 
#' \code{constr_eigval} generates inputs for Stan (model string or external file)
#' for sampling eigenvalues from the Gaussian orthogonal/unitary/symplectic ensembles
#' with prescriped range for the index (number of positive eigenvalues).
#' 
#' @param beta Dyson index specifying the underlying (skew-) field:
#'             \describe{
#'               \item{\code{beta==1}:}{real numbers}
#'               \item{\code{beta==2}:}{complex numbers}
#'               \item{\code{beta==4}:}{quaternion numbers}
#'             }
#' @param n size of matrix.
#' @param ind_low lower bound for index
#' @param ind_upp upper bound for index
#' @param filename filename for output
#' @param overwrite logical; determines whether the output should overwrite an existing file
#' 
#' @return The output of \code{constr_eigval} is a list containing the following elements:
#' \itemize{
#'   \item \code{model}: a string that forms the description of the Stan model,
#'   \item \code{data}: a data list containing the prepared data to be used
#'                    for defining a Stan model object (\code{beta}, number of positive, free,
#'                    negative eigenvalues).
#' }
#' If \code{filename!=NA} then the model string will also be written to the file with
#' the specified name.
#' 
#' @section See also:
#' \code{\link[symconivol]{constr_eigval_to_bcbsq}}, 
#' \code{\link[symconivol]{prepare_em_cm}}, 
#' \code{\link[symconivol]{estim_em_cm}}
#' 
#' Package: \code{\link[symconivol]{symconivol}}
#' 
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(rstan)
#' 
#' filename <- "tmp.stan"
#' M <- constr_eigval( beta=2, n=12, ind_low=8, ind_upp=9, filename=filename )
#' stan_samp <- stan( file = filename, data = M$data,
#'                    chains = 1, warmup = 1e3, iter = 1e5, cores = 2, refresh = 1e4 )
#' file.remove(filename)
#' 
#' tib_ep <- rstan::extract(stan_samp)$ep %>% as_tibble() %>% gather() %>% add_column(type="0")
#' tib_ef <- rstan::extract(stan_samp)$ef %>% as_tibble() %>% gather() %>% add_column(type="1")
#' tib_en <- rstan::extract(stan_samp)$en %>% as_tibble() %>% gather() %>% add_column(type="2")
#' tib <- bind_rows(tib_ep, tib_ef, tib_en)
#' 
#' ggplot() +
#'     geom_density(data=tib, aes(x=value, y=..count.., group=type, color=type, fill=type), alpha=0.5, bw=0.03) + theme_bw() +
#'     theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")
#' }
#' 
#' @export
#'
constr_eigval <- function(beta, n, ind_low, ind_upp, filename=NA, overwrite=FALSE) {
    if (!(n>0 & ind_low>=0 & ind_low<=ind_upp & ind_upp<=n)) stop("\n Invalid parameters (n>0, 0<=ind_low<=ind_upp<=n).")
    pos  <- ind_low
    free <- ind_upp-ind_low
    neg  <- n-ind_upp
    modelinfo <- paste0("
// Model for sampling eigenvalues from the Gaussian orthogonal/unitary/symplectic ensemble
// constrained on the index, the number of positive eigenvalues.
// Model is created with the R-method 'constr_eigval' from the 'symconivol' package (v.",
packageVersion("symconivol"),").
// See 'https://github.com/damelunx/symconivol' for more information.
\n")
    
    # data definition
    # 
    model_string <- paste0(modelinfo, "data {
    int <lower=1> beta ;                                        // Dyson index")
    if (pos) model_string <- paste0(model_string,"
    int <lower=1> np ;                                          // number of positive eigenvalues")
    if (free) model_string <- paste0(model_string,"
    int <lower=1> nf ;                                          // number of free eigenvalues")
    if (neg) model_string <- paste0(model_string,"
    int <lower=1> nn ;                                          // number of negative eigenvalues")
    model_string <- paste0(model_string,"\n}")
    
    # parameter definition
    # 
    model_string <- paste0(model_string,"\n\nparameters{")
    if (pos)                     model_string <- paste0(model_string,"
    positive_ordered[np] ep ;                                   // positive eigenvalues")
    if (neg)                     model_string <- paste0(model_string,"
    positive_ordered[nn] en_abs ;                               // absolute of negative eigenvalues")
    if (free & pos & neg)        model_string <- paste0(model_string,"
    simplex[nf+1] ef_gaps ;                                     // gaps between free eigenvalues")
    else if (free & pos & !neg)  model_string <- paste0(model_string,"
    positive_ordered[nf] ef_ep ;                                // smallest positive minus free eigenvalues")
    else if (free & !pos & neg)  model_string <- paste0(model_string,"
    positive_ordered[nf] ef_en ;                                // free minus smallest negative eigenvalues")
    else if (free & !neg & !pos) model_string <- paste0(model_string,"
    ordered[nf] ef ;                                            // free eigenvalues")
    model_string <- paste0(model_string,"\n}")
    
    # transformed parameters for free eigenvalues
    # 
    if (free & (neg | pos)) model_string <- paste0(model_string,"\n\ntransformed parameters {")
    if (free & pos & neg)        model_string <- paste0(model_string,"
    ordered[nf] ef ;                                            // free eigenvalues
    ef = -en_abs[1] + head(cumulative_sum(ef_gaps),nf) * (ep[1]+en_abs[1]) ;")
    else if (free & pos & !neg)  model_string <- paste0(model_string,"
    ordered[nf] ef_min ;                                        // minus free eigenvalues
    ef_min = ef_ep - ep[1] ;")
    else if (free & !pos & neg)  model_string <- paste0(model_string,"
    ordered[nf] ef ;                                            // free eigenvalues
    ef = ef_en - en_abs[1] ;")
    if (free & (neg | pos)) model_string <- paste0(model_string,"\n}")

    # model definition
    # 
    model_string <- paste0(model_string,"\n\nmodel {
    // Gaussian weight:
    target += -(0")
    if (pos)  model_string <- paste0(model_string,"+dot_self(ep)")
    if (free) {
        if (free & pos & !neg)
            model_string <- paste0(model_string,"+dot_self(ef_min)")
        else
            model_string <- paste0(model_string,"+dot_self(ef)")
    }
    if (neg)  model_string <- paste0(model_string,"+dot_self(en_abs)")
    model_string <- paste0(model_string,")/2 ;")
    
    if (pos) model_string <- paste0(model_string,"\n
    if (np>=2) {                                                // Vandermonde positive
        for (k in 2:np) {
            target += beta*sum(log(ep[k]-ep[1:(k-1)])) ;
        }
    }")
    if (free) {
        if (free & pos & !neg)
            model_string <- paste0(model_string,"\n
    if (nf>=2) {                                                // Vandermonde free
        for (k in 2:nf) {
            target += beta*sum(log(ef_min[k]-ef_min[1:(k-1)])) ;
        }
    }")
        else
            model_string <- paste0(model_string,"\n
    if (nf>=2) {                                                // Vandermonde free
        for (k in 2:nf) {
            target += beta*sum(log(ef[k]-ef[1:(k-1)])) ;
        }
    }")
    }
    if (neg) model_string <- paste0(model_string,"\n
    if (nn>=2) {                                                // Vandermonde negative
        for (k in 2:nn) {
            target += beta*sum(log(en_abs[k]-en_abs[1:(k-1)])) ;
        }
    }")
    
    if (neg & pos) model_string <- paste0(model_string,"\n
    for (k in 1:nn) {                                           // interaction negative-positive
        target += beta*sum(log(ep+en_abs[k])) ;
    }")
    if (neg & free) model_string <- paste0(model_string,"\n
    for (k in 1:nf) {                                           // interaction negative-free
        target += beta*sum(log(en_abs+ef[k])) ;
    }")
    if (pos & free & neg) model_string <- paste0(model_string,"\n
    for (k in 1:nf) {                                           // interaction positive-free
        target += beta*sum(log(ep-ef[k])) ;
    }")
    if (pos & free & !neg) model_string <- paste0(model_string,"\n
    for (k in 1:nf) {                                           // interaction positive-free
        target += beta*sum(log(ep+ef_min[k])) ;
    }")
    model_string <- paste0(model_string,"\n}")
    
    # generated quantities for negative eigenvalues
    # 
    if (neg) model_string <- paste0(model_string,"\n\ngenerated quantities {
    vector[nn] en ;                                            // negative eigenvalues
    en = -en_abs ;\n}")

    # generated quantities for free eigenvalues (if no negative eigenvalues)
    # 
    if (pos & free & !neg) model_string <- paste0(model_string,"\n\ngenerated quantities {
    vector[nf] ef ;                                            // free eigenvalues
    ef = -ef_min ;\n}")
    
    out <- list()
    out$data <- list( beta=beta, nn=neg, nf=free, np=pos ) 
    out$model <- model_string
    if ( !is.na(filename) && file.exists(filename) && overwrite==FALSE )
        stop("\n File with given filename exists and overwrite==FALSE.")
    else if ( !is.na(filename) ) {
        if (file.exists(filename))
            file.remove(filename)
        file_conn<-file(filename)
        writeLines(model_string, file_conn)
        close(file_conn)
    }
    return(out)
}


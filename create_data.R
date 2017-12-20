
library(conivol)
library(tidyverse)
library(rstan)

# SAMPLE EIGENVALUES -- DIVERSE FOR VIGNETTE
# 
# warmup <- 1e3
# num_samp <- 1e5+warmup
# format <- tribble( ~N, ~F, ~P,
#                     5, 10, 25,
#                     0, 15, 25,
#                    15,  0, 25,
#                    35,  5,  0,
#                    20,  0,  0,
#                     0, 40,  0,
#                     0,  0, 40
#                    )
# filename <- "tmp.stan"
# for (i in 1:dim(format)[1]) {
#     nn <- format[[i,1]]
#     nf <- format[[i,2]]
#     np <- format[[i,3]]
#     invisible( tmp <- constr_eigval( pos=(np>0), free=(nf>0), neg=(nn>0),
#                                      filename=filename, overwrite=TRUE ) )
#     for (beta in c(1,2,4)) {
#         data = list(beta=beta)
#         if (np>0) data$np <- np
#         if (nf>0) data$nf <- nf
#         if (nn>0) data$nn <- nn
#         stan_samp <- stan( file = filename, data = data, chains = 1, warmup = warmup,
#                            iter = num_samp, cores = 2, refresh = 1e4 )
#         samp <- list()
#         if (np>0) samp$ep <- rstan::extract(stan_samp)$ep
#         if (nf>0) samp$ef <- rstan::extract(stan_samp)$ef
#         if (nn>0) samp$en <- rstan::extract(stan_samp)$en
#         save(samp, file=str_c("data/beta_nn_nf_np=",beta,"_",nn,"_",nf,"_",np,"_wmup=",warmup,".RData"))
#     }
# }


# SAMPLE EIGENVALUES -- UNCONSTRAINED
# 
# warmup <- 1e3
# nmin <- 3
# nmax <- 13
# num_samp <- 1e6+warmup
# filename <- "data/unconstrained.stan"
# for (beta in c(1,2,4)) {
#     for (n in nmin:nmax) {
#         stan_inp <- list( n=n, beta=beta )
#         stan_samp <- stan( file=filename, data=stan_inp,
#                            chains=1, warmup=warmup, iter=num_samp, cores=2, refresh=1e4 )
#         samp <- rstan::extract(stan_samp)$e
#         save(samp, file=str_c("data/beta=",beta,"/unc_beta=",beta,"_n=",n,"_wmup=",warmup,".RData"))
#     }
# }


# SAMPLE EIGENVALUES -- CONSTRAINED
# 
# warmup <- 1e3
# nmin <- 14
# nmax <- 20
# num_samp <- 1e5+warmup
# filename <- "data/constrained.stan"
# for (beta in c(1,2,4)) {
# # beta <- 1
#     for (n in nmin:nmax) {
#         for (np in ceiling(n/2):(n-1)) {
#             stan_inp <- list( np=np, nn=n-np, beta=beta )
#             stan_samp <- stan( file=filename, data=stan_inp,
#                                chains=1, warmup=warmup, iter=num_samp, cores=2, refresh=1e4 )
#             samp <- list( pos=rstan::extract(stan_samp)$ep, neg=rstan::extract(stan_samp)$en )
#             save(samp, file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",np,"_wmup=",warmup,".RData"))
#         }
#     }
# }


# SAMPLE EIGENVALUES -- CONSTRAINED - ALL POSITIVE
# 
# warmup <- 1e3
# nmin <- 14
# nmax <- 20
# num_samp <- 1e5+warmup
# filename <- "data/constrained_allpos.stan"
# for (beta in c(1,2,4)) {
# # beta <- 1
#     for (n in nmin:nmax) {
#         np <- n
#         stan_inp <- list( np=np, beta=beta )
#         stan_samp <- stan( file=filename, data=stan_inp,
#                            chains=1, warmup=warmup, iter=num_samp, cores=2, refresh=1e4 )
#         samp <- rstan::extract(stan_samp)$ep
#         save(samp, file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",np,"_wmup=",warmup,".RData"))
#     }
# }


# RECONSTRUCT INTRINSIC VOLUMES (100 STEPS OF EM)
# 
# warmup <- 1e3
# nmin <- 3
# nmax <- 13
# for (beta in c(1,2,4)) {
#     for (n in nmin:nmax) {
#         d <- n + beta*choose(n,2)
#         Rec_iv <- list( beta=beta, n=n )
#         load(file=str_c("data/beta=",beta,"/unc_beta=",beta,"_n=",n,"_wmup=",warmup,".RData"))
#         bcb_unconstr <- samp %>%
#             apply( 1 ,
#                    function(row) return(c(sum(row[row>0]^2), sum(row[row<0]^2))) ) %>%
#             t() %>% as_tibble()
#         bcb_unconstr_mat <- as.matrix(bcb_unconstr)
#         est <- estim_em( d, bcb_unconstr_mat, N=100, init_mode=0, no_of_lcc_projections=1  )
#         Rec_iv$EM_its <- est$iterates
#         save(Rec_iv, file=str_c("data/beta=",beta,"/IV_beta=",beta,"_n=",n,".RData"))
#     }
# }


# RECONSTRUCT (RANK) NORMALIZED CURVATURE MEASURES (100 STEPS OF EM)
# 
# warmup <- 1e3
# nmin <- 14
# nmax <- 20
# for (beta in c(1,2,4)) {
#     for (n in nmin:nmax) {
#         d <- n + beta*choose(n,2)
#         Rec_cm_rknorm <- list( beta=beta, n=n )
#         EM_its <- vector(mode="list", length=(n-1-ceiling(n/2)+1) )
#         names(EM_its) <- sapply(as.character(ceiling(n/2):(n-1)),
#                                 function(x) return(str_c("np=",x)))
#         for (np in ceiling(n/2):(n-1)) {
#             load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",np,"_wmup=",warmup,".RData"))
#             bcb_constr <- rowSums(samp$pos) %>% as_tibble() %>%
#                 add_column(rowSums(samp$neg)) %>% `colnames<-`(c("V1","V2"))
#             bcb_constr_mat <- as.matrix(bcb_constr)
#             est <- estim_em( d, bcb_constr_mat, N=100, init_mode=0, no_of_lcc_projections=1  )
#             EM_its[[np-ceiling(n/2)+1]] <- est$iterates
#         }
#         Rec_cm_rknorm$EM_its <- EM_its
#         save(Rec_cm_rknorm, file=str_c("data/beta=",beta,"/CM_rknorm_beta=",beta,"_n=",n,".RData"))
#     }
# }
















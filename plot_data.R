
library(conivol)
library(tidyverse)
library(rstan)
library(Rmisc)

# EIGENVALUES -- DIVERSE FOR VIGNETTE
# 
# warmup <- 1e3
# format <- tribble( ~N, ~F, ~P,
#                     5, 10, 25,
#                     0, 15, 25,
#                    15,  0, 25,
#                    35,  5,  0,
#                    20,  0,  0,
#                     0, 40,  0,
#                     0,  0, 40 )
# nbins <- 500
# for (i in 1:dim(format)[1]) {
#     nn <- format[[i,1]]
#     nf <- format[[i,2]]
#     np <- format[[i,3]]
#     plotsPB <- list()
#     j <- 0
#     for (beta in c(1,2,4)) {
#         j <- j+1
#         load(file=str_c("data/beta_nn_nf_np=",beta,"_",nn,"_",nf,"_",np,"_wmup=",warmup,".RData"))
#         plotsPB[[j]] <- ggplot()
#         if (np>0) {
#             tmp <- samp$ep %>% as_tibble() %>% gather(factor_key=TRUE)
#             plotsPB[[j]] <- plotsPB[[j]] +
#                 geom_histogram(data=tmp, aes(x=value, y=..density..), bins=nbins, fill="green", alpha=0.8)
#         }
#         if (nf>0) {
#             tmp <- samp$ef %>% as_tibble() %>% gather(factor_key=TRUE)
#             plotsPB[[j]] <- plotsPB[[j]] +
#                 geom_histogram(data=tmp, aes(x=value, y=..density..), bins=nbins, fill="blue", alpha=0.5)
#         }
#         if (nn>0) {
#             tmp <- samp$en %>% as_tibble() %>% gather(factor_key=TRUE)
#             plotsPB[[j]] <- plotsPB[[j]] +
#                 geom_histogram(data=tmp, aes(x=value, y=..density..), bins=nbins, fill="red", alpha=0.8)
#         }
#         plotsPB[[j]] <- plotsPB[[j]] +
#             scale_y_continuous(labels = scales::percent) + theme_bw() + xlab(str_c("beta=",beta)) +
#             theme(axis.title.y=element_blank(), #axis.title.y=element_blank(),
#                   axis.text.y=element_blank(), axis.ticks.y=element_blank(),
#                   panel.border = element_blank(), panel.grid.major = element_blank(),
#                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#     }
#     plotFileName <-str_c("vignettes/curv-meas-from-constr-evals_figures/",
#                          "nn_nf_np=",nn,"_",nf,"_",np,".png")
#     ggsave(plotFileName, plot=multiplot(plotlist = plotsPB, cols = 3), device="png")
# }
# multiplot(plotlist = plotsPB, cols = 3)

# EIGENVALUES -- UNCONSTRAINED
# 
# warmup <- 1e3
# beta <- 1
# n <- 10
# 
# load(file=str_c("data/beta=",beta,"/unc_beta=",beta,"_n=",n,"_wmup=",warmup,".RData"))
# 
# S <- samp %>% as_tibble() %>% gather(factor_key=TRUE)
# ggplot(S, aes(value)) +
#     geom_histogram(aes(y=..density..), bins=300) +
#     scale_y_continuous(labels = scales::percent) + theme_bw() +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#           axis.text.y=element_blank(), axis.ticks.y=element_blank())


# EIGENVALUES -- CONSTRAINED
# 
# warmup <- 1e3
# beta <- 2
# n <- 10
# np <- 3
# nbins <- 800
# 
# if (np==n) {
#     load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",n,"_wmup=",warmup,".RData"))
#     Sp <- samp %>% as_tibble() %>% gather(factor_key=TRUE)
#     ggplot() +
#         geom_histogram(data=Sp, aes(x=value, y=..density..), bins=nbins) +
#         scale_y_continuous(labels = scales::percent) + theme_bw() +
#         theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#               axis.text.y=element_blank(), axis.ticks.y=element_blank())
# } else if (np==0) {
#     load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",n,"_wmup=",warmup,".RData"))
#     Sn <- samp %>% as_tibble() %>% gather(factor_key=TRUE)
#     ggplot() +
#         geom_histogram(data=Sn, aes(x=-value, y=..density..), bins=nbins) +
#         scale_y_continuous(labels = scales::percent) + theme_bw() +
#         theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#               axis.text.y=element_blank(), axis.ticks.y=element_blank())
# } else if (np>=ceiling(n/2) && np<n) {
#     load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",np,"_wmup=",warmup,".RData"))
#     Sp <- samp$pos %>% as_tibble() %>% gather(factor_key=TRUE)
#     Sn <- samp$neg %>% as_tibble() %>% gather(factor_key=TRUE)
#     ggplot() +
#         geom_histogram(data=Sp, aes(x=value, y=..density..), bins=nbins) +
#         geom_histogram(data=Sn, aes(x=-value, y=..density..), bins=nbins) +
#         scale_y_continuous(labels = scales::percent) + theme_bw() +
#         theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#               axis.text.y=element_blank(), axis.ticks.y=element_blank())
# } else if (np>0) {
#     load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",n-np,"_wmup=",warmup,".RData"))
#     Sn <- samp$pos %>% as_tibble() %>% gather(factor_key=TRUE)
#     Sp <- samp$neg %>% as_tibble() %>% gather(factor_key=TRUE)
#     ggplot() +
#         geom_histogram(data=Sp, aes(x=value, y=..density..), bins=nbins) +
#         geom_histogram(data=Sn, aes(x=-value, y=..density..), bins=nbins) +
#         scale_y_continuous(labels = scales::percent) + theme_bw() +
#         theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#               axis.text.y=element_blank(), axis.ticks.y=element_blank())
# } else print("Invalid value for np.")


# BIVARIATE CHI-BAR SQUARED DATA -- UNCONSTRAINED
# 
# warmup <- 1e3
# beta <- 1
# n <- 10
# 
# load(file=str_c("data/beta=",beta,"/unc_beta=",beta,"_n=",n,"_wmup=",warmup,".RData"))
# 
# bcb <- samp %>%
#     apply( 1 , function(row) return(c(sum(row[row>0]^2), sum(row[row<0]^2))) ) %>%
#     t() %>% as_tibble()
# 
# ggplot( bcb, aes(V1,V2)) + geom_point(alpha=.02) +
#     theme_bw() +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank())


# BIVARIATE CHI-BAR SQUARED DATA -- CONSTRAINED
# 
# warmup <- 1e3
# beta <- 2
# n <- 10
# np <- 4
# nbins <- 200
# N_plot <- 1e2
# 
# if (np %in% c(0,n)) {
#     load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",n,"_wmup=",warmup,".RData"))
#     csq <- rowSums(samp^2)
#     d <- n+beta*choose(n,2)
#     x <- seq(d-5*sqrt(d),d+7*sqrt(d),length.out=N_plot)
#     tib_chisq <- tibble(x=x, y=dchisq(x,d))
#     ggplot() +
#         geom_histogram(data=tibble(x=csq), aes(x=x, y=..density..), bins=nbins) +
#         geom_line(data=tib_chisq, aes(x=x,y=y), color="red",size=2,alpha=0.5) +
#         theme_bw() +
#         theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#               axis.text.y=element_blank(), axis.ticks.y=element_blank())
# } else if (np>=ceiling(n/2) && np<n) {
#     load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",np,"_wmup=",warmup,".RData"))
#     bcb <- rowSums(samp$pos) %>% as_tibble() %>%
#         add_column(rowSums(samp$neg)) %>% `colnames<-`(c("V1","V2"))
#     ggplot( bcb, aes(V1,V2)) + geom_point(alpha=.02) +
#         theme_bw() +
#         theme(axis.title.x=element_blank(),axis.title.y=element_blank())
# } else if (np>0) {
#     load(file=str_c("data/beta=",beta,"/con_beta=",beta,"_n=",n,"_np=",n-np,"_wmup=",warmup,".RData"))
#     bcb <- rowSums(samp$neg) %>% as_tibble() %>%
#         add_column(rowSums(samp$pos)) %>% `colnames<-`(c("V1","V2"))
#     ggplot( bcb, aes(V1,V2)) + geom_point(alpha=.02) +
#         theme_bw() +
#         theme(axis.title.x=element_blank(),axis.title.y=element_blank())
# } else print("Invalid value for np.")


# INTRINSIC VOLUMES - EM ITERATES
# 
# beta <- 1
# n <- 12
# d <- n + beta*choose(n,2)
# load(file=str_c("data/beta=",beta,"/IV_beta=",beta,"_n=",n,".RData"))
# tib_plot <- as_tibble( t(Rec_iv$EM_its[1+12*(0:8), ]) ) %>%
#     `colnames<-`(paste0("s_",12*(0:8))) %>%
#     add_column(k=0:d,.before=1) %>% gather(step,value,2:10)
# tib_plot$step <- factor(tib_plot$step, levels = paste0("s_",12*(0:8)))
# ggplot(tib_plot,aes(x=k,y=value,color=step)) +
#     geom_line() + theme_bw() +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank())


# CURVATURE MEASURES - RANK NORMALIZED - VIOLIN PLOTS
# 
# beta <- 1
# n <- 6
# d <- n + beta*choose(n,2)
# load(file=str_c("data/beta=",beta,"/CM_rknorm_beta=",beta,"_n=",n,".RData"))
# num_cm <- length(Rec_cm_rknorm$EM_its)
# N <- 1e5
# tib_tmp <- tibble()
# for (i in 1:num_cm) {
#     w <- Rec_cm_rknorm$EM_its[[i]][101,]
#     tmp <- numeric()
#     for (k in 0:d) {
#         tmp <- c(tmp,rep(k, round(w[k+1]*N)))
#     }
#     tib_tmp2 <- tibble(r=ceiling(n/2)+i-1, k=tmp)
#     tib_tmp <- bind_rows(tib_tmp,tib_tmp2)
#     if (ceiling(n/2)+i-1 != floor(n/2)-i+1) {
#         tib_tmp2 <- tibble(r=floor(n/2)-i+1, k=d-tmp)
#         tib_tmp <- bind_rows(tib_tmp,tib_tmp2)
#     }
# }
# ggplot(tib_tmp, aes(factor(r),k)) + geom_violin(bw=0.6) + coord_flip() +
#     theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())


# CURVATURE MEASURES - RANK NORMALIZED
# 
# beta <- 1
# n <- 13
# d <- n + beta*choose(n,2)
# load(file=str_c("data/beta=",beta,"/CM_rknorm_beta=",beta,"_n=",n,".RData"))
# tib_shape <- tibble()
# for (r in 0:n) {
#     pat_low <- r+beta*choose(r,2)
#     pat_upp <- r+beta*(choose(r,2)+r*(n-r))
#     tib_r <- tibble( r=r, k=0:d, pat = 0:d>=pat_low & 0:d<=pat_upp )
#     if (r==0)
#         tib_r <- add_column( tib_r, cm=c(1,rep(0,d)) )
#     else if (r==n)
#         tib_r <- add_column( tib_r, cm=c(rep(0,d),1) )
#     else if (r>=ceiling(n/2))
#         tib_r <- add_column( tib_r, cm=Rec_cm_rknorm$EM_its[[r-ceiling(n/2)+1]][101,] )
#     else
#         tib_r <- add_column( tib_r, cm=rev(Rec_cm_rknorm$EM_its[[floor(n/2)-r+1]][101,]) )
#     tib_r$cm <- tib_r$cm/max(tib_r$cm)
#     tib_shape <- bind_rows(tib_shape, tib_r)
# }
# y_pat <- seq(0,1,length.out=1e2)
# pat_bd_low <- tibble(patX=d*y_pat^2, patY=n*y_pat)
# pat_bd_upp <- tibble(patX=d*(1-y_pat^2), patY=n*(1-y_pat))
# tib_diag <- tibble(x=c(0,d),y=c(0,n))
# ggplot() +
#     geom_line(data=pat_bd_low, aes(patX, patY), linetype=2) +
#     geom_line(data=pat_bd_upp, aes(patX, patY), linetype=2) +
#     geom_point(data=tib_shape, aes(k,r,color=pat,size=cm,alpha=cm)) +
#     scale_color_manual(values=c("gray80", "black")) +
#     theme_bw() + theme(legend.position="none") +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#           panel.border = element_blank(), panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#     geom_hline(yintercept=0.5+0:(n-1), linetype=3) +
#     geom_line(data=tib_diag, aes(x=x,y=y), color="red",size=2,alpha=0.5)


# CURVATURE MEASURES - RENORMALIZED
# 
# beta <- 1
# n <- 11
# d <- n + beta*choose(n,2)
# load(file=str_c("data/beta=",beta,"/IV_beta=",beta,"_n=",n,".RData"))
# iv <- Rec_iv$EM_its[101, ]
# cm_rk <- matrix(0,n+1,d+1)
# load(file=str_c("data/beta=",beta,"/CM_rknorm_beta=",beta,"_n=",n,".RData"))
# for (r in 0:n) {
#     if (r==0)
#         cm_rk[r+1,] <- c(1,rep(0,d))
#     else if (r==n)
#         cm_rk[r+1,] <- c(rep(0,d),1)
#     else if (r>=ceiling(n/2))
#         cm_rk[r+1,] <- Rec_cm_rknorm$EM_its[[r-ceiling(n/2)+1]][101,]
#     else
#         cm_rk[r+1,] <- rev(Rec_cm_rknorm$EM_its[[floor(n/2)-r+1]][101,])
# }
# rk <- lsfit(t(cm_rk),iv,intercept=FALSE)$coefficients
# names(rk) <- str_c("r=",0:n)
# cm <- rk * cm_rk
# tib_shape <- tibble()
# for (r in 0:n) {
#     pat_low <- r+beta*choose(r,2)
#     pat_upp <- r+beta*(choose(r,2)+r*(n-r))
#     tib_r <- tibble( r=r, k=0:d, pat = 0:d>=pat_low & 0:d<=pat_upp ) %>%
#         add_column( cm=cm[r+1,] )
#     tib_shape <- bind_rows(tib_shape, tib_r)
# }
# y_pat <- seq(0,1,length.out=1e2)
# pat_bd_low <- tibble(patX=d*y_pat^2, patY=n*y_pat)
# pat_bd_upp <- tibble(patX=d*(1-y_pat^2), patY=n*(1-y_pat))
# tib_diag <- tibble(x=c(0,d),y=c(0,n))
# ggplot() +
#     geom_line(data=pat_bd_low, aes(patX, patY), linetype=2) +
#     geom_line(data=pat_bd_upp, aes(patX, patY), linetype=2) +
#     geom_point(data=tib_shape, aes(k,r,color=pat,size=cm,alpha=cm)) +
#     scale_color_manual(values=c("gray80", "black")) +
#     theme_bw() + theme(legend.position="none") +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#           panel.border = element_blank(), panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#     geom_line(data=tib_diag, aes(x=x,y=y), color="red",size=2,alpha=0.5)


# CURVATURE MEASURES - INTRINSIC VOLUMES NORMALIZED
# # 
# beta <- 1
# n <- 11
# d <- n + beta*choose(n,2)
# load(file=str_c("data/beta=",beta,"/IV_beta=",beta,"_n=",n,".RData"))
# iv <- Rec_iv$EM_its[101, ]
# cm_rk <- matrix(0,n+1,d+1)
# load(file=str_c("data/beta=",beta,"/CM_rknorm_beta=",beta,"_n=",n,".RData"))
# for (r in 0:n) {
#     if (r==0)
#         cm_rk[r+1,] <- c(1,rep(0,d))
#     else if (r==n)
#         cm_rk[r+1,] <- c(rep(0,d),1)
#     else if (r>=ceiling(n/2))
#         cm_rk[r+1,] <- Rec_cm_rknorm$EM_its[[r-ceiling(n/2)+1]][101,]
#     else
#         cm_rk[r+1,] <- rev(Rec_cm_rknorm$EM_its[[floor(n/2)-r+1]][101,])
# }
# rk <- lsfit(t(cm_rk),iv,intercept=FALSE)$coefficients
# names(rk) <- str_c("r=",0:n)
# cm <- rk * cm_rk
# tib_shape <- tibble()
# for (r in 0:n) {
#     pat_low <- r+beta*choose(r,2)
#     pat_upp <- r+beta*(choose(r,2)+r*(n-r))
#     tib_r <- tibble( r=r, k=0:d, pat = 0:d>=pat_low & 0:d<=pat_upp ) %>%
#         add_column( cm=cm[r+1,] )
#     tib_shape <- bind_rows(tib_shape, tib_r)
# }
# y_pat <- seq(0,1,length.out=1e2)
# pat_bd_low <- tibble(patX=d*y_pat^2, patY=n*y_pat)
# pat_bd_upp <- tibble(patX=d*(1-y_pat^2), patY=n*(1-y_pat))
# for (k in 0:d) {
#     I <- which(tib_shape$k==k)
#     tib_shape$cm[I] <- tib_shape$cm[I] / max(tib_shape$cm[I])
# }
# tib_diag <- tibble(x=c(0,d),y=c(0,n))
# ggplot() +
#     geom_line(data=pat_bd_low, aes(patX, patY), linetype=2) +
#     geom_line(data=pat_bd_upp, aes(patX, patY), linetype=2) +
#     geom_point(data=tib_shape, aes(k,r,color=pat,size=cm,alpha=cm)) +
#     scale_color_manual(values=c("gray80", "black")) +
#     theme_bw() + theme(legend.position="none") +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#           panel.border = element_blank(), panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#     geom_vline(xintercept=0.5+0:(d-1), linetype=3) +
#     geom_line(data=tib_diag, aes(x=x,y=y), color="red",size=2,alpha=0.5)












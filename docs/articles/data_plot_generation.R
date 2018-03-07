
# THIS SCRIPT CONTAINS ALL COMMANDS TO CREATE THE DATA AND PLOTS SHOWN IN THE VIGNETTE


# LOAD SOME LIBRARIES
# 
library(conivol)
library(tidyverse)
library(rstan)
library(Rmisc)
library(png)
library(latex2exp)

warmup <- 1e3
num_samp <- 1e5
C <- c(0.6, 0.8)
N <- c(10, 40)
format <- tribble( ~N, ~F, ~P,
                   5, 10, 25,
                   0, 15, 25,
                   16,  0, 24,
                   35,  5,  0,
                   20,  0,  0,
                   0,  0, 40,
                   0, 40,  0,
                   0, 10,  0,
                   (1-C[1])*N[1],  0,  C[1]*N[1],
                   (1-C[1])*N[2],  0,  C[1]*N[2],
                   (1-C[2])*N[1],  0,  C[2]*N[1],
                   (1-C[2])*N[2],  0,  C[2]*N[2]
)
hues = seq(15, 375, length = 3 + 1)
colors <- hcl(h = hues, l = 65, c = 100)[c(2,3,1)]
ymax <- c(0.18,0.16,0.16,0.24,0.6,0.5,0.083,0,rep(0.9,2),rep(1.3,2))

# CREATE AND SAVE PLOTS OF EMPIRICAL EIGENVALUE DISTRIBUTION
# # 
for (i in 1:7) {
    nn <- format[[i,1]]
    nf <- format[[i,2]]
    np <- format[[i,3]]
    n <- nn+nf+np
    plotsPB <- list()
    j <- 0
    for (beta in c(1,2,4)) {
        j <- j+1
        load(file=str_c("data/beta_nn_nf_np=",beta,"_",nn,"_",nf,"_",np,"_wmup=",warmup,".RData"))
        tmp <- tibble()
        if (np>0) tmp <- bind_rows(tmp,samp$ep %>% as_tibble() %>% gather(factor_key=TRUE) %>% add_column(type="0"))
        if (nf>0) tmp <- bind_rows(tmp,samp$ef %>% as_tibble() %>% gather(factor_key=TRUE) %>% add_column(type="1"))
        if (nn>0) tmp <- bind_rows(tmp,samp$en %>% as_tibble() %>% gather(factor_key=TRUE) %>% add_column(type="2"))
        this_colors <- colors[c(np>0,nf>0,nn>0)]
        plotsPB[[j]] <- ggplot() +
            geom_density(data=tmp, aes(x=value, y=..count.., group=type, color=type, fill=type), alpha=0.5, bw=0.03) +
            theme_bw() + xlab(TeX(sprintf("$\\beta = %d$", beta))) +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                  panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
            coord_cartesian(ylim = c(0, ymax[i]/sqrt(beta)*num_samp*n)) +
            theme(legend.position="none") +
            scale_fill_manual(values=this_colors) + scale_color_manual(values=this_colors)
    }
    plotFileName <-str_c("vignettes/curv-meas_figures/",
                         "nn_nf_np=",nn,"_",nf,"_",np,".png")
    ggsave(plotFileName, plot=multiplot(plotlist = plotsPB, cols = 3), device="png",
           width=178, height=80, units="mm", dpi=300)
}


# CREATE AND SAVE PLOTS OF CONVERGENCE TO SEMICIRCLE
# 
# t <- seq(-pi/2,pi/2,length=100)
# C <- tibble(x=sqrt(2)*sin(t),y=sqrt(2)/pi*cos(t))
# for (i in 7:8) {
#     n <- format[[i,2]]
#     plotsPB <- list()
#     j <- 0
#     for (beta in c(1,2,4)) {
#         j <- j+1
#         filename <- str_c("data/beta_nn_nf_np=",beta,"_0_",n,"_0_wmup=",warmup,".RData")
#         load(file=filename)
#         tmp <- (1/sqrt(beta*n)*samp$ef) %>% as_tibble() %>% gather(factor_key=TRUE) %>% add_column(type="1")
#         plotsPB[[j]] <- ggplot() +
#             geom_density(data=tmp, aes(x=value, group=type, color=type, fill=type), alpha=0.5, bw=0.03) +
#             theme_bw() + xlab(str_c("beta=",beta)) +
#             theme(axis.title.y=element_blank(),
#                   axis.text.y=element_blank(), axis.ticks.y=element_blank(),
#                   panel.border = element_blank(), panel.grid.major = element_blank(),
#                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#             coord_cartesian(ylim = c(0, 1.2*sqrt(2)/pi)) +
#             theme(legend.position="none") + scale_color_manual(values=colors[2]) +
#             scale_fill_manual(values=colors[2]) +
#             geom_line(data=C, aes(x=x,y=y), color="red",size=2,alpha=0.3)
#     }
#     plotFileName <-str_c("vignettes/curv-meas-from-constr-evals_figures/",
#                          "semic_n=",n,".png")
#     ggsave(plotFileName, plot=multiplot(plotlist = plotsPB, cols = 3), device="png",
#            width=178, height=80, units="mm", dpi=300)
# }


# CREATE AND SAVE PLOTS OF CONVERGENCE TO INDEX CONSTRAINED LIMIT DISTRIBUTION
# 
c_a <- function(a) {
    integrand <- function(x) return(sqrt((1-x)/x*(x^2+x+(a-1)/a^2)))
    return(2/pi/(1-(a-1)/a^2)*integrate(integrand, lower=0, upper=1)$value)
}
L_a <- function(a) return(a*sqrt(2/(a^2-a+1)))
bnd_a <- function(a) {
    L <- L_a(a)
    return(list(neg_low=-L/a, neg_upp=-L*(1-1/a), pos_low=0, pos_upp=L))
}
rho_a <- function(a,moment=0) {
    L <- L_a(a)
    return( function(x) x^moment/pi*sqrt( (L-x)/x * (x+L/a) * (x+(1-1/a)*L) ) )
}
a <- seq(1,2,length.out=1e3)
c <- sapply(a, c_a)

for (i in 9:12) {
    nn <- format[[i,1]]
    np <- format[[i,3]]
    n <- nn+np
    c0 <- np/n
    ind <- which.min((c-c0)^2)
    a0 <- a[ind]
    bnd <- bnd_a(a0)
    x_neg <- seq(bnd$neg_low,bnd$neg_upp,length.out=1e3)
    x_pos <- seq(bnd$pos_low,bnd$pos_upp,length.out=1e3)
    data_lim_neg <- tibble(x=x_neg, y=sapply(x_neg, rho_a(a0))*n*num_samp)
    data_lim_pos <- tibble(x=x_pos, y=sapply(x_pos, rho_a(a0))*n*num_samp)
    plotsPB <- list()
    j <- 0
    for (beta in c(1,2,4)) {
        j <- j+1
        load(file=str_c("data/beta_nn_nf_np=",beta,"_",nn,"_0_",np,"_wmup=",warmup,".RData"))
        data_pos <- (1/sqrt(beta*n)*samp$ep) %>% as_tibble() %>% gather()
        data_neg <- (1/sqrt(beta*n)*samp$en) %>% as_tibble() %>% gather()
        plotsPB[[j]] <- ggplot() +
            geom_density(data=data_pos, aes(x=value, y=..count..), color=colors[1], fill=colors[1], alpha=0.5, bw=0.01) +
            geom_density(data=data_neg, aes(x=value, y=..count..), color=colors[3], fill=colors[3], alpha=0.5, bw=0.01) +
            geom_line(data=data_lim_neg, aes(x=x,y=y),color=colors[2],size=1,alpha=0.7) +
            geom_line(data=data_lim_pos, aes(x=x,y=y),color=colors[2],size=1,alpha=0.7) +
            theme_bw() + xlab(str_c("beta=",beta)) +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                  panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
            coord_cartesian(ylim = c(0, ymax[i]*n*num_samp)) +
            theme(legend.position="none")
    }
    plotFileName <-str_c("vignettes/curv-meas-from-constr-evals_figures/",
                         "indlim_n_cPerc=",n,"_",round(c0*100),".png")
    ggsave(plotFileName, plot=multiplot(plotlist = plotsPB, cols = 3), device="png",
           width=178, height=80, units="mm", dpi=300)
}




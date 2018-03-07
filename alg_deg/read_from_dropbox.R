library(tidyverse)
# install.packages("rdrop2")
# install.packages("httpuv")
library(rdrop2)
library(Rmisc)
library(latex2exp)
library(symconivol)

# outputDir <- "algebr_deg"
# getFileFromDropbox <- function(filename) {
#     return(drop_read_csv(str_c(outputDir,"/",filename),header=FALSE,skip=1))
# }
# saveFileToDropbox <- function(filename,data) {
#     filePath <- file.path(tempdir(), filename)
#     write.csv(data, filePath, row.names = FALSE, quote = TRUE)
#     drop_upload(filePath, path = outputDir)
# }
# 
# 
# nmin <- 2
# nmax <- 14
# alg_deg <- list()
# for (n in nmin:nmax) {
#     data <- getFileFromDropbox(str_c("n=",n,".txt"))
#     i <- 0
#     dat_list <- list()
#     for (r in 1:(n-1)) {
#         m_low <- choose(n-r+1,2)
#         m_upp <- choose(n+1,2)-choose(r+1,2)
#         dat_vec <- rep(0,m_upp-m_low+1)
#         for (m in m_low:m_upp) {
#            i <- i+1
#            dat_vec[m-m_low+1] <- data[i,1]
#         }
#         names(dat_vec) <- str_c("m=",m_low:m_upp)
#         dat_list[[r]] <- dat_vec
#     }
#     names(dat_list) <- str_c("r=",1:(n-1))
#     alg_deg[[n-nmin+1]] <- dat_list
# }
# names(alg_deg) <- str_c("n=",nmin:nmax)
# save(alg_deg,file=str_c("data/algdegSDP_nmin=",nmin,"_nmax=",nmax,".RData"))


# (has been put in symconivol package)
# alg_deg <- function(m,n,r,nmin=2,nmax=14) {
#     if (any( m<0 | n<0 | r<0 | m<choose(n-r+1,2) | m>choose(n+1,2)-choose(r+1,2) )) {
#         return(0)
#     } else if ( n<nmin | n>nmax ){
#         return(stop("n outside available range"))
#     } else {
#         filename <- str_c("data/algdegSDP_nmin=",nmin,"_nmax=",nmax,".RData")
#         e <- new.env()
#         load(filename, envir = e)
#         n_ind <- str_c("n=",n)
#         r_ind <- str_c("r=",r)
#         m_ind <- str_c("m=",m)
#         return(unname(e$alg_deg[[n_ind]][[r_ind]][m_ind]))
#     }
# }
# 
# alg_deg(3,6,5)
# 
# n <- 6
# r <- 2
# pat <- pat_bnd(1,n)
# alg_deg((pat$d-pat$k_upp(r)):(pat$d-pat$k_low(r)),n,r)
# 
# 
plotsPB <- list()
j <- 0
# for (n in c(3,6,10)) {
for (n in 3:14) {
    beta <- 1
    j <- j+1

    pat <- pat_bnd(beta,n)
    d <- pat$d
    # construct matrix of algebraic degrees (easier to work with)
    AlgDeg <- matrix(0,d+1,n+1)
    AlgDeg[1,1] <- 1            # (this is just for the plots)
    AlgDeg[d+1,n+1] <- 1        # (this is just for the plots)
    for (r in 1:(n-1)) {
        for (k in pat$k_low(r):pat$k_upp(r)) {
            AlgDeg[k+1,r+1] <- alg_deg(d-k,n,r)
        }
    }
    indMax <- apply(AlgDeg, 2, max)
    dimMax <- apply(AlgDeg, 1, max)

    n_grid <- 0:n
    d_grid <- 0:d
    pts <- expand.grid(d_grid,n_grid)
    colnames(pts) <- c("k","r")
    pts$AlgDeg <- 0
    for (i in 1:dim(pts)[1]) {
        r <- pts$r[i]
        k <- pts$k[i]
        if (r<ceiling(n/2)) {
            r <- n-r
            k <- d-k
        }
        if (k==d & r==n)
            pts$AlgDeg[i] <- 1
        else if (k >= pat$k_low(r) & k <= pat$k_upp(r)) {
            # m <- pat$d-k
            pts$AlgDeg[i] <- AlgDeg[k+1,r+1] /
                                    # indMax[r+1]
                                    dimMax[k+1]
        }
    }
    # plot of Pataki bounds
    y_pat <- seq(0,1,length.out=1e2)
    pat_bd_low <- tibble(patX=d*y_pat^2, patY=n*y_pat)
    pat_bd_upp <- tibble(patX=d*(1-y_pat^2), patY=n*(1-y_pat))
    #
    tmp <- leigh(1e3)
    leighscurve <- tmp$table
    names(leighscurve) <- c("c","y")
    leighscurve$c <- leighscurve$c*n
    leighscurve$y <- leighscurve$y*d
    plotsPB[[j]] <- ggplot() +
        geom_line(data=leighscurve,aes(x=c,y=y), linetype=2, alpha=0.5) +
        geom_point(data=pts, aes(x=r, y=k, color=AlgDeg), alpha=0.5) +
        theme_bw() +
        xlab(TeX(sprintf("algebraic degree: $n=%d$", n, beta, n+beta*choose(n,2)))) +
        theme(legend.position="none",
              axis.title.y=element_blank(),
              axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")
        ) +
        scale_colour_gradient2() +
        geom_line(data=pat_bd_low, aes(patY, patX), linetype=2, alpha=0.5) +
        geom_line(data=pat_bd_upp, aes(patY, patX), linetype=2, alpha=0.5)
}
multiplot(plotlist = plotsPB, cols = 4)




plotsPB_Phi <- list()
j <- 0
for (beta in c(1,2,4)) {
    for (n in c(3,6,10)) {
        j <- j+1
        
        pat <- pat_bnd(beta,n)
        d <- pat$d
        n_grid <- 0:n
        d_grid <- 0:d
        pts <- expand.grid(d_grid,n_grid)
        colnames(pts) <- c("k","r")
        pts$PhiInd <- 0
        load(file=str_c("data/n=",n,"/PhiInd_beta=",beta,".RData"))
        for (i in 1:dim(pts)[1]) {
            r <- pts$r[i]
            k <- pts$k[i]
            if (r<ceiling(n/2)) {
                r <- n-r
                k <- d-k
            }
            if (k==d & r==n)
                pts$PhiInd[i] <- 1
            else if (k >= pat$k_low(r) & k <= pat$k_upp(r)) {
                colPhi <- str_c("r=",r,",s=0")
                pts$PhiInd[i] <- PhiInd[[colPhi]][k-pat$k_low(r)+1]/max(PhiInd[[colPhi]])
            }
        }
        # plot of Pataki bounds
        y_pat <- seq(0,1,length.out=1e2)
        pat_bd_low <- tibble(patX=d*y_pat^2, patY=n*y_pat)
        pat_bd_upp <- tibble(patX=d*(1-y_pat^2), patY=n*(1-y_pat))
        # 
        leighscurve <- leigh(1e3)$table
        leighscurve$c <- leighscurve$rho*n
        leighscurve$y <- leighscurve$kappa*d
        plotsPB_Phi[[j]] <- ggplot() + 
            geom_line(data=leighscurve,aes(x=c,y=y), linetype=2, alpha=0.5) +
            geom_point(data=pts, aes(x=r, y=k, color=PhiInd), alpha=0.5) +
            theme_bw() +
            xlab(TeX(sprintf("curvature measures: $n=%d,\\; \\beta = %d\\; (d=%d)$", n, beta, n+beta*choose(n,2)))) +
            theme(legend.position="none", 
                  axis.title.y=element_blank(),
                  axis.text.x=element_blank(), axis.text.y=element_blank(),
                  axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                  panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black")
            ) +
            scale_colour_gradient2() +
            geom_line(data=pat_bd_low, aes(patY, patX), linetype=2, alpha=0.5) +
            geom_line(data=pat_bd_upp, aes(patY, patX), linetype=2, alpha=0.5)
    }
}
multiplot(plotlist = plotsPB_Phi, cols = 3)


plotsPB_both <- list()
for (i in 1:3) {
    plotsPB_both[[i]] <- plotsPB_Phi[[i]]
    plotsPB_both[[i+3]] <- plotsPB[[i]]
}
multiplot(plotlist = plotsPB_both, cols = 2)





# USE ALGEBRAIC DEGREE TO PREDICT THE OTHER CURVE
n <- 14
beta <- 1

pat <- pat_bnd(beta,n)
d <- pat$d
# construct matrix of algebraic degrees (easier to work with)
AlgDeg <- matrix(0,d+1,n+1)
AlgDeg[1,1] <- 1            # (this is just for the plots)
AlgDeg[d+1,n+1] <- 1        # (this is just for the plots)
for (r in 1:(n-1)) {
    for (k in pat$k_low(r):pat$k_upp(r)) {
        AlgDeg[k+1,r+1] <- alg_deg(d-k,n,r)
    }
}
indMax <- apply(AlgDeg, 2, max)
dimMax <- apply(AlgDeg, 1, max)

n_grid <- 0:n
d_grid <- 0:d
pts <- expand.grid(d_grid,n_grid)
colnames(pts) <- c("k","r")
pts$AlgDeg <- 0
for (i in 1:dim(pts)[1]) {
    r <- pts$r[i]
    k <- pts$k[i]
    if (r<ceiling(n/2)) {
        r <- n-r
        k <- d-k
    }
    if (k==d & r==n)
        pts$AlgDeg[i] <- 1
    else if (k >= pat$k_low(r) & k <= pat$k_upp(r)) {
        # m <- pat$d-k
        pts$AlgDeg[i] <- AlgDeg[k+1,r+1] /
            # indMax[r+1]
            dimMax[k+1]
    }
}
# plot of Pataki bounds
y_pat <- seq(0,1,length.out=1e2)
pat_bd_low <- tibble(patX=d*y_pat^2, patY=n*y_pat)
pat_bd_upp <- tibble(patX=d*(1-y_pat^2), patY=n*(1-y_pat))
#
tmp <- leigh(1e3)
leighscurve <- tmp$table
names(leighscurve) <- c("c","y")
leighscurve$c <- leighscurve$c*n
leighscurve$y <- leighscurve$y*d
ggplot() +
    geom_line(data=leighscurve,aes(x=c,y=y), linetype=2, alpha=0.5) +
    geom_point(data=pts, aes(x=r, y=k, color=AlgDeg), alpha=0.5) +
    theme_bw() +
    xlab(TeX(sprintf("algebraic degree: $n=%d$", n, beta, n+beta*choose(n,2)))) +
    theme(legend.position="none",
          axis.title.y=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")
    ) +
    scale_colour_gradient2() +
    geom_line(data=pat_bd_low, aes(patY, patX), linetype=2, alpha=0.5) +
    geom_line(data=pat_bd_upp, aes(patY, patX), linetype=2, alpha=0.5)

################################################################################


mass_dim <- diag(1/rowSums(AlgDeg)) %*% AlgDeg
n_samp <- 1e3
samp_smooth <- tibble(k,r)
for (k in 0:d) {
    samp_smooth <- bind_rows(samp_smooth,
                             tibble( k=k, r=sample(0:n, n_samp, replace=TRUE, prob = mass_dim[k+1,]) ))
}
# dev.new(width=5, height=5)
P <- ggplot(samp_smooth, aes(x=r,y=k)) +
    # geom_point(alpha=0.01) +
    geom_smooth(se=FALSE, linetype=2, alpha=0.5, color="black", size=1) +
    geom_line(data=leighscurve,aes(x=c,y=y),size=1) +
    geom_line(data=pat_bd_low, aes(patY, patX),size=1) +
    geom_line(data=pat_bd_upp, aes(patY, patX),size=1) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
          axis.line.x=element_blank(), axis.line.y=element_blank())
ggsave("pat_leigh.png", plot=P, device="png",
       width=100, height=100, units="mm", dpi=300)

# mass_ind <- AlgDeg %*% diag(1/colSums(AlgDeg))
# n_samp <- 1e3
# samp_smooth <- tibble(k,r)
# for (r in 0:n) {
#     samp_smooth <- bind_rows(samp_smooth,
#                              tibble( k=sample(0:d, n_samp, replace=TRUE, prob = mass_ind[,r+1]), r=r ))
# }
# ggplot(samp_smooth, aes(x=r,y=k)) +
#     geom_point(alpha=0.01) +
#     geom_smooth(se=FALSE) +
#     geom_line(data=leighscurve,aes(x=c,y=y), linetype=2, alpha=0.5)























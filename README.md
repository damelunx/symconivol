
<!-- README.md is generated from README.Rmd. Please edit that file -->
symconivol:<img src="logo.png" width="125" height="125" align="right"/> An R package for curvature measures of symmetric cones
==============================================================================================================================

This R package provides functions for analyzing the curvature measures of symmetric cones.

Installation
------------

You can install `symconivol` from github with:

``` r
# install.packages("devtools")
devtools::install_github("damelunx/symconivol")
```

Motivation
----------

Intrinsic volumes form interesting and useful characteristics of convex cones. The cones of positive semidefinite matrices have an additional property, which is its decomposition into the rank strata. This decomposition lead to a decomposition of its intrinsic volumes into the so-called curvature measures. These curvature measures can be studied through the eigenvalue distribution of the [Gaussian orthogonal/unitary/symplectic ensemble](https://en.wikipedia.org/wiki/Random_matrix). The functions provided by this package facilitate studying this connection. See the accompanying vignette for such a study.

One application of the curvature measures is that they can be used for understanding the distribution of the rank of the solution of a random [semidefinite program](https://en.wikipedia.org/wiki/Semidefinite_programming). Concretely, if the semidefinite program is of the form
\begin{align*}
    \underset{X\in\mathcal{S}^n}{\text{min}} \quad & \langle C,X\rangle_{\mathcal{S}^n}
\\ \text{subject to} \quad & \langle A_k,X\rangle_{\mathcal{S}^n}=b_k ,\quad k=1,\ldots,m
\\ & X\succeq 0 ,
\end{align*}
that is, the SDP optimizes a linear functional over the intersection of the cone of positive semidefinite matrices aith an affine linear subspace of (generically) codimension *m*, then the rank of the solution (assuming that the SDP has a solution) can be predicted by:

``` r
n <- 100
m <- 17
pat <- pat_bnd(1,n)
d <- pat$d
pred_rank_sol <- round(n*mu()$lkup_rho(m/d))
```

See the corresponding [section](articles/curv_meas.html#appl_SDP) in the accompanying vignette for more details.

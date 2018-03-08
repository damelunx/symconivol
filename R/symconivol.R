#' symconivol: A package for curvature measures of symmetric cones
#'
#' The symconivol package provides functions for analyzing intrinsic volumes and
#' curvature measures of symmetric cones, as well as
#' the Gaussian orthogonal ensembles conditioned on the index
#' function.
#'
#' \code{symconivol} provides functions for analyzing intrinsic volumes and
#' curvature measures of symmetric cones (positive semidefinite
#' real/complex/quaternion matrices). These quantities can be estimated through
#' the eigenvalue distribution of the Gaussian ensembles conditioned on the
#' index, that is, the number of positive eigenvalues. The package provides
#' functions for sampling from these conditioned eigenvalue distributions
#' via Stan, and for reconstructing the curvature measures via MOSEK
#' (second-order program). The package also provides several convenient functions
#' for studiying these quantities, as well as a table of the algebraic degree
#' of semidefinite programming.
#' See the accompanying vignette for more information on how to use these functions.
#' 
#' @section Functions:
#' \itemize{
#'   \item \code{\link[symconivol]{curv_meas_exact}}: gives the exact curvature
#'         measures for \code{n=1,2,3}
#'   
#'   \item \code{\link[symconivol]{pat_bnd}}: provides the Pataki inequalities
#'         for given \code{beta} and \code{n}
#'   
#'   \item \code{\link[symconivol]{leigh}}: produces a table and lookup functions
#'         for Leigh's curve (see vignette for definition)
#'   
#'   \item \code{\link[symconivol]{constr_eigval}}: generates inputs for Stan
#'         (model string or external file) for sampling from the Gaussian
#'         orthogonal/unitary/symplectic ensemble conditioned on the index,
#'         the number of positive eigenvalues
#'   
#'   \item \code{\link[symconivol]{constr_eigval_to_bcbsq}}: converts a sample
#'         of eigenvalues produced by \code{constr_eigval} to a sample of the
#'         corresponding bivariate chi-bar-squared distribution
#'   
#'   \item \code{\link[symconivol]{prepare_em_cm}}: evaluates the sample data of
#'         the bivariate chi-bar-squared data (find the corresponding
#'         chi-squared density values)
#'   
#'   \item \code{\link[symconivol]{estim_em_cm}}: produces EM-type iterates to
#'         estimate the (normalized) curvature measures from a sample of the
#'         bivariate chi-bar-squared distribution
#'   
#'   \item \code{\link[symconivol]{indnorm_to_unnorm}}: estimates the unnormalized
#'         curvature measures from estimates of the dimension normalized ones and
#'         further estimates (see vignette for details)
#'   
#'   \item \code{\link[symconivol]{alg_deg}}: looks up the algebraic degree of
#'         semidefinite programming from a table
#' }
#'
#' @section See also:
#' \describe{
#'   \item{manual}{\href{https://damelunx.github.io/symconivol/}{damelunx.github.io/symconivol}}
#'   \item{sources}{\href{https://github.com/damelunx/symconivol}{github.com/damelunx/symconivol}}
#'   \item{vignette}{\href{https://damelunx.github.io/symconivol/articles/curv_meas.html}{Studying curvature measures of symmetric cones}:
#'         introduces curvature measures of symmetric cones, their relation to
#'         the Gaussian orthogonal/unitary/symplectic ensemble conditioned on
#'         the index function, explains the algorithms involved for estimating
#'         the curvature measures, gives some background and estimates involving
#'         limiting distributions and the algebraic degree of semidefinite programming}
#'   \item{vignette}{\href{https://damelunx.github.io/symconivol/articles/curv_meas_tech.html}{Studying curvature measures of symmetric cones - Technical details}:
#'         a technical note to accompany the other vignette to give the commands
#'         for producing some figures in the main vignette, which are
#'          not computed on the fly}
#' }
#'
#' @docType package
#' @name symconivol
NULL















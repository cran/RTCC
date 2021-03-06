#' RTCC: Detecting trait clustering in environmental gradients
#' with the Randomized Trait Community Clustering method
#'
#'
#' A set of functions which allows to determine if the observed traits
#' present clustering/overdispersion patterns on the observed samples,
#' and if so, to stablish if the observed pattern is linked to the effect
#' of an environmental gradient.
#'
#' @details The study of phenotypic similarities and differences within species along
#'   environmental gradients might be used as a powerful tool complementing
#'   taxon-based approaches when assesing the contribution of stochastic and
#'   deterministic processes in community assembly. For this, this package allows an
#'   easy implementation of a method for detecting clustering/overdispersion patterns
#'   along an environmental gradient (Triado-Margarit et al., 2019). A first function
#'   assesses if the observed traits exhibit a clustering/overdispersion pattern on the
#'   tested samples. If positive, two subsequent functions determine whether the
#'   observed pattern is linked to the effect of an environmental varible and its
#'   statistical significance.
#'
#' @section Data entry: The data consists on presence-absence observations
#'    along a measured environmental gradient and trait quantitative
#'    information of the observed organisms.
#'
#' @references Triado-Margarit, X., Capitan, J.A., Menendez-Serra, M. et al. (2019)
#'   A Randomized Trait Community Clustering approach to unveil consistent environmental
#'   thresholds in community assembly. \emph{ISME J} \bold{13}, 2681–2689 .
#'   \url{https://doi.org/10.1038/s41396-019-0454-4}
#'
#' @docType package
#' @name RTCC
NULL

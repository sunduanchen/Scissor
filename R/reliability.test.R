#' Reliability Significance Test
#'
#' Performs the reliability significance test to determine whether the inferred phenotype-to-cell
#' associations are reliable (statistical p-value less than 0.05) or are false positives.
#'
#' This statistical test can determine whether the inferred phenotype-to-cell associations are
#' reliable (statistical p-value less than 0.05) or are false positives.
#'
#' The motivation for the reliability significance test is: if the chosen single-cell and bulk data are not suitable for
#' the phenotype-to-cell associations, the correlations would be less informative and not well associated with the phenotype
#' labels. Thus, the corresponding prediction performance would be poor and not be significantly distinguishable from the
#' randomly permutated labels.
#'
#' The evaluation measurements used in the reliability significance test are the mean squared error (MSE) for linear regression,
#' the area under the ROC curve (AUC) for classification, and the concordance index (c-index) for Cox regression.
#'
#' @param X Scissor calculated correlation matrix for each pair of cells and bulk samples.
#' @param Y Scissor calculated response variable in the regression model.
#' @param network Scissor calculated cell-cell similarity network.
#' @param alpha Same parameter alpha used in Scissor.
#' @param family Same parameter family used in Scissor.
#' @param cell_num The number of the Scissor selected cells.
#' @param n Permutation times.
#' @param nfold The fold number in cross-validation.
#'
#' @return A list containing the following components:
#'   \item{statistic}{The test statistic.}
#'   \item{p}{The test p-value.}
#'   \item{Measurement_test_real}{The evaluation measurement on each fold calculated using the real label.}
#'   \item{Measurement_test_back}{A list with each component contains the evaluation measurements calculated using the permutated labels.}
#'
#' @import progress Matrix
#' @export
reliability.test <- function(X, Y, network, alpha, family = c("gaussian","binomial","cox"), cell_num, n = 100, nfold = 10){

    library(progress)
    library(Matrix)

    if (family == 'gaussian'){
        result <- test_lm(X, Y, network, alpha, cell_num, n, nfold)
    }
    if (family == 'binomial'){
        library(pROC)
        result <- test_logit(X, Y, network, alpha, cell_num, n, nfold)
    }
    if (family == 'cox'){
        library(survival)
        result <- test_cox(X, Y, network, alpha, cell_num, n, nfold)
    }
    return(result)
}

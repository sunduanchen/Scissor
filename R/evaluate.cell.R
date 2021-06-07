#' Cell level evaluations on the Scissor selected cells
#'
#' This function perfroms evaluations for each Scissor selected cell.
#'
#' This function uses two main strategies to evaluate each Scissor selected cell. Correlation check was first performed for each cell. Then,
#' we used the nonparametric bootstrap strategy to assess the sampling distribution of coefficient for each Scissor selected cell.
#'
#' @param Load_file File name for loading the preprocessed regression inputs.
#' @param Scissor_result Output variable from Scissor.
#' @param FDR_cutoff FDR cutoff in the correlation test for each pair of bulk sample and single-cell. The default value is 0.05.
#' @param bootstrap_n Bootstrap resampling times.
#'
#' @return This function returns a \code{data.frame} with rows are the Scissor selected cells. The columns contain the following attributes:
#'   \item{Mean correlation}{The mean correlation of a cell with all bulk samples.}
#'   \item{Correlation > 0}{The percentage of positive correlations of a cell with all bulk samples.}
#'   \item{Correlation < 0}{The percentage of negative correlations of a cell with all bulk samples.}
#'   \item{Significant Correlation}{The percentage of significant correlations (FDR < \code{FDR_cutoff}) of a cell with all bulk samples.}
#'   \item{Coefficient}{The regression coefficient beta calculated in Scissor.}
#'   \item{Beta 0\%}{The minimum value of bootstrap coefficient.}
#'   \item{Beta 25\%}{The lower quartile of bootstrap coefficient.}
#'   \item{Beta 50\%}{The median value of bootstrap coefficient.}
#'   \item{Beta 75\%}{TThe upper quartile of bootstrap coefficient.}
#'   \item{Beta 100\%}{The maximum value of bootstrap coefficient.}
#'   \item{Probability of zero}{The probability of a cell is not selected in bootstrap resamplings.}
#'
#' @import progress scales
#' @export
evaluate.cell <- function(Load_file, Scissor_result, FDR_cutoff = 0.05, bootstrap_n = 100){

    library(progress)
    library(scales)
    load(Load_file)     # X, Y, network, Expression_bulk, Expression_cell

    selected_cell <- c(Scissor_result$Scissor_pos, Scissor_result$Scissor_neg)
    m <- ncol(Expression_bulk)
    n <- length(selected_cell)
    evaluate_summary <- as.data.frame(matrix(0, n, 11, dimnames = list(selected_cell,
                        c("Mean correlation","Correlation > 0","Correlation < 0","Significant Correlation","Coefficient",
                          "Beta 0%","Beta 25%","Beta 50%","Beta 75%","Beta 100%","Probability of zero"))))

    print("|**************************************************|")
    print("Performing correlation check for each selected cell")
    pb1 <- progress_bar$new(total = m)
    cor_test_p <- matrix(0, m, n)
    for (i in 1:m){
        #pb1$tick()
        Sys.sleep(1 / 100)
        for (j in 1:n){
            cor_test_p[i,j] <- cor.test(Expression_bulk[,i], Expression_cell[,selected_cell[j]])$p.value
        }
    }
    cor_test_FDR <- matrix(p.adjust(as.numeric(cor_test_p), method = "fdr"), m)
    for (j in 1:n){
        evaluate_summary[j,1] <- mean(X[,selected_cell[j]])
        evaluate_summary[j,2] <- percent(sum(X[,selected_cell[j]] > 0)/m)
        evaluate_summary[j,3] <- percent(sum(X[,selected_cell[j]] < 0)/m)
        evaluate_summary[j,4] <- percent(sum(cor_test_FDR[,j] < FDR_cutoff)/m)
    }
    names(Scissor_result$Coefs) <- colnames(X)
    evaluate_summary[,5] <- Scissor_result$Coefs[selected_cell]
    cat("Finished!\n")

    print("|**************************************************|")
    print("Performing nonparametric bootstrap")
    pb2 <- progress_bar$new(total = bootstrap_n)
    beta_bootstrap <- matrix(0, bootstrap_n, n, dimnames = list(paste0("Bootstrap_",  1:bootstrap_n), selected_cell))
    for (i in 1:bootstrap_n){
        set.seed(i)
        alpha  <- Scissor_result$para$alpha
        lambda <- Scissor_result$para$lambda
        family <- Scissor_result$para$family
        index  <- sample(m, size = m, replace = TRUE)
        X2 <- X[index,]
        if (family == "cox"){
            Y2 <- Y[index,]
        }else{
            Y2 <- Y[index]
        }
        bootstrap_fit <- NULL
        while (is.null(bootstrap_fit$fit)){
            bootstrap_fit <- APML1(X2, Y2, family = family, penalty = "Net", alpha = alpha, Omega = network, lambda = lambda)
        }
        if (family == "binomial"){
            Coefs_tmp <- as.numeric(bootstrap_fit$Beta[2:(ncol(X2)+1),])
        }else{
            Coefs_tmp <- as.numeric(bootstrap_fit$Beta)
        }
        names(Coefs_tmp) <- colnames(X2)
        beta_bootstrap[i,] <- Coefs_tmp[selected_cell]
        #pb2$tick()
        Sys.sleep(1 / 100)
        if (i == bootstrap_n) cat("Finished!\n")
    }
    for (j in 1:n){
        tmp <- beta_bootstrap[beta_bootstrap[,j] != 0, j]
        evaluate_summary[j,6:10] <- fivenum(tmp)
        evaluate_summary[j,11] <- sum(beta_bootstrap[,j] == 0)/bootstrap_n
    }
    return(evaluate_summary)
}


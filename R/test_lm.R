test_lm <- function(X, Y, network, alpha, cell_num, n = 100, nfold = 10){

    set.seed(1)
    m <- nrow(X)
    index0 <- sample(cut(seq(m), breaks = nfold, labels = F))

    print("|**************************************************|")
    print("Perform cross-validation on X with true label")
    MSE_test_real <- NULL
    pb1 <- progress_bar$new(total = nfold)
    for (j in 1:nfold){
        c_index <- which(index0 == j)
        X_train <- X[-c_index,]
        Y_train <- Y[-c_index]
        fit <- NULL
        while (is.null(fit$fit)){
            set.seed(123)
            fit <- APML1(X_train, Y_train, family = "gaussian", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
        }
        index <- which.min(abs(fit$fit$nzero - cell_num))
        Coefs <- as.numeric(fit$Beta[,index])
        Cell1 <- Coefs[which(Coefs > 0)]
        Cell2 <- Coefs[which(Coefs < 0)]

        X_test <- X[c_index,]
        Y_test <- Y[c_index]
        MSE_test_real[j] <- mean((Y_test - X_test%*%Coefs)^2)

        #pb1$tick()
        Sys.sleep(1 / 100)
        if (j == nfold) cat("Finished!\n")
    }

    print("|**************************************************|")
    print("Perform cross-validation on X with permutated label")
    MSE_test_back <- list()
    pb2 <- progress_bar$new(total = n)
    for (i in 1:n){
        set.seed(i+100)
        MSE_test_back[[i]] <- matrix(0, nfold, 1, dimnames = list(paste0("Testing_", 1:nfold), "MSE"))
        Y2 <- Y[sample(m)]
        for (j in 1:nfold){
            c_index <- which(index0 == j)
            X_train <- X[-c_index,]
            Y_train <- Y2[-c_index]
            fit <- NULL
            while (is.null(fit$fit)){
                set.seed(123)
                fit <- APML1(X_train, Y_train, family = "gaussian", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
            }
            index <- which.min(abs(fit$fit$nzero - cell_num))
            Coefs <- as.numeric(fit$Beta[,index])
            Cell1 <- Coefs[which(Coefs > 0)]
            Cell2 <- Coefs[which(Coefs < 0)]

            X_test <- X[c_index,]
            Y_test <- Y2[c_index]
            MSE_test_back[[i]][j] <- mean((Y_test - X_test%*%Coefs)^2)
        }
        #pb2$tick()
        Sys.sleep(1 / 100)
        if (i == n) cat("Finished!\n")
    }
    statistic  <- mean(MSE_test_real)
    background <- NULL
    for (i in 1:n){
        background[i] <- mean(MSE_test_back[[i]][,1])
    }
    p <- sum(background < statistic)/n

    print(sprintf("Test statistic = %s", formatC(statistic, format = "f", digits = 3)))
    print(sprintf("Reliability significance test p = %s", formatC(p, format = "f", digits = 3)))

    return(list(statistic = statistic,
                p = p,
                MSE_test_real = MSE_test_real,
                MSE_test_back = MSE_test_back))
}


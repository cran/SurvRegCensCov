SurvRegCens <- function(time, event, CovariateNonCens = NULL, 
                        CovariateCens, Density, initial, conf.level = 0.95, 
                        intlimit = 10^-10, namCens = "VarCens", trace = 0, 
                        reltol = 10^-8){
    
    nam1 <- colnames(CovariateNonCens)
    if (is.null(CovariateNonCens) == FALSE){CovariateNonCens <- as.matrix(CovariateNonCens)}

    CovariateCens[is.na(CovariateCens[,1])==FALSE & is.na(CovariateCens[,2])==TRUE,2] <- Inf
    CovariateCens2 <- CovariateCens
    CovariateCens2[is.na(CovariateCens[, 1]) == FALSE & is.na(CovariateCens[, 2]) == FALSE & CovariateCens[, 1] == CovariateCens[, 2], 1] <- NA
    VectorR <- matrix(1, nrow = nrow(CovariateCens), ncol = 1)
    VectorR[is.na(CovariateCens[, 1]) == TRUE & is.na(CovariateCens[, 2]) == FALSE] <- 0
    VectorR[is.na(CovariateCens[, 1]) == FALSE & is.na(CovariateCens[, 2]) == FALSE & CovariateCens[, 1] != CovariateCens[, 2]] <- 0
    VectorLB <- CovariateCens2[, 1]

    signif.level <- 1 - conf.level
	  sample_size <- length(time)
    
	  result_value <- optim(initial, LoglikWeibullSurvRegCens, data_y = time, data_cov_noncens = CovariateNonCens, data_cov_cens = CovariateCens[, 2], density = Density, data_r_loglik=VectorR, data_lowerbound = VectorLB, data_delta_loglik = event, intlimit=intlimit, control = list(maxit = 5000, fnscale = -1, trace = trace, reltol = reltol), hessian = FALSE)
	
    if(is.null(CovariateNonCens) == FALSE){
        NumberNonCens <- ncol(as.matrix(CovariateNonCens))
        NamesBeta <- matrix(nrow = 1, ncol = NumberNonCens)
     }
    
	  Estimation_lambda <- result_value$par[1]
	  Estimation_gamma <- result_value$par[2]
    if(is.null(CovariateNonCens) == FALSE){Estimation_betaNonCens <- result_value$par[3:(length(result_value$par) - 1)]}
	  Estimation_betaCens <- result_value$par[length(result_value$par)]
    
    
    # standard errors
    HessianMatrix <- numDeriv::hessian(func = LoglikWeibullSurvRegCens, x = result_value$par, data_y = time, data_cov_noncens = CovariateNonCens, data_cov_cens = CovariateCens[,2], density = Density, data_r_loglik = VectorR, data_lowerbound = VectorLB, data_delta_loglik = event, intlimit=intlimit)
    SEs <- sqrt(-diag(solve(HessianMatrix)))

    CIlambda <-  cbind(Estimation_lambda - qnorm(1 - signif.level / 2) * SEs[1], Estimation_lambda + qnorm(1 - signif.level / 2) * SEs[1])
	  CIgamma <-  cbind(Estimation_gamma - qnorm(1 - signif.level / 2) * SEs[2], Estimation_gamma + qnorm(1 - signif.level / 2) * SEs[2])
	  if(is.null(CovariateNonCens) == FALSE){
        for(ii in 1:NumberNonCens){
            assign(paste("CIBetaNonCens", ii, sep = ""), cbind(Estimation_betaNonCens[ii] - qnorm(1 - signif.level / 2) * SEs[ii+2], Estimation_betaNonCens[ii] + qnorm(1 -     signif.level / 2) * SEs[ii+2]))
            NamesBeta[ii] <- paste("Beta", ii, sep = "")
        }
    }
	
    CIbetaCens <- cbind(Estimation_betaCens - qnorm(1 - signif.level / 2) * SEs[length(SEs)], Estimation_betaCens + qnorm(1 - signif.level / 2) * SEs[length(SEs)])
    
    if(is.null(CovariateNonCens) == FALSE){
        CIbetaNonCens <- matrix(nrow = 2 ,ncol = NumberNonCens)
        for(ii in 1:NumberNonCens){
            CIbetaNonCens[,ii] <- get(paste("CIBetaNonCens", ii, sep = ""))
    }
    
    if(NumberNonCens >= 2){
        Estimation_betaNonCens <- t(as.matrix(Estimation_betaNonCens))
        colnames(Estimation_betaNonCens) <- NamesBeta
        colnames(CIbetaNonCens) <- NamesBeta
        }
    }
    
  # collect results
	table <- as.data.frame(matrix(nrow = length(result_value$par), ncol = 8))
	tNamesRow <- matrix(nrow = nrow(table), ncol = 1)
	tNamesRow[1] <- "lambda"
	tNamesRow[2] <- "gamma"
    if(is.null(CovariateNonCens) == FALSE){
        if(ncol(as.matrix(CovariateNonCens)) > 1){
            tNamesRow[3:(nrow(tNamesRow) - 1)] <- nam1  ## paste("BetaNonCens", seq(from = 1, to = ncol(as.matrix(CovariateNonCens)), by = 1), sep = "")
        }
        
        if(ncol(as.matrix(CovariateNonCens)) == 1){tNamesRow[3] <- nam1}
    
    } 
    
	
    tNamesRow[nrow(tNamesRow)] <- namCens
    rownames(table) <- tNamesRow
    colnames(table) <- c("Estimator", "Std. Error", "CI.low", "CI.up", "p-value", "exp(Estimator)", "exp(CI.low)", "exp(CI.up)")
    
    table[, 1] <- result_value$par
    table[, 2] <- SEs
    table[1, 3:4] <- CIlambda
    table[2, 3:4] <- CIgamma
  
    if(is.null(CovariateNonCens) == FALSE){table[3:(2+ncol(CIbetaNonCens)), 3:4] <- t(CIbetaNonCens)}
    table[nrow(table), 3:4] <- CIbetaCens
    table[3:nrow(table), 6] <- exp(table[3:nrow(table), 1])
    table[3:nrow(table), 7:8] <- exp(table[3:nrow(table), 3:4])
    
    if(is.null(CovariateNonCens) == FALSE){
        p.value.BetaNonCens <- matrix(nrow = 1, ncol = NumberNonCens)
        for(ii in 1:NumberNonCens){
            test.BetaNonCens_ii <- (Estimation_betaNonCens[ii] - 0) / SEs[ii + 2]
            p.value.BetaNonCens[ii] <- 1 - pchisq(test.BetaNonCens_ii ^ 2, df = 1)
        }                                                                
        
        if(NumberNonCens >= 2){colnames(p.value.BetaNonCens) <- NamesBeta}
    }
    
	test.BetaCens <- (Estimation_betaCens - 0) / SEs[length(SEs)]
	p.value.BetaCens <- 1 - pchisq(test.BetaCens ^ 2, df = 1)
    
    if(is.null(CovariateNonCens) == FALSE){table[3:(nrow(table) - 1), 5] <- as.vector(p.value.BetaNonCens)}
     
	
    table[nrow(table), 5] <- p.value.BetaCens
    table[, 5] <- format.pval(table[, 5])

    percentage <- (1 - (sum(VectorR) / length(VectorR))) * 100
    value_loglik <- result_value$value

    ReturnList <- list(coeff = table, percent.cens = percentage, loglik= value_loglik, info.converg = result_value$convergence, info.converg.message = result_value$message)
    return(ReturnList)
}











#

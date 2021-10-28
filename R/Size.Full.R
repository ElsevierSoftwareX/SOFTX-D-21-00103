Size.Full <- function(factor.lev, interaction=FALSE, delta_type=1, delta=c(1, 0, 1), alpha=0.05, beta=0.2, maxsize=1000)
{
    if (!all(factor.lev > 0)) stop("The numbers of levels for each factor must be positive.\n")
    if (!any(delta_type == c(1, 2))) stop("The input argument 'delta_type' must be 1 or 2.\n")
    if (!interaction & !all(delta[-2] > 0)) stop("The effect sizes must be positive.\n")
    if (interaction & !all(delta > 0)) stop("The effect sizes must be positive.\n")
    if (alpha <= 0) stop("Type I error must be posive.\n")
    if (beta <= 0) stop("Type II error must be posive.\n")
  
    nfactor <- length(factor.lev)
    if (nfactor == 1 & interaction) {
        cat("There is just one factor in the model. Two-way interactions can't be included.\n")
        return(list(model=NULL, n=NULL, Delta=NULL))
    }
    
    prodlev <- prod(factor.lev)
    
    main_n <- two_n <- nn <- 0
    Delta <- NULL

    if (!interaction) {
        v <- factor.lev - 1
        tmpcoeff <- prodlev / factor.lev
    } else {
        v <- (factor.lev - 1) %*% t(factor.lev - 1)
        v <- c(factor.lev - 1, v[upper.tri(v, diag=FALSE)])      
        tmpcoeff <- prodlev / c(factor.lev, (factor.lev %*% t(factor.lev))[upper.tri((factor.lev) %*% t(factor.lev), diag=FALSE)])
        nfactorp1 <- nfactor + 1
        nterm <- length(v)
    }
    sumvp1 <- sum(v) + 1
    
    deltaratio13 <- delta[1]/delta[3]
    deltaratio23 <- delta[2]/delta[3]
    
    for (n in 2:maxsize) {
        v.denom <- prodlev * n - sumvp1
        coeff <- tmpcoeff * n 
        if (!interaction) {
            for (i in 1:nfactor)
                Delta[i] <- fsize(alpha, beta, v[i], v.denom, coeff[i], delta_type, FALSE)

            if (max(Delta) <= deltaratio13) {
                nn <- n  
                break
            }
        } else {
            if (main_n == 0) {
                for (i in 1:nfactor)
                    Delta[i] <- fsize(alpha, beta, v[i], v.denom, coeff[i], delta_type, FALSE)
                
                if (max(Delta[1:nfactor]) <= deltaratio13) 
                    main_n <- n
            }

            if (two_n == 0) {
                for (i in nfactorp1:nterm) 
                    Delta[i] <- fsize(alpha, beta, v[i], v.denom, coeff[i], delta_type, TRUE)
                
                if (max(Delta[nfactorp1:nterm]) <= deltaratio23) 
                    two_n <- n
            }
            
            if (main_n > 0 & two_n > 0) {
                nn <- max(main_n, two_n)   
                for (i in 1:nterm)
                    Delta[i] <- fsize(alpha, beta, v[i], v.denom, coeff[i], delta_type, FALSE)
  
                for (i in nfactorp1:nterm) 
                    Delta[i] <- fsize(alpha, beta, v[i], v.denom, coeff[i], delta_type, TRUE)
                break
            }
        }
    }

    if (nn == 0) {
        cat(paste0("The optimal sample size will be greater than ", maxsize, ".\n", "You may increase the input argument 'maxsize'.\n"))
        list(model=NULL, n=NULL, Delta=NULL)
    } else {
        model <- sizelist(nfactor, interaction)
        names(Delta) <- unlist(strsplit(model, "[+]"))
        list(model=model, n=nn, Delta=Delta)      
    }
}

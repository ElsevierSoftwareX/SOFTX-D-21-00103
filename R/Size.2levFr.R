Size.2levFr <- function(nfactor, nfraction, interaction=FALSE, delta_type=1, delta=c(1, 0, 1), alpha=0.05, beta=0.2, maxsize=1000)
{
    if (nfactor <= 0) stop("The number of factor must be positive.\n")
    if (nfraction <= 0) stop("The number of fraction must be positive.\n")
    if (!any(delta_type == c(1, 2))) stop("The input argument 'delta_type' must be 1 or 2.\n")
    if (!interaction & !all(delta[-2] > 0)) stop("The effect sizes must be positive.\n")
    if (interaction & !all(delta > 0)) stop("The effect sizes must be positive.\n")
    if (alpha <= 0) stop("Type I error must be posive.\n")
    if (beta <= 0) stop("Type II error must be posive.\n")  
  
    if (2^(nfactor - nfraction) - 1 - (nfactor + nfactor * (nfactor - 1) / 2) < 0 & interaction) {
        cat("Two-way interactions can't be included.\n") 
        return(list(model=NULL, n=NULL, Delta=NULL))
    }

    if (2^(nfactor - nfraction) - 1 - nfactor < 0 & !interaction) {
        cat("All of main effects can't be estimated.\n")   
        return(list(model=NULL, n=NULL, Delta=NULL))
    }
  
    factor.lev <- rep(2, nfactor)
    prodlev <- 2^nfactor
    prodfraclev <- 2^nfraction

    if (!interaction) {
        v <- factor.lev - 1     
        tmpcoeff <- prodlev / prodfraclev / factor.lev
        tmp2 <- -1 - nfactor
    } else {
        v <- (factor.lev - 1) %*% t(factor.lev - 1)
        v <- c(factor.lev - 1, v[upper.tri(v, diag=FALSE)])      
        tmpcoeff <- (prodlev / prodfraclev) /
            c(factor.lev, (factor.lev %*% t(factor.lev))[upper.tri(factor.lev %*% t(factor.lev), diag = FALSE)])
        tmp2 <- -1 - nfactor - nfactor * (nfactor - 1) / 2
    }

    deltaratio13 <- delta[1]/delta[3]
    deltaratio23 <- delta[2]/delta[3]
    
    nfactorp1 <- nfactor + 1
    tmp1 <- 2^(nfactor - nfraction)
 
    main_n <- two_n <- nn <- 0
    Delta <- NULL
    
    for (n in 2:maxsize) {
        v.denom <- tmp1 * n + tmp2
        coeff <- tmpcoeff * n 

        if (!interaction) {
            Delta[1] <- fsize(alpha, beta, v[1], v.denom, coeff[1], delta_type, FALSE)
            if (Delta[1] <= deltaratio13) {
                nn <- n  
                break
            }
        } else {
            if (main_n == 0) {
                Delta[1] <- fsize(alpha, beta, v[1], v.denom, coeff[1], delta_type, FALSE)
                if (Delta[1] <= deltaratio13) 
                    main_n <- n
            }
            if (two_n == 0) {
                Delta[2] <- fsize(alpha, beta, v[nfactorp1], v.denom, coeff[nfactorp1], delta_type, TRUE)
                if (Delta[2] <= deltaratio23) 
                    two_n <- n
            }
            if (main_n > 0 & two_n > 0) 
                nn <- max(main_n, two_n)   

            if (nn > 0) {
                if (!interaction) 
                    Delta[1] <- fsize(alpha, beta, v[1], v.denom, coeff[1], delta_type, FALSE)
                else {
                    Delta[1] <- fsize(alpha, beta, v[1], v.denom, coeff[1], delta_type, FALSE)
                    Delta[2] <- fsize(alpha, beta, v[nfactor+1], v.denom, coeff[nfactorp1], delta_type, TRUE)
                }
                break
            }       
        }
    }

    if (nn == 0) {
        cat(paste0("The optimal sample size will be greater than ", maxsize, ".\n", "You may increase the input argument 'maxsize'.\n"))
        list(model=NULL, n=NULL, Delta=NULL)
    } else {
        model <- sizelist(nfactor, interaction)
        terms <- "Main"
        if (interaction) terms <- c(terms, "Interaction")
        names(Delta) <- terms
        list(model=model, n=nn, Delta=Delta)      
    }
}
  

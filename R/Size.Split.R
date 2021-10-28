Size.Split <- function(whole.factor.lev, split.factor.lev, interaction=FALSE, 
                      delta_type=1, delta=c(1, 0, 1, 1), alpha=0.05, beta=0.2, maxsize=1000) 
{
    if (!all(whole.factor.lev > 0)) stop("The numbers of levels for each whole factor must be positive.\n")
    if (!all(split.factor.lev > 0)) stop("The numbers of levels for each split factor must be positive.\n")    
    if (!any(delta_type == c(1, 2))) stop("The input argument 'delta_type' must be 1 or 2.\n")
    if (!interaction & !all(delta[-2] > 0)) stop("The effect sizes must be positive.\n")
    if (interaction & !all(delta > 0)) stop("The effect sizes must be positive.\n")
    if (alpha <= 0) stop("Type I error must be posive.\n")
    if (beta <= 0) stop("Type II error must be posive.\n")
    
    nfactor <- length(whole.factor.lev)
    nsplit <- length(split.factor.lev)

    if (nfactor > 2 || nsplit > 2) {
        cat("The number of whole plot factors and split plot factors are up to 2 in the current package vesion.\n")
        return(list(model=NULL, n=NULL, Delta=NULL))
    } 
    
    prodwhole <- prod(whole.factor.lev)
    prodsplit <- prod(split.factor.lev)
    
    prodwholem1 <- prodwhole - 1
    prodsplitm1whole <- (prodsplit - 1) * prodwhole
    
    sqrtprodsplitp1 <- sqrt(prodsplit + 1)
    
    main_n <- two_n <- nn <- 0
    whole.Delta <- split.Delta <- NULL

    deltaratio13 <- delta[1]/delta[3]
    deltaratio14 <- delta[1]/delta[4]
    
    if (!interaction) {
        v.whole <- whole.factor.lev - 1
        tmpwcoeff <- prodwhole * prodsplit / whole.factor.lev
        
        v.split <- split.factor.lev - 1
        tmpscoeff <- prodwhole * prodsplit / split.factor.lev       

        wv_flag <- rep(0, nfactor)
        sv_flag <- rep(0, nsplit)
    } else {
        v.whole <- (whole.factor.lev - 1) %*% t(whole.factor.lev - 1)
        v.whole <- c(whole.factor.lev - 1, v.whole[upper.tri(v.whole, diag=FALSE)])
        tmpwcoeff <- prodwhole * prodsplit / c(whole.factor.lev, (whole.factor.lev %*% t(whole.factor.lev))[upper.tri((whole.factor.lev) %*% t(whole.factor.lev), diag=FALSE)])
          
        v.split <- (split.factor.lev - 1) %*% t(split.factor.lev - 1)
        v.split.temp <- (whole.factor.lev - 1) %*% t(split.factor.lev - 1)
        v.split <- c(split.factor.lev - 1, v.split[upper.tri(v.split, diag=FALSE)], as.vector(t(v.split.temp)))
        
        tmpscoeff <- prodwhole * prodsplit / c(split.factor.lev, (split.factor.lev %*% t(split.factor.lev))[upper.tri(split.factor.lev %*% t(split.factor.lev), diag=FALSE)], 
                                               as.vector(t(whole.factor.lev %*% t(split.factor.lev))))
        
        nfactorp1 <- nfactor + 1
        nsplitp1 <- nsplit + 1

        deltaratio23 <- delta[2]/delta[3]
        deltaratio24 <- delta[2]/delta[4]
        
        wv_flag <- c(rep(0, nfactor), rep(1, nfactor * (nfactor - 1) / 2))
        sv_flag <- c(rep(0, nsplit), rep(1, nsplit * (nsplit - 1) / 2), rep(1, nfactor * nsplit))
    }
    sumv.whole <- sum(v.whole)
    sumv.split <- sum(v.split)
    nvwhole <- length(v.whole)
    nvsplit <- length(v.split)

    for (n in 2:maxsize) {
        # whole
        v.whole.denom <- prodwholem1 * n - sumv.whole
        c.whole <- tmpwcoeff * n      
        
        # split
        v.split.denom <- prodsplitm1whole * n - sumv.split
        c.split <- tmpscoeff * n
        
        if (!interaction) {
            for (i in 1:nfactor) 
                whole.Delta[i] <- fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i], delta_type, FALSE) * sqrtprodsplitp1
        
            for (i in 1:nsplit)
                split.Delta[i] <- fsize(alpha, beta, v.split[i], v.split.denom, c.split[i], delta_type, FALSE)
        
            if (max(whole.Delta) <= deltaratio13 & max(split.Delta) <= deltaratio14)
                nn <- n
        } else {
            if (main_n == 0) {
                for (i in 1:nfactor)
                    whole.Delta[i] <- fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i], delta_type, FALSE) * sqrtprodsplitp1
                  
                for (i in 1:nsplit)
                    split.Delta[i] <- fsize(alpha, beta, v.split[i], v.split.denom, c.split[i], delta_type, FALSE)
                  
                if (max(whole.Delta[1:nfactor]) <= deltaratio13 & max(split.Delta[1:nsplit]) <= deltaratio14) 
                    main_n <- n
            }

            if (two_n == 0) {
                for (i in nsplitp1:nvsplit)
                    split.Delta[i] <- fsize(alpha, beta, v.split[i], v.split.denom, c.split[i], delta_type, TRUE)
                  
                if (nfactor > 1) {
                    for (i in nfactorp1:nvwhole)
                        whole.Delta[i] <- fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i], delta_type, TRUE) * sqrtprodsplitp1
                    
                    if (max(whole.Delta[nfactorp1:nvwhole]) <= deltaratio23 & max(split.Delta[nsplitp1:nvsplit]) <= deltaratio24) 
                        two_n <- n
                } else if (max(split.Delta[nsplitp1:nvsplit]) <= deltaratio24) 
                    two_n <- n 
            }
          
            if(main_n > 0 & two_n > 0) 
                nn <- max(main_n, two_n)
        }

        if (nn > 0) {
            for (i in 1:nvwhole)
                whole.Delta[i] <- fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i], delta_type, wv_flag[i]) * sqrtprodsplitp1

            for (i in 1:nvsplit)
                split.Delta[i] <-  fsize(alpha, beta, v.split[i], v.split.denom, c.split[i], delta_type, sv_flag[i])

            Delta <- c(whole.Delta, split.Delta)
            break
        }
    }

    if (nn == 0) {
        cat(paste0("The optimal sample size will be greater than ", maxsize, ".\n", "You may increase the input argument 'maxsize'.\n"))
        list(model=NULL, n=NULL, Delta=NULL)
    } else {
        model <- sizelist.split(nfactor, nsplit, interaction)
        names(Delta) <- unlist(strsplit(model, "[+]"))
        list(model=model, n=nn, Delta=Delta)     
    }
}

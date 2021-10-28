plots.Full <- function(factor.lev, interaction=FALSE, delta_type=1, delta=c(1, 0, 1), deltao=NULL, alpha=0.05, beta=0.2, type=1, maxsize=1000)
{
    if (!is.null(deltao) & type == 3)
        if (deltao <= 0) stop("The minimal detectable standardized effect size must be positive.\n")
    if (is.null(deltao) & type == 3)
        stop("When 'type=3', 'deltao' must be specified to draw the plot of power versus the 
             sample size acquiring the minimal detectable standardized effect size given by 'delato'.") 
    if (!any(type == c(1, 2, 3))) stop("The type of graph must be 1, 2 or 3.\n")

    nfactor <- length(factor.lev)
    prodlev <- prod(factor.lev)
    
    FF <- Size.Full(factor.lev, interaction, delta_type, delta, alpha, beta, maxsize)
    if (is.null(FF$model))
        return()
    
    n.choose <- FF$n
    
    terms <- unlist(strsplit(FF$model, "[+]"))
    ndata <- 1000
    
    power <- round(seq(0, 1, length.out=ndata+1), 3)
    
    if (!interaction) {
        v <- factor.lev - 1
        tmpcoeff <- prodlev / factor.lev
        v_flag <- rep(0, nfactor)
    } else {
        v <- (factor.lev - 1) %*% t(factor.lev - 1)
        v <- c(factor.lev - 1, v[upper.tri(v, diag=FALSE)])      
        tmpcoeff <- prodlev / c(factor.lev, (factor.lev %*% t(factor.lev))[upper.tri((factor.lev) %*% t(factor.lev), diag=FALSE)])
        v_flag <- c(rep(0, nfactor), rep(1, nfactor * (nfactor - 1) / 2))
    }
    nterm <- length(v)
    sumvp1 <- sum(v) + 1

    uniqueindex <- NULL
    uv <- unique(v[1:nfactor])
    if (length(uv) != nfactor) {
        for (i in 1:(length(uv))) {
            tmpindex <- (1:nfactor)[uv[i] == v[1:nfactor]]
            terms[tmpindex] <- paste(terms[tmpindex], collapse=", ")         
            uniqueindex <- c(uniqueindex, tmpindex[1])
        }

        if (interaction) {
            uv <- unique(v[(nfactor+1):nterm])
            for (i in 1:(length(uv))) {
                tmpindex <- ((nfactor+1):nterm)[uv[i] == v[(nfactor+1):nterm]]
                terms[tmpindex] <- paste(terms[tmpindex], collapse=", ")
                uniqueindex <- c(uniqueindex, tmpindex[1])
            }
        }
    } else 
        uniqueindex <- 1:nterm
    nunique <- length(uniqueindex)
    
    factor_type <- rep(terms[uniqueindex], each=ndata)        

    Delta1 <- power1 <- NULL
    if (type == 1) {
        v.denom <- prodlev * n.choose - sumvp1
        coeff <- tmpcoeff * n.choose
        tmp1 <- tmp2 <- tmp3 <- NULL 
        for (j in uniqueindex) 
            for (ind in 1:ndata) 
                if (alpha + 1 - power[ind] < 0.9999 & 1 - power[ind] > 0.0001 & 1 - power[ind] < 0.9999) {
                    tmp1 <- c(tmp1, terms[j])
                    tmp2 <- c(tmp2, fsize(alpha, 1-power[ind], v[j], v.denom, coeff[j], delta_type, v_flag[j])) 
                    tmp3 <- c(tmp3, power[ind])
                }
 
        data1 <- data.frame(factor_type=tmp1, Delta1=tmp2, power1=tmp3)
        data1$factor_type <- as.character(data1$factor_type)
                            
        if (nunique != 1) {
            tmpindex <- 1:nunique       
            for (j in 1:(nunique-1)) {
                tmp0 <- data1[data1[, 1] == terms[uniqueindex][j], 2]
                for (k in (j+1):nunique) 
                    if (sum(abs((data1[data1[, 1] == terms[uniqueindex][k], 2] - tmp0))) < 0.001 * stats::sd(tmp0))
                        tmpindex[k] <- tmpindex[j]
            }   
            uniqueindex2 <- unique(tmpindex)    
            if (length(uniqueindex2) != nunique) 
                for (i in 1:length(uniqueindex2)) {
                    tmpindex0 <- uniqueindex[uniqueindex2[i] == tmpindex]
                    for (j in 1:length(tmpindex0))
                        data1[data1[, 1] == terms[tmpindex0][j], 1] <- paste(terms[tmpindex0], collapse=", ") 
                }              
        }    

        data1$factor_type <- factor(data1$factor_type, levels = unique(data1$factor_type)) #terms[uniqueindex])
        gr <- ggplot2::ggplot(data1, ggplot2::aes(x=Delta1, y=power1, color=factor_type)) + 
            ggplot2::geom_line(size=1.5) + ggplot2::labs(title = "Power vs Delta") + 
            ggplot2::ylab("Power") + ggplot2::xlab("Delta") + 
            ggplot2::geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = c(1.0, 1.5), linetype = "dashed")
    } else if (type == 2) {
        tmp2 <- NULL
        for (j in uniqueindex)
            for (n in 1:ndata + 1) {
                v.denom <- prodlev * n - sumvp1
                coeff <- tmpcoeff * n 
                tmp2 <- c(tmp2, fsize(alpha, beta, v[j], v.denom, coeff[j], delta_type, v_flag[j]))       
            }
        data2 <- data.frame(n=1:ndata + 1, factor_type=factor_type, Delta1=tmp2)
        data2$factor_type <- as.character(data2$factor_type)
        
        if (nunique != 1) {
            tmpindex <- 1:nunique       
            for (j in 1:(nunique-1)) {
                tmp0 <- data2[data2[, 2] == terms[uniqueindex][j], 3]
                for (k in (j+1):nunique) 
                    if (sum(abs((data2[data2[, 2] == terms[uniqueindex][k], 3] - tmp0))) < 0.001 * stats::sd(tmp0))
                        tmpindex[k] <- tmpindex[j]
            }   
            uniqueindex2 <- unique(tmpindex)    
            if (length(uniqueindex2) != nunique)
                for (i in 1:length(uniqueindex2)) {
                    tmpindex0 <- uniqueindex[uniqueindex2[i] == tmpindex]
                    for (j in 1:length(tmpindex0))
                        data2[data2[, 2] == terms[tmpindex0][j], 2] <- paste(terms[tmpindex0], collapse=", ") 
                }              
        }  
        
        data2$factor_type <- factor(data2$factor_type, levels = unique(data2$factor_type)) #terms[uniqueindex])
        gr <- ggplot2::ggplot(data2[1:ndata + 1 <= 2*n.choose, ], ggplot2::aes(x=n, y=Delta1, color=factor_type)) +
            ggplot2::geom_line(size=1.5) + ggplot2::labs(title = "Delta vs Sample size") + 
            ggplot2::ylab("Delta") + ggplot2::xlab("Sample size") + 
            ggplot2::geom_hline(yintercept = c(1.0, 1.5), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = n.choose, linetype = "dashed")
    } else if (type == 3) {
        deltao2 <- deltao^2
        tmp3 <- NULL
        for (j in uniqueindex)
            for (n in 1:ndata + 1) {
                v.denom <- prodlev * n - sumvp1
                coeff <- tmpcoeff * n 
                tmp3 <- c(tmp3, round(1-stats::pf(stats::qf(1-alpha, v[j], v.denom), v[j], v.denom, 
                                            ncp = deltao2 * coeff[j] * ifelse(delta_type == 1, v[j], ifelse(v_flag[j] == 1, 1, 0.5))), 3))      
            }
        data3 <- data.frame(n=1:ndata + 1, factor_type=factor_type, power1=tmp3)    
        data3$factor_type <- as.character(data3$factor_type)
        
        if (nunique != 1) {
            tmpindex <- 1:nunique       
            for (j in 1:(nunique-1)) {
                tmp0 <- data3[data3[, 2] == terms[uniqueindex][j], 3]
                for (k in (j+1):nunique) 
                    if (sum(abs((data3[data3[, 2] == terms[uniqueindex][k], 3] - tmp0))) < 0.001 * stats::sd(tmp0))
                        tmpindex[k] <- tmpindex[j]
            }   
            uniqueindex2 <- unique(tmpindex)    
            if (length(uniqueindex2) != nunique)
                for (i in 1:length(uniqueindex2)) {
                    tmpindex0 <- uniqueindex[uniqueindex2[i] == tmpindex]
                    for (j in 1:length(tmpindex0))
                        data3[data3[, 2] == terms[tmpindex0][j], 2] <- paste(terms[tmpindex0], collapse=", ") 
                }              
        }  
        
        data3$factor_type <- factor(data3$factor_type, unique(data3$factor_type)) #terms[uniqueindex])
        gr <- ggplot2::ggplot(data3[1:ndata + 1 <= 2*n.choose, ], ggplot2::aes(x=n, y=power1, color=factor_type)) +
            ggplot2::geom_line(size=1.5) + ggplot2::labs(title = "Power vs Sample size") + 
            ggplot2::ylab("Power") + ggplot2::xlab("Sample size") + 
            ggplot2::geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = n.choose, linetype = "dashed")        
    }    
    
    # The palette with grey:
    cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # The palette with black:
    cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    gr + ggplot2::theme(text = ggplot2::element_text(size = 20)) + ggplot2::scale_colour_manual(values=cbp2)   
}

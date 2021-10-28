plots.Split <- function(whole.factor.lev, split.factor.lev, interaction=FALSE, 
                        delta_type=1, delta=c(1, 0, 1, 1), deltao=NULL, alpha=0.05, beta=0.2, type=1, maxsize=1000)
{
    if (!is.null(deltao) & type == 3)
        if (deltao <= 0) stop("The minimal detectable standardized effect size must be positive.\n")
    if (is.null(deltao) & type == 3)
        stop("When 'type=3', 'deltao' must be specified to draw the plot of power versus the 
              sample size acquiring the minimal detectable standardized effect size given by 'delato'.") 
    if (!any(type == c(1, 2, 3))) stop("The type of graph must be 1, 2 or 3.\n")
  
    nfactor <- length(whole.factor.lev)
    nsplit <- length(split.factor.lev)

    prodwhole <- prod(whole.factor.lev)
    prodsplit <- prod(split.factor.lev)
    
    prodwholem1 <- prodwhole - 1
    prodsplitm1whole <- (prodsplit - 1) * prodwhole
    
    sqrtprodsplitp1 <- sqrt(prodsplit + 1)

    FF <- Size.Split(whole.factor.lev, split.factor.lev, interaction, delta_type, delta, alpha, beta, maxsize)
    if (is.null(FF$model))
        return()
    
    n.choose <- FF$n
    
    terms <- unlist(strsplit(FF$model, "[+]"))
    ndata <- 1000
    
    power <- round(seq(0, 1, length.out=ndata+1), 3)
    
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

        wv_flag <- c(rep(0, nfactor), rep(1, nfactor * (nfactor - 1) / 2))
        sv_flag <- c(rep(0, nsplit), rep(1, nsplit * (nsplit - 1) / 2), rep(1, nfactor * nsplit))
    }
    
    sumv.whole <- sum(v.whole)
    sumv.split <- sum(v.split)
    nvwhole <- length(v.whole)
    nvsplit <- length(v.split)

    uniquewindex <- uniquesindex <- NULL
    uv <- unique(v.whole[1:nfactor])
    if (length(uv) != nfactor) {
        for (i in 1:length(uv)) {
            tmpindex <- (1:nfactor)[uv[i] == v.whole[1:nfactor]]
            terms[tmpindex] <- paste(terms[tmpindex], collapse=", ")     
            uniquewindex <- c(uniquewindex, tmpindex[1])
        }

        if (interaction) {
            uv <- unique(v.whole[(nfactor+1):nvwhole])
            for (i in 1:length(uv)) {
                tmpindex <- ((nfactor+1):nvwhole)[uv[i] == v.whole[(nfactor+1):nvwhole]]
                terms[tmpindex] <- paste(terms[tmpindex], collapse=", ")  
                uniquewindex <- c(uniquewindex, tmpindex[1])                
            }
        }
    } else 
        uniquewindex <- 1:nvwhole
    
    uv <- unique(v.split[1:nsplit])
    if (length(uv) != nsplit) {
        for (i in 1:length(uv)) {
            tmpindex <- (nvwhole + 1:nsplit)[uv[i] == v.split[1:nsplit]]
            terms[tmpindex] <- paste(terms[tmpindex], collapse=", ")     
            uniquesindex <- c(uniquesindex, (tmpindex-nvwhole)[1])          
        }
    } else 
      uniquesindex <- 1:nsplit 

    if (interaction) {
        uv <- unique(v.split[(nsplit+1):nvsplit])     
        if (length(uv) != length(v.split[(nsplit+1):nvsplit])) 
            for (i in 1:length(uv)) {
                tmpindex <- (nvwhole+(nsplit+1):nvsplit)[uv[i] == v.split[(nsplit+1):nvsplit]]
                terms[tmpindex] <- paste(terms[tmpindex], collapse=", ")       
                uniquesindex <- c(uniquesindex, (tmpindex-nvwhole)[1])
            }
        else
            uniquesindex <- c(uniquesindex, (nsplit+1):nvsplit)
    }
    
    uniqueindex <- c(uniquewindex, uniquesindex + nvwhole)
    nunique <- length(uniqueindex)
    
    factor_type <- rep(terms[uniqueindex], each=ndata)     
     
    Delta1 <- power1 <- NULL    
    if (type == 1) {
        v.whole.denom <- prodwholem1 * n.choose - sumv.whole
        c.whole <- tmpwcoeff * n.choose
        v.split.denom <- prodsplitm1whole * n.choose - sumv.split
        c.split <- tmpscoeff * n.choose
        
        tmp1 <- tmp2 <- tmp3 <- NULL 
        for (i in uniquewindex) 
            for (ind in 1:ndata) 
                if (alpha + 1 - power[ind] < 0.9999 & 1 - power[ind] > 0.0001 & 1 - power[ind] < 0.9999) {
                    tmp1 <- c(tmp1, terms[i])
                    tmp2 <- c(tmp2, fsize(alpha, 1-power[ind], v.whole[i],  v.whole.denom, c.whole[i], delta_type, wv_flag[i]) * sqrtprodsplitp1)
                    tmp3 <- c(tmp3, power[ind])
                }
        for (i in uniquesindex) 
            for (ind in 1:ndata) 
                if (alpha + 1 - power[ind] < 0.9999 & 1 - power[ind] > 0.0001 & 1 - power[ind] < 0.9999) {
                    tmp1 <- c(tmp1, terms[nvwhole+i])
                    tmp2 <- c(tmp2, fsize(alpha, 1-power[ind], v.split[i], v.split.denom, c.split[i], delta_type, sv_flag[i]))
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
            ggplot2::geom_line(size=1.5) + ggplot2::labs(title = "Powervs Delta") +
            ggplot2::ylab("Power") + ggplot2::xlab("Delta") + 
            ggplot2::geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = c(1.0, 1.5), linetype = "dashed")
    } else if (type == 2) {
        tmp2 <- NULL
        for (i in uniquewindex) 
            for (n in 1:ndata + 1) {
                # whole
                v.whole.denom <- prodwholem1 * n - sumv.whole
                c.whole <- tmpwcoeff * n      
                
                # split
                v.split.denom <- prodsplitm1whole * n - sumv.split
                c.split <- tmpscoeff * n
                
                tmp2 <- c(tmp2, fsize(alpha, beta, v.whole[i], v.whole.denom, c.whole[i], delta_type, wv_flag[i]) * sqrtprodsplitp1)
            }
        for (i in uniquesindex)  
            for (n in 1:ndata + 1) {
                # whole
                v.whole.denom <- prodwholem1 * n - sumv.whole
                c.whole <- tmpwcoeff * n      
              
                # split
                v.split.denom <- prodsplitm1whole * n - sumv.split
                c.split <- tmpscoeff * n
                tmp2 <- c(tmp2, fsize(alpha, beta, v.split[i], v.split.denom, c.split[i], delta_type, sv_flag[i]))
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
      deltao22 <- deltao2 / (prodsplit + 1)
      tmp3 <- NULL
      for (i in uniquewindex) 
          for (n in 1:ndata + 1) {
            # whole
            v.whole.denom <- prodwholem1 * n - sumv.whole
            c.whole <- tmpwcoeff * n      
            
            # split
            v.split.denom <- prodsplitm1whole * n - sumv.split
            c.split <- tmpscoeff * n
            tmp3 <- c(tmp3, 1 - stats::pf(stats::qf(1-alpha, v.whole[i], v.whole.denom), v.whole[i], v.whole.denom,
                        ncp = deltao22 * c.whole[i] * ifelse(delta_type == 1, v.whole[i],
                            ifelse(wv_flag[i] == 1, 1, 0.5)))) # (Delta^2)*(c*nu1)
        }
        for (i in uniquesindex) 
            for (n in 1:ndata + 1) {
            # whole
            v.whole.denom <- prodwholem1 * n - sumv.whole
            c.whole <- tmpwcoeff * n      
        
            # split
            v.split.denom <- prodsplitm1whole * n - sumv.split
            c.split <- tmpscoeff * n
            tmp3 <- c(tmp3, 1 - stats::pf(stats::qf(1-alpha, v.split[i], v.split.denom), v.split[i], v.split.denom,
                        ncp = deltao2 * c.split[i] * ifelse(delta_type == 1, v.split[i],
                            ifelse(sv_flag[i] == 1, 1, 0.5)))) #(Delta^2)*(c*nu1)
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

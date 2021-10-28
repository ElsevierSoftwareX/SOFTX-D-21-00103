sizelist.split <- function(nfactor, nsplit, interaction=FALSE)
{
    terms <- terms1 <- letters[1:nfactor]
    terms2 <- LETTERS[nfactor + 1:nsplit]

    if (interaction & nfactor != 1)  
        for (i in 1:(nfactor - 1)) 
            for (j in 1:(nfactor - i)) 
                terms <- c(terms, paste0(terms1[i], "*" , terms1[i+j]))
    
    terms <- c(terms, terms2)

    if (interaction & nsplit != 1)      
        for (i in 1:(nsplit - 1)) 
            for (j in 1:(nsplit - i)) 
                terms <- c(terms, paste0(terms2[i], "*" , terms2[i+j]))
  
    if (interaction)   
        for (i in 1:nfactor) 
            for (j in 1:nsplit) 
                terms <- c(terms, paste0(terms1[i], "*" , terms2[j]))  
  
    paste(terms, collapse="+")
}


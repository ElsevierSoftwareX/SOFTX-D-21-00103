sizelist <- function(nfactor, interaction=FALSE)
{
    if (nfactor == 1 & interaction)
        stop("There is just one factor in the model. Two-way interactions can't be included.")
    
    terms <- LETTERS[1:nfactor]

    if (interaction) 
        for (i in 1:(nfactor - 1)) 
            for (j in 1:(nfactor - i)) 
                terms <- c(terms, paste0(terms[i], "*" , terms[i+j]))

    paste(terms, collapse="+")
}

fsize <- function(alpha, beta, nu1, nu2, coeff, delta_type=1, interaction=FALSE)
{
    fl <- fpow::ncparamF(alpha, beta, nu1, nu2) / 2
    if (delta_type == 1) {
        Delta <- sqrt(2 * fl / (coeff * nu1)) #ncp <- Delta^2 * coeff * nu1
    } else if (delta_type == 2 & !interaction) {
        Delta <- sqrt(4 * fl / coeff) #ncp <- Delta^2 * coeff /2
    } else if (delta_type == 2 & interaction) 
        Delta <- sqrt(2 * fl / coeff) #ncp <- Delta^2 * coeff /2
    
    Delta
}

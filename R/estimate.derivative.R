### Code:

estimate.derivative <- function(dt, alpha_hat, tau = 3,
                                parameter = "auxiliary", 
                                h = length(unique(dt[, id]))^{-1/6},
                                parallelize.Z = 50, ...) {

    if (alpha_hat < h) {
        h <- alpha_hat*3/4
    }
    
    psi_plus  <- alpha.est.fun(dt, alpha = alpha_hat + h, tau = tau,
                               parameter = parameter, parallelize.Z = parallelize.Z,
                               ...)
    
    psi_minus <- alpha.est.fun(dt, alpha = alpha_hat - h, tau = tau,
                               parameter = parameter, parallelize.Z = parallelize.Z,
                               ...)

    return((psi_plus["est"] - psi_minus["est"]) / (2 * h))
        
}

######################################################################


### Code:

calibration.fun <- function(theta = 0.5, rho = NULL, browse = FALSE,
                            tau = 3, parallelize.Z = 50, output.eic = FALSE,
                            test.derivatives = FALSE, ...) {

    if (length(rho)>0) {
        est.aux.1 <- alpha.est.fun(dt, verbose = TRUE, tau = 3, alpha = 1,
                                   parameter = "auxiliary", parallelize.Z = parallelize.Z,
                                   output.eic = TRUE)
        theta <- rho*est.aux.1[["estimate"]]["est"]
    }

    if (browse) browser()
    
    est.alpha <- estimate.alpha.tmle(dt, theta = theta, tau = tau, parallelize.Z = parallelize.Z)

    est.deriv.z <- estimate.derivative(dt, est.alpha$alpha.hat, tau = tau,
                                       h = 0.3*length(unique(dt[, id]))^{-1/6}*sqrt(mean(est.alpha$eic^2)),
                                       parallelize.Z = parallelize.Z)

    if (test.derivatives) {
        est.deriv.n.z <- estimate.derivative(dt, est.alpha$alpha.hat, tau = tau,
                                             h =  length(unique(dt[, id]))^{-1/6},
                                             parallelize.Z = parallelize.Z)
        est.deriv.logn.z <- estimate.derivative(dt, est.alpha$alpha.hat, tau = tau,
                                                h = length(unique(dt[, id]))^{-1/6}/log(length(unique(dt[, id]))),
                                                parallelize.Z = parallelize.Z)
    }
   
    if (length(rho)>0) {
        alpha.eic <- (1/est.deriv.z)*(rho*est.aux.1$eic - est.alpha$eic)
        if (test.derivatives) {
            alpha.n.eic <- (1/est.deriv.n.z)*(rho*est.aux.1$eic - est.alpha$eic)
            alpha.logn.eic <- (1/est.deriv.logn.z)*(rho*est.aux.1$eic - est.alpha$eic)            
        }
    } else {
        alpha.eic <- (-1/est.deriv.z)*est.alpha$eic
        if (test.derivatives) {
            alpha.n.eic <- (-1/est.deriv.n.z)*est.alpha$eic
            alpha.logn.eic <- (-1/est.deriv.logn.z)*est.alpha$eic            
        }
    }
    
    alpha.se <- sqrt(mean(alpha.eic^2)/length(dt[, unique(id)]))

    if (test.derivatives) {
        alpha.n.se <- sqrt(mean(alpha.n.eic^2)/length(dt[, unique(id)]))
        alpha.logn.se <- sqrt(mean(alpha.logn.eic^2)/length(dt[, unique(id)]))            
    }
    
    target.est <- alpha.est.fun(dt, verbose = TRUE, tau = tau, alpha = est.alpha$alpha.hat,
                                output.eic = TRUE, parallelize.Z = parallelize.Z)

    target.se.crude <- target.est[["estimate"]]["se"]

    est.deriv.1 <- estimate.derivative(dt, est.alpha$alpha.hat,
                                       parameter = "target", tau = tau,
                                       h = 0.3*length(unique(dt[, id]))^{-1/6}*sqrt(mean(est.alpha$eic^2)),
                                       parallelize.Z = parallelize.Z)

    if (test.derivatives) {
        est.deriv.n.1 <- estimate.derivative(dt, est.alpha$alpha.hat,
                                             parameter = "target", tau = tau,
                                             h = length(unique(dt[, id]))^{-1/6},
                                             parallelize.Z = parallelize.Z)
        est.deriv.logn.1 <- estimate.derivative(dt, est.alpha$alpha.hat,
                                                parameter = "target", tau = tau,
                                                h = length(unique(dt[, id]))^{-1/6}/log(length(unique(dt[, id]))),
                                                parallelize.Z = parallelize.Z)
    }

    target.se <- sqrt(mean((target.est[["eic"]] + est.deriv.1*alpha.eic)^2)/length(dt[, unique(id)]))

    if (test.derivatives) {
        target.n.se <- sqrt(mean((target.est[["eic"]] + est.deriv.n.1*alpha.n.eic)^2)/length(dt[, unique(id)]))
        target.logn.se <- sqrt(mean((target.est[["eic"]] + est.deriv.logn.1*alpha.logn.eic)^2)/length(dt[, unique(id)]))            
    }

    out.estimate <- c(alpha.est = est.alpha$alpha.hat, alpha.se = alpha.se,
                      est.deriv.z = est.deriv.z,
                      est.alpha.converged = est.alpha$converged, est.alpha.dist = est.alpha$dist, est.alpha.cn = est.alpha$c.n,
                      target.est = as.numeric(target.est[["estimate"]]["est"]), target.se = target.se,
                      est.deriv.1 = est.deriv.1,
                      target.se.crude = as.numeric(target.se.crude))

    if (test.derivatives) {
        out.estimate <- c(out.estimate,
                          alpha.n.se = alpha.n.se,
                          alpha.logn.se = alpha.logn.se,
                          est.deriv.n.z = est.deriv.n.z,
                          est.deriv.logn.z = est.deriv.logn.z,
                          target.n.se = target.n.se,
                          target.logn.se = target.logn.se,
                          est.deriv.n.1 = est.deriv.n.1,
                          est.deriv.logn.1 = est.deriv.logn.1)
    }

    if (output.eic) {
        return(list(estimate = out.estimate,
                    alpha.eic = alpha.eic, # est.alpha$eic, #
                    target.eic = target.est[["eic"]],
                    eic = target.est[["eic"]] + est.deriv.1*alpha.eic))
    } else {
        return(out.estimate)
    }

}

######################################################################


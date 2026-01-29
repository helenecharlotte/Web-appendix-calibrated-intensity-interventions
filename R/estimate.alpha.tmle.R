### Code:

estimate.alpha.tmle <- function(dt,
                                theta,
                                tau = 3,
                                c_n = length(unique(dt[, id]))^{-1/2}/log(length(unique(dt[, id]))),
                                alpha_init = 1,
                                expand_up = 1.25,
                                expand_down = 0.8,
                                max_iter = 100,
                                alpha_min = 1e-3,
                                alpha_max = 100,
                                parallelize.Z = 50,
                                effect_alpha = NULL, 
                                ...) {

    # storage
    alpha_hist   <- numeric()
    psi_hist     <- numeric()
    se_hist      <- numeric()
    eic          <- numeric()

    # helper with caching
    if (length(effect_alpha) == 0) {
        eval_psi <- function(alpha) {
            idx <- which(abs(alpha_hist - alpha) < 1e-12)
            if (length(idx) > 0) return(psi_hist[idx[1]])

            psi <- alpha.est.fun(dt, tau = tau, alpha = alpha, parameter = "auxiliary",
                                 parallelize.Z = parallelize.Z, output.eic = TRUE, ...)
            alpha_hist <<- c(alpha_hist, alpha)
            psi_hist   <<- c(psi_hist, psi[["estimate"]]["est"])
            se_hist    <<- c(se_hist, psi[["estimate"]]["se"])
            eic        <<- psi[["eic"]]
            psi[["estimate"]]["est"]
        }
    } else {
         eval_psi <- function(alpha) {
            idx <- which(abs(alpha_hist - alpha) < 1e-12)
            if (length(idx) > 0) return(psi_hist[idx[1]])

            psi <- as.numeric(effect_alpha(tau = tau, alpha = alpha)["truth_z"])
            alpha_hist <<- c(alpha_hist, alpha)
            psi_hist   <<- c(psi_hist, psi)
            psi
         }
    }

    # step 0
    alpha_prev <- NA
    alpha_curr <- alpha_init
    psi_curr   <- eval_psi(alpha_curr)

    if (abs(psi_curr - theta) <= c_n) {
        return(list(
            alpha.hat = alpha_curr,
            converged = 1,
            dist = abs(psi_curr - theta),
            c.n = c_n,
            grid = data.frame(alpha = alpha_hist,
                              psi   = psi_hist,
                              se    = se_hist),
            eic = eic
        ))
    }

    # step 1
    if (psi_curr < theta - c_n) {
        alpha_next <- max(0, expand_up * alpha_curr)
    } else {
        alpha_next <- max(0, expand_down * alpha_curr)
    }

    # main loop
    for (m in 1:max_iter) {

        alpha_next <- max(min(alpha_next, alpha_max), alpha_min)

        psi_next   <- eval_psi(alpha_next)

        if (abs(psi_next - theta) <= c_n) break

        # update rule
        if (psi_next < theta - c_n) {
            alpha_new <- if (!is.na(alpha_prev) && alpha_next < alpha_curr)
                             (alpha_next + alpha_curr) / 2
                         else
                             expand_up * alpha_next
        } else {
            alpha_new <- if (!is.na(alpha_prev) && alpha_next > alpha_curr)
                             (alpha_next + alpha_curr) / 2
                         else
                             expand_down * alpha_next
        }

        if (alpha_new < 0) alpha_new <- 0

        alpha_prev <- alpha_curr
        alpha_curr <- alpha_next
        alpha_next <- alpha_new
    }

    # ordered grid
    ord <- order(alpha_hist)

    dist <- abs(psi_hist - theta)
    best_idx <- which.min(dist)
    alpha_best <- alpha_hist[best_idx]

    if (min(dist) <= c_n) {
        converged <- TRUE
    } else {
        converged <- FALSE
        eval_psi(alpha_best)
    }

    if (length(effect_alpha) == 0) {
        list(
            alpha.hat = alpha_best,
            converged = converged,
            dist = min(dist),
            c.n = c_n,
            grid = data.frame(alpha = alpha_hist[ord],
                              psi   = psi_hist[ord],
                              se    = se_hist[ord]),
            eic = eic
        )
    } else {
        list(
            alpha.hat = alpha_best,
            converged = converged,
            dist = min(dist),
            c.n = c_n,
            grid = data.frame(alpha = alpha_hist[ord],
                              psi   = psi_hist[ord])
        )
    }
}

######################################################################


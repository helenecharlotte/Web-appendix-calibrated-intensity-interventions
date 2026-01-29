### output.int.effect.alpha.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Dec 19 2025 (14:01) 
## Version: 
## Last-Updated: Jan 28 2026 (21:15) 
##           By: Helene
##     Update #: 558
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


library(data.table)
library(survival)
library(ggplot2)
library(xtable)
library(gridExtra)
library(zoo)
library(nleqslv)
library(foreach)
library(doParallel)
library(parallel)

#-------------------------------------------------------------------------------------------#
## alpha-specific output 
#-------------------------------------------------------------------------------------------#

get.output.fun <- function(alpha = 1, target = "target",
                           alpha2 = NULL, target2 = target,
                           M = 500, N = 1000, firstM = M, 
                           misspecify = c(), setting = 0,
                           truncate.weights = 0,
                           print.output = FALSE,
                           output.weights = FALSE,
                           browse = FALSE) {

    print(file.info(paste0("./output/",
                           ifelse(setting == 1, "setting1-", ""),
                           "estimation-", target, "-alpha", alpha,
                           ifelse("1" %in% misspecify, "-misspecify-1", ""),
                           ifelse("z" %in% misspecify, "-misspecify-z", ""),
                           ifelse("l" %in% misspecify, "-misspecify-l", ""),
                           ifelse(truncate.weights>0, paste0("-truncate-weights-", truncate.weights), ""),
                           "-M", M, "-N", N, 
                           ".rds"))$mtime)

    get.output <- readRDS(file = paste0("./output/",
                                        ifelse(setting == 1, "setting1-", ""),
                                        "estimation-", target, "-alpha", alpha,
                                        ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                        ifelse("z" %in% misspecify, "-misspecify-z", ""),
                                        ifelse("l" %in% misspecify, "-misspecify-l", ""),
                                        ifelse(truncate.weights>0, paste0("-truncate-weights-", truncate.weights), ""),
                                        "-M", M, "-N", N, 
                                        ".rds"))

    if (any(sapply(get.output, function(rep) "eic" %in% names(rep)))) {
        get.output.eic <- lapply(get.output, function(rep) {
            if ("eic" %in% names(rep)) {
                rep[["eic"]]
            } else {
                NULL
            }
        })
        get.output <- lapply(get.output, function(rep) {
            if ("eic" %in% names(rep)) {
                rep[["estimate"]]
            } else {
                rep
            }
        })
    }

    if (output.weights) {
        which.weights <- grep(".weight", names(get.output[[1]]), value = TRUE)
        out.weights <- data.table(do.call("cbind", lapply(which.weights, function(which.weight) {
            sapply(get.output, function(est) est[which.weight])
        })))
        names(out.weights) <- names(sapply(which.weights, function(which.weight) {
            gsub(strsplit(which.weight, ".weight")[[1]][2], "", which.weight)
        }))
        return(out.weights)
    }
   
    truth <- readRDS(file = paste0("./output/",
                                   ifelse(setting == 1, "setting1-", ""),
                                   "true-estimands-alpha-", alpha,
                                   ".rds"))[ifelse(target == "target", "truth_1", "truth_z")]

    if (browse) browser()

    if (length(alpha2)>0) {
        get.output2 <- readRDS(file = paste0("./output/",
                                             ifelse(setting == 1, "setting1-", ""),
                                             "estimation-", target2, "-alpha", alpha2,
                                             ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                             ifelse("z" %in% misspecify, "-misspecify-z", ""),
                                             ifelse("l" %in% misspecify, "-misspecify-l", ""),
                                             ifelse(truncate.weights>0, paste0("-truncate-weights-", truncate.weights), ""),
                                             "-M", M, "-N", N, 
                                             ".rds"))

        if ("eic" %in% names(get.output2[[1]])) {
            get.output2.eic <- lapply(get.output2, function(rep) rep[["eic"]])
            get.output2 <- lapply(get.output2, function(rep) rep[["estimate"]])
            eic <- TRUE
        } else {
            eic <- FALSE
        }
        
        truth2 <- readRDS(file = paste0("./output/",
                                        ifelse(setting == 1, "setting1-", ""),
                                        "true-estimands-alpha-", alpha2,
                                        ".rds"))[ifelse(target2 == "target", "truth_1", "truth_z")]

        est.output <- sapply(get.output, function(rep) rep["est"])
        est.output2 <- sapply(get.output2, function(rep) rep["est"])

        se.output <- sapply(get.output, function(rep) rep["se"])
        se.output2 <- sapply(get.output2, function(rep) rep["se"])

        est.diff <- est.output-est.output2

        if (eic) {
            se.diff <- sapply(1:length(get.output2.eic), function(rep.ii) {
                sqrt(mean((get.output.eic[[rep.ii]] - get.output2.eic[[rep.ii]])^2/N))
            })
        } else {
            se.diff <- sqrt(se.output^2+se.output2^2)
        }

        ci.lower <- est.diff-1.96*se.diff
        ci.upper <- est.diff+1.96*se.diff

        print(summary(ci.lower))
        print(summary(ci.upper))

        return(data.table(comparison = paste0("alpha = ", alpha, " vs alpha = ", alpha2),
                          est.diff = mean(est.diff), truth = truth-truth2,
                          sd.diff = sd(est.diff), se.diff = mean(se.diff),
                          cov = mean(ci.lower <= truth-truth2 & ci.upper >= truth-truth2),
                          oracle.cov = mean(est.diff-1.96*sd(est.diff) <= truth-truth2 &
                                            est.diff+1.96*sd(est.diff) >= truth-truth2),
                          ci.lower = mean(ci.lower), ci.upper = mean(ci.upper),
                          power = mean((ci.lower>0 & ci.upper>0) |
                                       (ci.lower<0 & ci.upper<0))))
    }

    print(length(get.output))

    if (firstM<length(get.output)) {
        get.output <- get.output[1:firstM]
        print(length(get.output))
    }

    if (print.output) {
        print(get.output)
    }

    g.est <- mean(g.est.output <- sapply(get.output, function(rep) rep["g.est"]))
    est <- mean(est.output <- sapply(get.output, function(rep) rep["est"]))
    sd <- sd(sapply(get.output, function(rep) rep["est"]))
    se <- mean(se.output <- sapply(get.output, function(rep) rep["se"]))
    hist(sapply(get.output, function(rep) rep["est"])); abline(v = truth)

    cov <- mean(est.output-1.96*se.output <= truth & est.output+1.96*se.output >= truth)
    oracle.cov <- mean(est.output-1.96*sd <= truth & est.output+1.96*sd >= truth)

    mse.fun <- function(x) mean((x - as.numeric(truth))^2)

    out <- data.table(alpha = alpha, g.est = g.est, est = est, truth = as.numeric(truth),
                      bias = est-as.numeric(truth), mse = mse.fun(est.output), sd = sd, se = se,
                      cov = cov, oracle.cov = oracle.cov)

    if ("q50.weight.50%" %in% names(get.output[[1]])) {
        out <- cbind(out,
                     data.table(
                         q50.weight = mean(sapply(get.output, function(est) est["q50.weight.50%"])),
                         q60.weight = mean(sapply(get.output, function(est) est["q60.weight.60%"])),
                         q70.weight = mean(sapply(get.output, function(est) est["q70.weight.70%"])),
                         q80.weight = mean(sapply(get.output, function(est) est["q80.weight.80%"])),
                         q90.weight = mean(sapply(get.output, function(est) est["q90.weight.90%"])),
                         q100.weight = mean(sapply(get.output, function(est) est["q100.weight.100%"])),
                         q50.weight.max = max(sapply(get.output, function(est) est["q50.weight.50%"])),
                         q60.weight.max = max(sapply(get.output, function(est) est["q60.weight.60%"])),
                         q70.weight.max = max(sapply(get.output, function(est) est["q70.weight.70%"])),
                         q80.weight.max = max(sapply(get.output, function(est) est["q80.weight.80%"])),
                         q90.weight.max = max(sapply(get.output, function(est) est["q90.weight.90%"])),
                         q100.weight.max = max(sapply(get.output, function(est) est["q100.weight.100%"]))))
    }
    
    if ("q50.weight" %in% names(get.output[[1]])) {
        out <- cbind(out,
                     data.table(
                         q50.weight = mean(sapply(get.output, function(est) est["q50.weight"])),
                         q60.weight = mean(sapply(get.output, function(est) est["q60.weight"])),
                         q70.weight = mean(sapply(get.output, function(est) est["q70.weight"])),
                         q80.weight = mean(sapply(get.output, function(est) est["q80.weight"])),
                         q90.weight = mean(sapply(get.output, function(est) est["q90.weight"])),
                         q100.weight = mean(sapply(get.output, function(est) est["q100.weight"])),
                         q50.weight.max = max(sapply(get.output, function(est) est["q50.weight"])),
                         q60.weight.max = max(sapply(get.output, function(est) est["q60.weight"])),
                         q70.weight.max = max(sapply(get.output, function(est) est["q70.weight"])),
                         q80.weight.max = max(sapply(get.output, function(est) est["q80.weight"])),
                         q90.weight.max = max(sapply(get.output, function(est) est["q90.weight"])),
                         q100.weight.max = max(sapply(get.output, function(est) est["q100.weight"]))))
    }
    
    return(out)
}

#-------------------------------------------------------------------------------------------#
## calibration output 
#-------------------------------------------------------------------------------------------#

get.calibration.output.fun <- function(theta = 0.5, rho = NULL,
                                       setting = 0, 
                                       M = 500, N = 1000,
                                       firstM = M,
                                       plot = FALSE, browse = FALSE,
                                       alpha2 = NULL,
                                       theta2 = NULL,
                                       rho2 = NULL,
                                       sign.contrast = 1,
                                       print.all.derivatives = output.all.derivatives,
                                       print.derivatives = print.all.derivatives | output.derivatives,
                                       output.derivatives = output.all.derivatives,
                                       output.all.derivatives = FALSE) {

    print(file.info(paste0("./output/",
                           ifelse(setting == 1, "setting1-", ""),
                           "estimation-calibration-",
                           ifelse(length(rho)>0,
                                  paste0("rho-", rho),
                                  paste0("theta-", theta)),
                           "-M", M, "-N", N, 
                           ".rds"))$mtime)

    est.calibration <- readRDS(file = paste0("./output/",
                                             ifelse(setting == 1, "setting1-", ""),
                                             "estimation-calibration-",
                                             ifelse(length(rho)>0,
                                                    paste0("rho-", rho),
                                                    paste0("theta-", theta)),
                                             "-M", M, "-N", N, 
                                             ".rds"))

    if (firstM<length(est.calibration)) {
        est.calibration <- est.calibration[1:firstM]
    }

    M1 <- length(est.calibration)

    if ("est.alpha.converged" %in% names(est.calibration[[M1]])) {
        print(mean(sapply(est.calibration, function(est) est["est.alpha.converged"])))
    }

    if ("eic" %in% names(est.calibration[[M1]])) {
        est.calibration.eic <- lapply(est.calibration, function(rep) rep[grep("eic", names(rep))])
        est.calibration <- lapply(est.calibration, function(rep) rep[["estimate"]])
        eic <- TRUE
    } else {
        eic <- FALSE
    }

    if ("est.deriv.z.est" %in% names(est.calibration[[M1]])) {
        deriv.z <- mean(deriv.z.est <- sapply(est.calibration, function(rep) rep[["est.deriv.z.est"]]))
        deriv.1 <- mean(deriv.1.est <- sapply(est.calibration, function(rep) rep[["est.deriv.1.est"]]))

        true.derivatives <- readRDS(file = paste0("./output/",
                                                  ifelse(setting == 1, "setting1-", ""),
                                                  "true-derivatives-",
                                                  ifelse(length(rho)>0,
                                                         paste0("rho-", rho),
                                                         paste0("theta-", theta)),
                                                  ".rds"))

        if (print.derivatives) {
        
            message("---------------------------------")
            print(paste0("deriv.z = ", deriv.z))
            print(paste0("true.deriv.z = ", true.derivatives["true.derivative.z"]))
            message("---------------------------------")
            print(paste0("sd(deriv.z) = ", sd(deriv.z.est)))
            message("---------------------------------")
            print(paste0("deriv.1 = ", deriv.1))
            print(paste0("true.deriv.1 = ", true.derivatives["true.derivative.1"]))
            message("---------------------------------")
            print(paste0("sd(deriv.1) = ", sd(deriv.1.est)))
            message("---------------------------------")

            if (output.derivatives) {
                out.derivatives <- data.table(deriv.1 = deriv.1, true.deriv.1 = true.derivatives["true.derivative.1"],
                                              deriv.z = deriv.z, true.deriv.z = true.derivatives["true.derivative.z"])
            }

        }

        if ("est.deriv.n.z.est" %in% names(est.calibration[[M1]]) & print.all.derivatives) {

            message("---------------------------------")

            deriv.n.z <- mean(deriv.n.z.est <- sapply(est.calibration, function(rep) rep[["est.deriv.n.z.est"]]))
            deriv.n.1 <- mean(deriv.n.1.est <- sapply(est.calibration, function(rep) rep[["est.deriv.n.1.est"]]))

            print(paste0("deriv.n.z = ", deriv.n.z))
            print(paste0("true.deriv.n.z = ", true.derivatives["true.derivative.z"]))
            message("---------------------------------")
            print(paste0("sd(deriv.n.z) = ", sd(deriv.n.z.est)))
            message("---------------------------------")
            print(paste0("deriv.n.1 = ", deriv.n.1))
            print(paste0("true.deriv.n.1 = ", true.derivatives["true.derivative.1"]))
            message("---------------------------------")
            print(paste0("sd(deriv.n.1) = ", sd(deriv.n.1.est)))
            message("---------------------------------")

            if (output.all.derivatives) {
                out.derivatives <- cbind(out.derivatives,
                                         data.table(deriv.n.1 = deriv.n.1,
                                                    deriv.n.z = deriv.n.z))
            }

        }

        if ("est.deriv.logn.z.est" %in% names(est.calibration[[M1]]) & print.all.derivatives) {

            message("---------------------------------")

            deriv.logn.z <- mean(deriv.logn.z.est <- sapply(est.calibration, function(rep) rep[["est.deriv.logn.z.est"]]))
            deriv.logn.1 <- mean(deriv.logn.1.est <- sapply(est.calibration, function(rep) rep[["est.deriv.logn.1.est"]]))

            print(paste0("deriv.logn.z = ", deriv.logn.z))
            print(paste0("true.deriv.logn.z = ", true.derivatives["true.derivative.z"]))
            message("---------------------------------")
            print(paste0("sd(deriv.logn.z) = ", sd(deriv.logn.z.est)))
            message("---------------------------------")
            print(paste0("deriv.logn.1 = ", deriv.logn.1))
            print(paste0("true.deriv.logn.1 = ", true.derivatives["true.derivative.1"]))
            message("---------------------------------")
            print(paste0("sd(deriv.logn.1) = ", sd(deriv.logn.1.est)))
            message("---------------------------------")

            if (output.all.derivatives) {
                out.derivatives <- cbind(out.derivatives,
                                         data.table(deriv.logn.1 = deriv.logn.1,
                                                    deriv.logn.z = deriv.logn.z))
            }

         }
    }

    if (browse) browser()

    print(length(est.calibration))

    true.effect <- readRDS(file = paste0("./output/",
                                         ifelse(setting == 1, "setting1-", ""),
                                         "true-estimand-",
                                         ifelse(length(rho)>0,
                                                paste0("rho-", rho),
                                                paste0("theta-", theta)),
                                         ".rds"))

    mean(est.alpha <- sapply(est.calibration, function(est) est["alpha.est"]))
    mean(est.target <- sapply(est.calibration, function(est) est[grep("target.est", names(est), value = TRUE)]))


    mean(se.alpha <- sapply(est.calibration, function(est) est["alpha.se"]))
    sd(est.alpha)
    mean(se.target <- sapply(est.calibration, function(est) est["target.se"]))
    sd(est.target)

    cov.alpha <- mean(est.alpha - 1.96*se.alpha <= true.effect["alpha"] &
                      est.alpha + 1.96*se.alpha >= true.effect["alpha"])

    oracle.cov.alpha <- mean(est.alpha - 1.96*sd(est.alpha) <= true.effect["alpha"] &
                             est.alpha + 1.96*sd(est.alpha) >= true.effect["alpha"])

    cov.target <- mean(est.target - 1.96*se.target <= true.effect["truth_1"] &
                       est.target + 1.96*se.target >= true.effect["truth_1"])

    oracle.cov.target <- mean(est.target - 1.96*sd(est.target) <= true.effect["truth_1"] &
                              est.target + 1.96*sd(est.target) >= true.effect["truth_1"])

    if ("est.deriv.n.z.est" %in% names(est.calibration[[M1]]) & print.all.derivatives) {

        message("---------------------------------")

        se.n.alpha <- sapply(est.calibration, function(est) est["alpha.n.se"])
        se.n.target <- sapply(est.calibration, function(est) est["target.n.se"])

        print(paste0("cov.alpha.n = ", cov.n.alpha <- mean(est.alpha - 1.96*se.n.alpha <= true.effect["alpha"] &
                                                           est.alpha + 1.96*se.n.alpha >= true.effect["alpha"])))

        print(paste0("cov.target.n = ", cov.n.target <- mean(est.target - 1.96*se.n.target <= true.effect["truth_1"] &
                                                             est.target + 1.96*se.n.target >= true.effect["truth_1"])))
        
        message("---------------------------------")

    }

    if ("est.deriv.logn.z.est" %in% names(est.calibration[[M1]]) & print.all.derivatives) {

        message("---------------------------------")

        se.logn.alpha <- sapply(est.calibration, function(est) est["alpha.logn.se"])
        se.logn.target <- sapply(est.calibration, function(est) est["target.logn.se"])

        print(paste0("cov.alpha.logn = ", cov.logn.alpha <- mean(est.alpha - 1.96*se.logn.alpha <= true.effect["alpha"] &
                                                           est.alpha + 1.96*se.logn.alpha >= true.effect["alpha"])))

        print(paste0("cov.target.logn = ", cov.logn.target <- mean(est.target - 1.96*se.logn.target <= true.effect["truth_1"] &
                                                                   est.target + 1.96*se.logn.target >= true.effect["truth_1"])))
        
        message("---------------------------------")

    }

    if (("est.deriv.z.est" %in% names(est.calibration[[M1]])) & print.derivatives) {

        if (FALSE) {
            se.alpha.true.deriv <- sapply(est.calibration.eic, function(est) {
                sqrt(mean(((-1/true.derivatives["true.derivative.z"])*est[["alpha.eic"]])^2)/N)
            })

            se.target.true.deriv <- sapply(est.calibration.eic, function(est) {
                sqrt(mean((est[["target.eic"]]+true.derivatives["true.derivative.1"]*
                           (-1/true.derivatives["true.derivative.z"])*est[["alpha.eic"]])^2)/N)
            })

            message("---------------------------------")
        
            print(paste0("se.alpha.true.deriv = ", mean(se.alpha.true.deriv)))
            print(paste0("se.alpha = ", mean(se.alpha)))
            print(paste0("se.target.true.deriv = ", mean(se.target.true.deriv)))
            print(paste0("se.target = ", mean(se.target)))

            message("---------------------------------")
        }
        
        alpha.eic <- lapply(1:length(est.calibration.eic), function(est.ii) {
            (est.calibration.eic[[est.ii]][["eic"]] - est.calibration.eic[[est.ii]][["target.eic"]]) /
                (est.calibration[[est.ii]]["est.deriv.1.est"] / est.calibration[[est.ii]]["est.deriv.z.est"])
        })

        se.alpha.true.deriv <- sapply(alpha.eic, function(est) {
            sqrt(mean(((1/true.derivatives["true.derivative.z"])*est)^2)/N)
        })

        se.target.true.deriv <- sapply(1:length(est.calibration.eic), function(est.ii) {
            sqrt(mean((est.calibration.eic[[est.ii]][["target.eic"]]+
                       true.derivatives["true.derivative.1"]*
                       (1/true.derivatives["true.derivative.z"])*alpha.eic[[est.ii]])^2)/N)
        })

        message("---------------------------------")
        
        print(paste0("se.alpha.true.deriv = ", mean(se.alpha.true.deriv)))
        print(paste0("se.alpha = ", mean(se.alpha)))
        print(paste0("se.target.true.deriv = ", mean(se.target.true.deriv)))
        print(paste0("se.target = ", mean(se.target)))

        message("---------------------------------")

        print(paste0("cov.alpha.true.deriv = ", cov.alpha.true.deriv <-
                                                    mean(est.alpha-1.96*se.alpha.true.deriv <= true.effect["alpha"] &
                                                         est.alpha+1.96*se.alpha.true.deriv >= true.effect["alpha"])))

        print(paste0("cov.target.true.deriv = ", cov.target.true.deriv <-
                                                     mean(est.target-1.96*se.target.true.deriv <= true.effect["truth_1"] &
                                                          est.target+1.96*se.target.true.deriv >= true.effect["truth_1"])))

        message("---------------------------------")

        if (output.derivatives) {
            out.derivatives <- cbind(out.derivatives,
                                     data.table(cov.alpha.true.deriv = cov.alpha.true.deriv, cov.alpha = cov.alpha,
                                                cov.target.true.deriv = cov.target.true.deriv, cov.target = cov.target))
        }
    }

    if (output.derivatives) {
        if (length(rho)>0) {
            return(cbind(data.table(which = paste0("rho = ", rho)), out.derivatives))
        } else {
            return(cbind(data.table(which = paste0("theta = ", theta)), out.derivatives))
        }
    }
    
    message("---------------------------------")
        
    print(paste0("se.target.crude = ", mean(se.target.crude <- sapply(est.calibration, function(est) est["target.se.crude"]))))
    print(paste0("cov.crude.se = ",  mean(est.target - 1.96*se.target.crude <= true.effect["truth_1"] &
                                          est.target + 1.96*se.target.crude >= true.effect["truth_1"])))

    message("---------------------------------")
    
    if (plot) {
        par(mfrow = c(1,2))
        hist(est.alpha); abline(v = as.numeric(true.effect["alpha"]))
        hist(est.target); abline(v = as.numeric(true.effect["truth_1"]))
    }

    if (length(alpha2)>0 | length(theta2)>0 | length(rho2)>0) {

        target2 <- "target"; misspecify <- c(); truncate.weights <- 0

        if (length(alpha2)>0) {
            get.output2 <- readRDS(file = paste0("./output/",
                                                 ifelse(setting == 1, "setting1-", ""),
                                                 "estimation-", target2, "-alpha", alpha2,
                                                 ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                                 ifelse("z" %in% misspecify, "-misspecify-z", ""),
                                                 ifelse("l" %in% misspecify, "-misspecify-l", ""),
                                                 ifelse(truncate.weights>0, paste0("-truncate-weights-", truncate.weights), ""),
                                                 "-M", M, "-N", N, 
                                                 ".rds"))
        } else {
            get.output2 <-  readRDS(file = paste0("./output/",
                                                  ifelse(setting == 1, "setting1-", ""),
                                                  "estimation-calibration-",
                                                  ifelse(length(rho2)>0,
                                                         paste0("rho-", rho2),
                                                         paste0("theta-", theta2)),
                                                  "-M", M, "-N", N, 
                                                  ".rds"))
        }

        get.output2 <- get.output2[1:min(length(est.calibration), length(get.output2))] 
        get.output <- est.calibration[1:min(length(est.calibration), length(get.output2))]

        print(length(get.output))

        if ("eic" %in% names(get.output2[[M1]])) {
            if (length(alpha2)>0) {
                if (eic) get.output.eic <- est.calibration.eic
                get.output2.eic <- lapply(get.output2, function(rep) rep[["eic"]])
                get.output2 <- lapply(get.output2, function(rep) rep[["estimate"]])
            } else {
                if (eic) get.output.eic <- est.calibration.eic
                get.output2.eic <- lapply(get.output2, function(rep) rep[grep("eic", names(rep))])
                get.output2 <- lapply(get.output2, function(rep) rep[["estimate"]])
            }
        } else {
            eic <- FALSE
        }

        truth <- as.numeric(true.effect["truth_1"])

        if (length(alpha2)>0) {
            truth2 <- readRDS(file = paste0("./output/",
                                            ifelse(setting == 1, "setting1-", ""),
                                            "true-estimands-alpha-", alpha2,
                                            ".rds"))[ifelse(target2 == "target", "truth_1", "truth_z")]
        } else {
            truth2 <- readRDS(file = paste0("./output/",
                                            ifelse(setting == 1, "setting1-", ""),
                                            "true-estimand-",
                                            ifelse(length(rho2)>0,
                                                   paste0("rho-", rho2),
                                                   paste0("theta-", theta2)),
                                            ".rds"))[ifelse(target2 == "target", "truth_1", "truth_z")]
        }

        est.output <- sapply(get.output, function(rep) rep["target.est"])
        se.output <- sapply(get.output, function(rep) rep["target.se"])

        if (length(alpha2)>0) {
            est.output2 <- sapply(get.output2, function(rep) rep["est"])
            se.output2 <- sapply(get.output2, function(rep) rep["se"])
        } else {
            est.output2 <- sapply(get.output2, function(rep) rep["target.est"])
            se.output2 <- sapply(get.output2, function(rep) rep["target.se"])
        }
        
        est.diff <- est.output-est.output2

        message("---------------------------------")        
        if (eic) {
            if (length(alpha2)>0) {
                se.diff <- sapply(1:length(get.output2.eic), function(rep.ii) {
                    sqrt(mean((get.output.eic[[rep.ii]]$eic - get.output2.eic[[rep.ii]])^2/N))
                })
            } else {
                se.diff <- sapply(1:length(get.output2.eic), function(rep.ii) {
                    sqrt(mean((get.output.eic[[rep.ii]]$eic - get.output2.eic[[rep.ii]]$eic)^2/N))
                })
            }
            print("se.diff based on eic for difference")
        } else {
            se.diff <- sqrt(se.output^2+se.output2^2)
            print("se.diff approximated via sum of variances")
        }
        message("---------------------------------")

        ci.lower <- est.diff-1.96*se.diff
        ci.upper <- est.diff+1.96*se.diff

        ci.oracle.lower <- est.diff-1.96*sd(est.diff)
        ci.oracle.upper <- est.diff+1.96*sd(est.diff)

        if (length(rho)>0) {
            if (length(rho2)>0) {
                if (sign.contrast == -1) {
                    comparison <- paste0("rho = ", rho2, " vs rho = ", rho)
                } else {
                    comparison <- paste0("rho = ", rho, " vs rho = ", rho2)
                }
            } else if (sign.contrast == -1) {
                comparison <- paste0("alpha = ", alpha2, " vs rho = ", rho)
            } else {
                comparison <- paste0("rho = ", rho, " vs alpha = ", alpha2)
            }
        } else {
            if (length(alpha2)>0) {
                if (sign.contrast == -1) {
                    comparison <- paste0("alpha = ", alpha2, " vs theta = ", theta)
                } else {
                    comparison <- paste0("theta = ", theta, " vs alpha = ", alpha2)
                }
            } else if (length(theta2)>0) {
                if (sign.contrast == -1) {
                    comparison <- paste0("theta = ", theta2, " vs theta = ", theta)
                } else {
                    comparison <- paste0("theta = ", theta, " vs theta = ", theta2)
                }
            } else if (length(rho2)>0) {
                if (sign.contrast == -1) {
                    comparison <- paste0("rho = ", rho2, " vs theta = ", theta)
                } else {
                    comparison <- paste0("theta = ", theta, " vs rho = ", rho2)
                }
            }
        }

        return(data.table(comparison = comparison,
                          est.diff = sign.contrast*mean(est.diff), truth = sign.contrast*(truth-truth2),
                          sd.diff = sd(est.diff), se.diff = mean(se.diff),
                          cov = mean(ci.lower <= truth-truth2 & ci.upper >= truth-truth2),
                          oracle.cov = mean(est.diff-1.96*sd(est.diff) <= truth-truth2 &
                                            est.diff+1.96*sd(est.diff) >= truth-truth2),
                          ci.lower = sign.contrast*mean(ci.lower), ci.upper = sign.contrast*mean(ci.upper),
                          power = mean((ci.lower>0 & ci.upper>0) |
                                       (ci.lower<0 & ci.upper<0))))
    }

    return(list(alpha = c(est = mean(est.alpha),
                          truth = as.numeric(true.effect["alpha"]), sd = sd(est.alpha),
                          se = mean(se.alpha), cov = cov.alpha, oracle.cov = oracle.cov.alpha), 
                target = c(est = mean(est.target), truth = as.numeric(true.effect["truth_1"]), sd = sd(est.target),
                           se = mean(se.target), cov = cov.target, oracle.cov = oracle.cov.target)))
}

    
######################################################################


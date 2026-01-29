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
library(devtools)
load_all('./simevent/')
library(simevent)

#-------------------------------------------------------------------------------------------#
## source relevant scripts
#-------------------------------------------------------------------------------------------#

source("./simulations/int.effect.alpha.R")
source("./R/alpha.scaling.estimation.fun.R")
source("./R/estimate.alpha.tmle.R")
source("./R/estimate.derivative.R")
source("./R/calibration.fun.R")

#-------------------------------------------------------------------------------------------#
## simulation study 
#-------------------------------------------------------------------------------------------#

M <- 500
N <- 1000

compute.true.derivatives <- compute.true <- FALSE

setting <- 1

if (setting == 0) {
    eta <- rep(0.1, 4)
} else if (setting == 1) {
    eta <- c(0.1, 0.025, 0.085, 0.1)
}

param.list <- list(
    c(theta = 0.25),
    c(theta = 1/3),
    c(theta = 0.5),
    c(theta = 2/3),
    c(theta = 0.75),
    c(rho = 0.6),
    c(rho = 0.8),
    c(rho = 1.3),
    c(rho = 1.5))

parallelize.Z <- 20

for (param in param.list) {

    theta <- param["theta"]
    rho <- param["rho"]

    if (is.na(rho)) {
        rho <- NULL
    }

    if (is.na(theta)) {
        theta <- NULL
    }
    
    est.list <- list()

    for (m in 1:M) {
    
        set.seed(m+33339)
    
        dt <- simTreatment(N = N, 
                           eta = eta, 
                           nu = rep(1.1, 4),
                           beta_L_A = 3, beta_L_D = 2.5,
                           beta_A_D = -0.5, beta_A_L = -2.5,
                           beta_L0_A = 1,
                           lower = 10^(-150),      
                           upper = 1e3
                           )

        setnames(dt, "ID", "id")
        setnames(dt, "Time", "time")
        setnames(dt, "Delta", "delta")

        est.list[[m]] <- calibration.fun(theta = theta, rho = rho, tau = 3,
                                         parallelize.Z = parallelize.Z,
                                         output.eic = TRUE,
                                         test.derivatives = FALSE)

        saveRDS(est.list,
                file = paste0("./output/",
                              ifelse(setting == 1, "setting1-", ""),
                              "estimation-calibration-",
                              ifelse(length(rho)>0,
                                     paste0("rho-", rho),
                                     paste0("theta-", theta)),
                              "-M", M, "-N", N, 
                              ".rds"))
    }
}


#-------------------------------------------------------------------------------------------#
## compute true values 
#-------------------------------------------------------------------------------------------#

if (compute.true) {
    
    for (param in param.list) {

        theta <- param["theta"]
        rho <- param["rho"]

        if (is.na(rho)) {
            rho <- NULL
        }

        if (is.na(theta)) {
            theta <- NULL
        }
    
        N2 <- 1e5

        effect_alpha <- function(alpha, tau, inner.rep = 25) {
            int_effect <- lapply(1:inner.rep, function(rep) {
                tmp.out <- try(int_effect_alpha(N = N2,
                                                tau = tau,
                                                alpha = alpha,
                                                eta = eta, 
                                                nu = rep(1.1, 4),
                                                beta_L_A = 3, beta_L_D = 2.5,
                                                beta_A_D = -0.5, beta_A_L = -2.5,
                                                beta_L0_A = 1 
                                                ))
                if (any(class(tmp.out) == "try-error")) {
                    print(paste0("alpha = ", alpha, ", rep = ", rep))
                    return(list(effect_z = NA,
                                effect_1 = NA))
                }
                return(tmp.out)
            })
            return(c(truth_1 = mean(sapply(int_effect, function(rep) rep$effect_1), na.rm = TRUE),
                     truth_z = mean(sapply(int_effect, function(rep) rep$effect_z), na.rm = TRUE)))
        }

        if (length(rho)>0) {
            effect_1 <- effect_alpha(1, tau = 3)
            true.alpha.calibrated <-
                estimate.alpha.tmle(NULL,
                                    theta = rho*as.numeric(effect_1["truth_z"]),
                                    effect_alpha = effect_alpha,
                                    c_n = 0.00001, max_iter = 300)
            true.effect <- c(alpha = true.alpha.calibrated$alpha.hat,
                             effect_alpha(true.alpha.calibrated$alpha.hat, tau = 3))
        } else {
            true.alpha.calibrated <-
                estimate.alpha.tmle(NULL, theta = theta,
                                    effect_alpha = effect_alpha,
                                    c_n = 0.00001, max_iter = 100)
            true.effect <- c(alpha = true.alpha.calibrated$alpha.hat,
                             effect_alpha(true.alpha.calibrated$alpha.hat, tau = 3))
        }

        saveRDS(true.effect,
                file = paste0("./output/",
                              ifelse(setting == 1, "setting1-", ""),
                              "true-estimand-",
                              ifelse(length(rho)>0,
                                     paste0("rho-", rho),
                                     paste0("theta-", theta)),
                              ".rds"))
    }
}

#-------------------------------------------------------------------------------------------#
## compute true values of derivatives 
#-------------------------------------------------------------------------------------------#

if (compute.true.derivatives) {

    for (param in param.list) {

        theta <- param["theta"]
        rho <- param["rho"]

        if (is.na(rho)) {
            rho <- NULL
        }

        if (is.na(theta)) {
            theta <- NULL
        }
    
        N2 <- 1e5

        h <- N2^{-1/6}

        effect_alpha <- function(alpha, tau, inner.rep = 25) {
            int_effect <- lapply(1:inner.rep, function(rep) {
                tmp.out <- try(int_effect_alpha(N = N2,
                                                tau = tau,
                                                alpha = alpha,
                                                eta = eta, 
                                                nu = rep(1.1, 4),
                                                beta_L_A = 3, beta_L_D = 2.5,
                                                beta_A_D = -0.5, beta_A_L = -2.5,
                                                beta_L0_A = 1 
                                                ))
                if (any(class(tmp.out) == "try-error")) {
                    print(paste0("alpha = ", alpha, ", rep = ", rep))
                    return(list(effect_z = NA,
                                effect_1 = NA))
                }
                return(tmp.out)
            })
            return(c(truth_1 = mean(sapply(int_effect, function(rep) rep$effect_1), na.rm = TRUE),
                     truth_z = mean(sapply(int_effect, function(rep) rep$effect_z), na.rm = TRUE)))
        }

        true.effect <- readRDS(file = paste0("./output/",
                                             ifelse(setting == 1, "setting1-", ""),
                                             "true-estimand-",
                                             ifelse(length(rho)>0,
                                                    paste0("rho-", rho),
                                                    paste0("theta-", theta)),
                                             ".rds"))

        true.alpha <- true.effect["alpha"]

        true.alpha.plus <- effect_alpha(true.alpha+h, tau = 3)

        if (true.alpha > h) {
            true.alpha.minus <- effect_alpha(true.alpha-h, tau = 3)
            true.derivative.z <- as.numeric((true.alpha.plus["truth_z"] - true.alpha.minus["truth_z"]) / (2*h))
            true.derivative.1 <- as.numeric((true.alpha.plus["truth_1"] - true.alpha.minus["truth_1"]) / (2*h))
        } else {
            true.derivative.z <- as.numeric((true.alpha.plus["truth_z"] - true.effect["truth_z"]) / h)
            true.derivative.1 <- as.numeric((true.alpha.plus["truth_1"] - true.effect["truth_1"]) / h)
        }

        print(true.derivatives <- c(true.derivative.z = true.derivative.z, true.derivative.1 = true.derivative.1))
        
        saveRDS(true.derivatives,
                file = paste0("./output/",
                              ifelse(setting == 1, "setting1-", ""),
                              "true-derivatives-",
                              ifelse(length(rho)>0,
                                     paste0("rho-", rho),
                                     paste0("theta-", theta)),
                              ".rds"))
        
    }
}


######################################################################


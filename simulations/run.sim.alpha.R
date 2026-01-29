### run.sim.intensity.intervention.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Dec 18 2025 (08:58) 
## Version: 
## Last-Updated: Jan 29 2026 (08:40) 
##           By: Helene
##     Update #: 256
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
library(devtools)
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

which.alphas <- sort(c(0, 0.1, 0.25, 2.25, 3, 0.5, 1, 1.5))

target.parameters <- c("target", "auxiliary")

parallelize.Z <- 30

est.lists <- lapply(which.alphas, function(alpha) list())
est.aux.lists <- lapply(which.alphas, function(alpha) list())

names(est.lists) <- which.alphas
names(est.aux.lists) <- which.alphas

setting <- 1

if (setting == 0) {
    eta <- rep(0.1, 4) 
} else if (setting == 1) {
    eta <- c(0.1, 0.025, 0.085, 0.1)
}

for (misspecify in list(c(), c("1", "l"))) {

    if ("1" %in% misspecify) {
        model.1 <- "Surv(tstart, tstop, delta == 1)~L0+A.1"
    } else {
        model.1 <- "Surv(tstart, tstop, delta == 1)~L0+A.1+L.1"
    }

    if ("z" %in% misspecify) {
        model.z <- "Surv(tstart, tstop, delta == 2)~L0"
    } else {
        model.z <- "Surv(tstart, tstop, delta == 2)~L0+L.1"
    }

    if ("l" %in% misspecify) {
        model.l <- "Surv(tstart, tstop, delta == 3)~L0"
    } else {
        model.l <- "Surv(tstart, tstop, delta == 3)~L0+A.1"
    }
   
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

        for (alpha in which.alphas) {

            if ("target" %in% target.parameters) {
                est.lists[[as.character(alpha)]][[m]] <- 
                    alpha.est.fun(dt, tau = 3, alpha = alpha,
                                  fit.type.1 = list(model = model.1,
                                                    fit = "cox"),
                                  fit.type.z = list(model = model.z,
                                                    fit = "cox",
                                                    atrisk = function(dt) (dt[["A.1"]] == 0)),
                                  fit.type.l = list(model = model.l,
                                                    fit = "cox",
                                                    atrisk = function(dt) (dt[["L.1"]] == 0)),
                                  output.eic = TRUE,
                                  parallelize.Z = parallelize.Z, verbose = TRUE)
            
                saveRDS(est.lists[[as.character(alpha)]],
                        file = paste0("./output/",
                                      ifelse(setting == 1, "setting1-", ""),
                                      "estimation-target-alpha", alpha,
                                      ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                      ifelse("z" %in% misspecify, "-misspecify-z", ""),
                                      ifelse("l" %in% misspecify, "-misspecify-l", ""),
                                      "-M", M, "-N", N, 
                                      ".rds"))
            }

            if ("auxiliary" %in% target.parameters) {
                est.aux.lists[[as.character(alpha)]][[m]] <- 
                    alpha.est.fun(dt, tau = 3, alpha = alpha, parameter = "auxiliary",
                                  fit.type.1 = list(model = model.1,
                                                    fit = "cox"),
                                  fit.type.z = list(model = model.z,
                                                    fit = "cox",
                                                    atrisk = function(dt) (dt[["A.1"]] == 0)),
                                  fit.type.l = list(model = model.l,
                                                    fit = "cox",
                                                    atrisk = function(dt) (dt[["L.1"]] == 0)),
                                  output.eic = TRUE,
                                  parallelize.Z = parallelize.Z, verbose = TRUE)

                saveRDS(est.aux.lists[[as.character(alpha)]],
                        file = paste0("./output/",
                                      ifelse(setting == 1, "setting1-", ""),
                                      "estimation-auxiliary-alpha", alpha,
                                      ifelse("1" %in% misspecify, "-misspecify-1", ""),
                                      ifelse("z" %in% misspecify, "-misspecify-z", ""),
                                      ifelse("l" %in% misspecify, "-misspecify-l", ""),
                                      "-M", M, "-N", N, 
                                      ".rds"))
            }
        }
    }
}

#-------------------------------------------------------------------------------------------#
## compute true values 
#-------------------------------------------------------------------------------------------#

N2 <- ifelse(local, 1e4, 1e5)

get.curve <- FALSE

if (!get.curve) {

    for (alpha in which.alphas) {

        int_effect <- lapply(1:50, function(rep) {
            tmp.out <- try(int_effect_alpha(N = N2,
                                            tau = 3,
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

        print(truth.out <- c(truth_1 = mean(sapply(int_effect, function(rep) rep$effect_1), na.rm = TRUE),
                             truth_z = mean(sapply(int_effect, function(rep) rep$effect_z), na.rm = TRUE)))

        saveRDS(truth.out,
                file = paste0("./output/",
                              ifelse(setting == 1, "setting1-", ""),
                              "true-estimands-alpha-", alpha,
                              ".rds"))
    }
} else {

    effect_alpha <- lapply(alpha_seq <- seq(0, 3, length = 100), function(alpha) {

        int_effect <- lapply(1:10, function(rep) int_effect_alpha(N = N2,
                                                                  alpha = alpha,
                                                                  eta = eta, 
                                                                  nu = rep(1.1, 4),
                                                                  beta_L_A = 3, beta_L_D = 2.5,
                                                                  beta_A_D = -0.5, beta_A_L = -2.5,
                                                                  beta_L0_A = 1 
                                                                  ))

        print(truth.out <- c(truth_1 = mean(sapply(int_effect, function(rep) rep$effect_1)),
                             truth_z = mean(sapply(int_effect, function(rep) rep$effect_z))))

        return(truth.out)

    })

    names(effect_alpha) <- alpha_seq

    saveRDS(truth.out,
            file = paste0("./output/",
                          ifelse(setting == 1, "setting1-", ""),
                          "true-estimands-alpha-", min(alpha_seq), "-", max(alpha_seq),
                          ".rds"))

}

######################################################################


### int.effect.alpha.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Dec 19 2025 (10:47) 
## Version: 
## Last-Updated: Jan 19 2026 (09:37) 
##           By: Helene
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

int_effect_alpha <- function(N2 = 1e4,
                             alpha = 0.5,
                             plot = FALSE,
                             beta_L_A = 3, beta_L_D = 2.5,
                             beta_A_D = -0.5, beta_A_L = -2.5,
                             beta_L0_A = 1, #???? beta_L0_L, beta_L0_D,
                             nu = rep(1.1, 4),
                             eta = c(rep(0.1, 4)[1], 0.025, 0.085, rep(0.1, 4)[4]),
                             tau = 3) {

    # Generate large data
    data_G1 <- simTreatment(N = N2,
                            cens = 0, 
                            eta = c(eta[1:2],eta[3]*alpha, eta[4]), 
                            nu = nu,
                            beta_L_A = beta_L_A, beta_L_D = beta_L_D,
                            beta_A_D = beta_A_D, beta_A_L = beta_A_L,
                            beta_L0_A = beta_L0_A,
                            lower = 10^(-150),      
                            upper = 1e3)
   
    if (plot) grid.arrange(plotEventData(data_G1[1:250]), nrow = 1)
    
    #Proportion of subjects dying before some time $\tau$ in treatment group
    prop_G1 <- data_G1[Delta == 1, mean(Delta == 1 & Time < tau)] # with intervention
    prop_G1_A <- mean(data_G1[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]]) # with intervention

    return(list(effect_z = prop_G1_A,
                effect_1 = prop_G1))
}

######################################################################
### int.effect.alpha.R ends here

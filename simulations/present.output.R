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

source("./simulations/get.output.R")

out.N1000.target <-
    rbind(
        get.output.fun(alpha = 0, target = "target", N = 1000, misspecify = c(), setting = 1),
        get.output.fun(alpha = 0.1, target = "target", N = 1000, misspecify = c(), setting = 1),
        get.output.fun(alpha = 0.25, target = "target", N = 1000, misspecify = c(), setting = 1),
        get.output.fun(alpha = 0.5, target = "target", N = 1000, misspecify = c(), setting = 1),
        get.output.fun(alpha = 1, target = "target", N = 1000, misspecify = c(), setting = 1),
        get.output.fun(alpha = 1.5, target = "target", N = 1000, misspecify = c(), setting = 1),
        get.output.fun(alpha = 2.25, target = "target", N = 1000, misspecify = c(), setting = 1),
        get.output.fun(alpha = 3, target = "target", N = 1000, misspecify = c(), setting = 1), fill = TRUE)

out.N1000.aux <-
    rbind(get.output.fun(alpha = 0, target = "auxiliary", N = 1000, misspecify = c(), setting = 1),
          get.output.fun(alpha = 0.1, target = "auxiliary", N = 1000, misspecify = c(), setting = 1),
          get.output.fun(alpha = 0.25, target = "auxiliary", N = 1000, misspecify = c(), setting = 1),
          get.output.fun(alpha = 0.5, target = "auxiliary", N = 1000, misspecify = c(), setting = 1),
          get.output.fun(alpha = 1, target = "auxiliary", N = 1000, misspecify = c(), setting = 1),
          get.output.fun(alpha = 1.5, target = "auxiliary", N = 1000, misspecify = c(), setting = 1),
          get.output.fun(alpha = 2.25, target = "auxiliary", N = 1000, misspecify = c(), setting = 1),
          get.output.fun(alpha = 3, target = "auxiliary", N = 1000, misspecify = c(), setting = 1), fill = TRUE)

out.N1000.target.miss <-
    rbind(get.output.fun(alpha = 0, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 0.1, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 0.25, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 0.5, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 1, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 1.5, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 2.25, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 3, target = "target", N = 1000, setting = 1, misspecify = c("1", "l")))

out.N1000.aux.miss <-
    rbind(get.output.fun(alpha = 0, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 0.1, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 0.25, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 0.5, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 1, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 1.5, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 2.25, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")),
          get.output.fun(alpha = 3, target = "auxiliary", N = 1000, setting = 1, misspecify = c("1", "l")))

pdat.N1000 <- rbind(
    out.N1000.target[, parameter := "target parameter"][, specification := "correctly specified nuisance"],
    out.N1000.aux[, parameter := "auxiliary parameter"][, specification := "correctly specified nuisance"],
    out.N1000.target.miss[, parameter := "target parameter"][, specification := "misspecified nuisance"],
    out.N1000.aux.miss[, parameter := "auxiliary parameter"][, specification := "misspecified nuisance"])

pdat.N1000[, cov.label := paste0("cov = ", round(cov*100,1), "%")]
pdat.N1000[alpha == 0 & parameter == "auxiliary parameter", cov.label := ""]

(p.N1000 <- ggplot(pdat.N1000) +
     theme_bw() +
     geom_point(aes(x = alpha, y = est)) +
     geom_line(aes(x = alpha, y = est)) +
     geom_point(aes(x = alpha, y = truth), alpha = 0.3) +
     geom_line(aes(x = alpha, y = g.est), linetype = "dashed", col = "red", alpha = 0.3) +
     geom_line(aes(x = alpha, y = truth), linetype = "dashed", alpha = 0.3) +
     geom_segment(aes(x = alpha, xend = alpha,
                      y = est-1.96*sd, yend = est+1.96*sd)) +
     geom_segment(aes(x = alpha-0.015, xend = alpha+0.015,
                      y = est-1.96*sd, yend = est-1.96*sd)) +
     geom_segment(aes(x = alpha-0.015, xend = alpha+0.015,
                      y = est+1.96*sd, yend = est+1.96*sd)) +
     geom_text(aes(x = alpha, y = est+1.96*sd, label = cov.label), vjust = -0.5, size = 5) +
     facet_grid(parameter~specification, scales = "free_y") +
     xlab(expression("value of "*alpha)) + ylab("estimate of intervention-specific parameter") + 
     theme(strip.text.x = element_text(size = 15),
           strip.text.y = element_text(size = 15),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 12), 
           strip.background = element_rect(fill = "white"),
           plot.title = element_text(size = 15, hjust = 0.5)) +
     labs(caption = "sample size n = 1000, M = 500 simulation repetitions") + 
     ggtitle(expression("Estimation of " * alpha * "-specific auxiliary and target parameters")))

ggsave("./figures/fig-sim-study-results-N1000-M500-alpha-all.pdf",
       width = 12, height = 9)


pdat.weights <- do.call("rbind", lapply(sort(c(0, 0.1, 0.25, 2.25, 3, 0.5, 1, 1.5)), function(alpha) {
    out.alpha <- get.output.fun(alpha = alpha, target = "target", N = 1000, misspecify = c(), setting = 1, output.weights = TRUE)
    out.alpha[, alpha := alpha]
    return(out.alpha)
}))

pdat.weights.long <- melt(pdat.weights, id.var = "alpha")

pdat.weights.long[, unique(variable)]

(p.weights.histogram <-
     ggplot(pdat.weights.long[variable %in% c("q100.weight", "q99.weight", "q975.weight")]) +
     theme_bw() +
     scale_y_log10() +
     ylab("log-scaled value of maximal weight") + 
     geom_boxplot(aes(group = alpha, y = value)) +
     theme(strip.text.x = element_text(size = 15),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 12), 
           strip.background = element_rect(fill = "white"),
           plot.title = element_text(size = 15, hjust = 0.5)) + 
     facet_wrap(.~variable, scales = "free_y", nrow = 1))

ggsave("./figures/fig-sim-study-results-N1000-M500-alpha-weights.pdf",
       width = 10, height = 4)
       

results.theta.N1000 <- do.call("rbind", lapply(c(0.5, 1/3, 0.25, 2/3, 0.75), function(theta) {
    out <- get.calibration.output.fun(N = 1000, theta = theta, setting = 1, plot = TRUE)
    return(data.table(rbind(c(theta = theta, which = "alpha", round(out$alpha, 4)),
                            c(theta = theta, which = "target", round(out$target, 4)))))
}))[, sd := as.numeric(sd)][, est := as.numeric(est)][, oracle.cov := as.numeric(oracle.cov)][, cov := as.numeric(cov)][, theta := as.numeric(theta)]

results.theta.N1000[which == "alpha", which := "'estimation of '*alpha^{theta}*'(P)'"]
results.theta.N1000[which == "target", which := "'estimation of composite parameter '*Psi['1']^{theta}*'(P)'"]

(p.theta.N1000 <- ggplot(results.theta.N1000) +
     geom_point(aes(x = theta, y = est)) +
     geom_point(aes(x = theta, y = as.numeric(truth)), alpha = 0.3## , shape = 4, size = 4
                ) +
     facet_wrap(.~which, scales = "free", labeller = label_parsed) +
     geom_segment(aes(x = theta, xend = theta,
                      y = est-1.96*sd, yend = est+1.96*sd)) +
     geom_segment(aes(x = theta-0.005, xend = theta+0.005,
                      y = est-1.96*sd, yend = est-1.96*sd)) +
     geom_segment(aes(x = theta-0.005, xend = theta+0.005,
                      y = est+1.96*sd, yend = est+1.96*sd)) +
     geom_text(aes(x = theta, y = est+1.96*sd, label = paste0("cov = ", cov*100, "%")), vjust = -0.5) +
     geom_text(aes(x = theta, y = est-1.96*sd, label = paste0("oracle = ", oracle.cov*100, "%")), vjust = 1.75) + 
     theme_bw() +
     xlab(expression("value of "*theta)) + ylab("estimate and true value") + 
     labs(caption = "sample size n = 1000, M = 500 simulation repetitions") + 
     ggtitle(expression("estimation of calibrated and composite parameter reaching target level "*theta)) + 
     theme(strip.text.x = element_text(size = 15),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 12), 
           strip.background = element_rect(fill = "white"),
           plot.title = element_text(size = 15, hjust = 0.5)))

ggsave("./figures/fig-sim-study-results-N1000-M500-theta.pdf",
       width = 10, height = 6)


results.rho.N1000 <- do.call("rbind", lapply(c(0.6, 0.8, 1.3, 1.5), function(rho) {
    out <- get.calibration.output.fun(N = 1000, rho = rho, setting = 1, plot = TRUE)
    return(data.table(rbind(c(rho = rho, which = "alpha", round(out$alpha, 4)),
                            c(rho = rho, which = "target", round(out$target, 4)))))
}))[, sd := as.numeric(sd)][, est := as.numeric(est)][, oracle.cov := as.numeric(oracle.cov)][, cov := as.numeric(cov)][, rho := as.numeric(rho)]

results.rho.N1000[which == "alpha", which := "'estimation of '*alpha^{rho}*'(P)'"]
results.rho.N1000[which == "target", which := "'estimation of composite parameter '*Psi['1']^{rho}*'(P)'"]

(p.rho.N1000 <- ggplot(results.rho.N1000) +
     geom_point(aes(x = rho, y = est)) +
     geom_point(aes(x = rho, y = as.numeric(truth)), alpha = 0.3## , shape = 4, size = 4
                ) +
     facet_wrap(.~which, scales = "free", labeller = label_parsed) +
     geom_segment(aes(x = rho, xend = rho,
                      y = est-1.96*sd, yend = est+1.96*sd)) +
     geom_segment(aes(x = rho-0.005, xend = rho+0.005,
                      y = est-1.96*sd, yend = est-1.96*sd)) +
     geom_segment(aes(x = rho-0.005, xend = rho+0.005,
                      y = est+1.96*sd, yend = est+1.96*sd)) +
     geom_text(aes(x = rho, y = est+1.96*sd, label = paste0("cov = ", cov*100, "%")), vjust = -0.5) +
     geom_text(aes(x = rho, y = est-1.96*sd, label = paste0("oracle = ", oracle.cov*100, "%")), vjust = 1.75) + 
     theme_bw() +
     xlab(expression("value of "*rho)) + ylab("estimate and true value") + 
     labs(caption = "sample size n = 1000, M = 500 simulation repetitions") + 
     ggtitle(expression("estimation of calibrated and composite parameter reaching target level "*rho)) + 
     theme(strip.text.x = element_text(size = 15),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 12), 
           strip.background = element_rect(fill = "white"),
           plot.title = element_text(size = 15, hjust = 0.5)))

ggsave("./figures/fig-sim-study-results-N1000-M500-rho.pdf",
       width = 10, height = 6)

pdat.test <- rbind(
    get.output.fun(alpha = 0, target = "target", N = 1000, misspecify = c(), setting = 1,
                   alpha2 = 1),
    get.calibration.output.fun(N = 1000, theta = 0.25, setting = 1, plot = TRUE, theta2 = 0.75),
    get.calibration.output.fun(N = 1000, rho = 0.6, setting = 1, plot = TRUE, alpha2 = 1),
    get.calibration.output.fun(N = 1000, rho = 1.5, setting = 1, plot = TRUE, alpha2 = 1, sign.contrast = -1))

pdat.test[, cov.label := paste0("cov = ", round(cov*100,1), "%")]
pdat.test[, power.label := paste0("power = ", round(power*100,1), "%")]

pdat.test[, comparison2 := c("a", "b", "c", "d")]
pdat.test[, comparison2 := paste0(comparison2, " ", comparison)]

(p.test <- ggplot(pdat.test) +
     theme_bw() +
     geom_hline(aes(yintercept = 0), color = "red", alpha = 0.3, linetype = "dashed") + 
     geom_point(aes(x = comparison2, y = est.diff)) +
     geom_point(aes(x = comparison2, y = truth), alpha = 0.3) +
     geom_segment(aes(x = comparison2, xend = comparison2,
                      y = est.diff-1.96*sd.diff, yend = est.diff+1.96*sd.diff)) +
     geom_text(aes(x = comparison2, y = est.diff+1.96*sd.diff, label = cov.label), vjust = -0.5) +
     geom_text(aes(x = comparison2, y = est.diff-1.96*sd.diff, label = power.label), vjust = 1.25) +
     xlab("") + ylab("estimate and true value of contrasts") +
     scale_x_discrete(labels = c("a alpha = 0 vs alpha = 1" = expression(alpha*"=0 vs "*alpha*"=1"),
                                 "d alpha = 1 vs rho = 1.5" = expression(alpha*"=1 vs "*rho*"=1.5"),
                                 "c rho = 0.6 vs alpha = 1" = expression(rho*"=0.6 vs "*alpha*"=1"),
                                 "b theta = 0.25 vs theta = 0.75" = expression(theta*"=0.25 vs "*theta*"=0.75"))) + 
     theme(strip.text.x = element_text(size = 15),
           axis.title.x = element_text(size = 15),
           axis.title.y = element_text(size = 12),
           axis.text.x = element_text(size = 14),
           strip.background = element_rect(fill = "white"),
           plot.title = element_text(size = 15, hjust = 0.5)) +
     labs(caption = "sample size n = 1000, M = 500 simulation repetitions") + 
     ggtitle(expression("Estimation of contrasts")))

ggsave("./figures/fig-sim-study-results-N1000-M500-contrasts.pdf",
       width = 8, height = 5) 



######################################################################


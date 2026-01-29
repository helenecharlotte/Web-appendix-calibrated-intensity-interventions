### Code:

alpha.est.fun <- function(dt,
                          alpha = 1,
                          intervention.A = 1,
                          tau = 1.2,
                          parameter = "target",
                          fit.type.1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+A.1+L.1",
                                            fit = "cox"),
                          fit.type.z = list(model = "Surv(tstart, tstop, delta == 2)~L0+L.1",
                                            fit = "cox",
                                            atrisk = function(dt) (dt[["A.1"]] == 0)),
                          fit.type.l = list(model = "Surv(tstart, tstop, delta == 3)~L0+A.1",
                                            fit = "cox",
                                            atrisk = function(dt) (dt[["L.1"]] == 0)),
                          fit.type.0 = list(model = "Surv(tstart, tstop, delta == 0)~L0+A.1+L.1",
                                            fit = "cox"),
                          fit.treatment = NULL,
                          output.eic = FALSE,
                          target = TRUE,
                          truncate.weights = 0,
                          use.exponential = FALSE,
                          verbose = FALSE, verbose2 = FALSE,
                          browse.hal = FALSE,
                          browse = FALSE,
                          browse2 = FALSE,
                          browse3 = FALSE,
                          #--- simpler estimators: 
                          baseline.tmle = FALSE,
                          naive.gcomp = FALSE, 
                          standard.np = FALSE,
                          max.iter = 100,
                          #--- HAL parameters:
                          lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                          min.no.of.ones = 0.01,
                          cut.one.way = 15,
                          cut.time.varying = 5,
                          cut.time = 35,
                          cut.two.way = 5,
                          penalize.time = FALSE,
                          reduce.seed.dependence = FALSE,
                          event.dependent.cv = FALSE,
                          cv.hal.fit = FALSE,
                          parallelize.Z = 1,
                          parallelize.cve = parallelize.Z,
                          parallelize.predict = parallelize.Z,
                          tmle.criterion = 1,
                          V = 10, 
                          cv.glmnet = FALSE,
                          seed.hal = NULL,
                          output.lambda.cvs = FALSE
                          ) {

    dt <- copy(dt)

    n <- length(unique(dt[["id"]]))

    #-- if no "fit" is specified, then specify it as cox: 
    if (!is.list(fit.type.1)) {
        fit.type.1 <- list(model = fit.type.1, fit = "cox")
    }

    if (!is.list(fit.type.z)) {
        fit.type.z <- list(model = fit.type.z, fit = "cox")
    }

    if (!is.list(fit.type.l)) {
        fit.type.l <- list(model = fit.type.l, fit = "cox")
    }

    if (!is.list(fit.type.0)) {
        fit.type.0 <- list(model = fit.type.0, fit = "cox")
    }

    #-- if no "fit" is specified for treatment model, then specify it as glm: 
    if (!is.list(fit.treatment)) {
        fit.treatment <- list(model = fit.treatment, fit = "glm")
    }
    
    #-- covariates/predictors extracted from models: 
    varnames <- unique(unlist(lapply(list(fit.type.1[["model"]], fit.type.z[["model"]], fit.type.l[["model"]], fit.type.0[["model"]], fit.treatment[["model"]]), function(fit.type) {
        if (length(fit.type)>0) return(strsplit(strsplit(fit.type, "~")[[1]][2], "\\+|\\*")[[1]])
    })))

    if (length(interaction.names <- grep(":", varnames, value = TRUE))>0) {
        varnames <- unique(c(varnames[!varnames %in% interaction.names],
                             str_split(interaction.names, ":")[[1]]))
    }

    #-- name of treatment variable:
    if (length(fit.treatment[["model"]])>0) {
        Aname <- strsplit(fit.treatment[["model"]], "~")[[1]][1]
    } else {
        Aname <- NULL 
    }
  
    if (typeof(dt[["time"]]) != "double") {
        warning("NB: the time variable is not numeric - will be converted")
        dt[, time := as.numeric(time)]
    }

    dt[, tstart := c(0, time[-.N]), by = "id"]
    dt[, tstop := time]

    dt[, A.1 := c(0, A[-.N]), by = "id"]
    dt[, L.1 := c(0, L[-.N]), by = "id"]

    #--------------------------------
    #-- "G-part"; for clever weight estimation:
    fit.cox.0 <- coxph(as.formula(fit.type.0[["model"]]), data = dt, 
                       control = coxph.control(timefix = FALSE))
    if (verbose) print(fit.cox.0)
    dt[, idN := 1:.N, by = "id"]
    if (fit.treatment[["fit"]] == "glm" & length(fit.treatment[["model"]])>0) {
        fitA <- glm(as.formula(fit.treatment[["model"]]), data=dt[idN == 1], family=binomial)
        if (verbose) print(summary(fitA))
    } else if (length(fit.treatment[["model"]])>0)  {
        print("NB: need to incorporate other estimations methods than glm for treatment")
    }
    if (length(fit.treatment[["model"]])>0) {
        dt[, probA := predict(fitA, newdata = dt, type = "response")]
    }
   
    #--------------------------------
    #-- "Q-part":
    fit.cox.1 <- coxph(as.formula(fit.type.1[["model"]]), data = dt, 
                       control = coxph.control(timefix = FALSE))
    fit.cox.z <- coxph(as.formula(fit.type.z[["model"]]), data = dt[fit.type.z[["atrisk"]](dt)], 
                       control = coxph.control(timefix = FALSE))
    fit.cox.l <- coxph(as.formula(fit.type.l[["model"]]), data = dt[fit.type.l[["atrisk"]](dt)], 
                       control = coxph.control(timefix = FALSE))
    if (verbose) print(fit.cox.1)
    if (verbose) print(fit.cox.z)
    if (verbose) print(fit.cox.l)

    dt2 <- copy(dt)

    #-- get baseline intensities: 
    tmp.type.1 <- suppressWarnings(setDT(basehaz(fit.cox.1, centered=TRUE)))[, dhazard.1 := c(hazard[1],diff(hazard))][, hazard.1 := hazard][, -"hazard", with = FALSE]
    tmp.type.z <- suppressWarnings(setDT(basehaz(fit.cox.z, centered=TRUE)))[, dhazard.z := c(hazard[1],diff(hazard))][, hazard.z := hazard][, -"hazard", with = FALSE]
    tmp.type.l <- suppressWarnings(setDT(basehaz(fit.cox.l, centered=TRUE)))[, dhazard.l := c(hazard[1],diff(hazard))][, hazard.l := hazard][, -"hazard", with = FALSE]
    tmp.type.0 <- suppressWarnings(setDT(basehaz(fit.cox.0, centered=TRUE)))[, dhazard.0 := c(hazard[1],diff(hazard))][, hazard.0 := hazard][, -"hazard", with = FALSE]

    #-- get all unique times; 
    unique.times <- sort(unique(dt2[["time"]]))
    unique.times <- unique.times[unique.times <= tau]
    all.times <- data.table(expand.grid(time = unique.times,
                                        id = unique(dt[["id"]])))

    #-- collect data with all time-points;
    tmp.inner <- merge(dt2[, time.obs := time], all.times, by = c("id", "time"), all = TRUE)[order(id, time.obs)]

    for (varname in varnames) {
        tmp.inner[, (varname) := na.locf(get(varname)), by = "id"]
    }

    #--------------------------------
    #-- expanded dataset to work with for TMLE:
    tmp3 <- merge(merge(merge(merge(tmp.inner[order(id, time)][is.na(delta), delta := 0][, -c("tstop"), with = FALSE],
                                    tmp.type.1[dhazard.1>0], by = "time", all = TRUE),
                              tmp.type.l[dhazard.l>0], by = "time", all = TRUE),
                        tmp.type.z[dhazard.z>0], by = "time", all = TRUE),
                  tmp.type.0[dhazard.0>0], by = "time", all = TRUE)[order(id,time)][!is.na(id)]

    #-- get last observed time for each line: 
    tmp3[delta == 1, time.obs2 := time.obs]
    tmp3[, time.obs.last := nafill(time.obs2, "locf"), by = "id"]
    tmp3[is.na(time.obs.last), time.obs.last := 0]
    
    tmp3[, time.obs := nafill(time.obs, "nocb"), by = "id"]
    tmp3[is.na(time.obs), time.obs := -Inf]

    tmp3[is.na(dhazard.0), dhazard.0 := 0]
    tmp3[is.na(dhazard.1), dhazard.1 := 0]
    tmp3[is.na(dhazard.z), dhazard.z := 0]
    tmp3[is.na(dhazard.l), dhazard.l := 0]

    tmp3[, A := cumsum(1*(delta == 2)), by = "id"]
    tmp3[, L := cumsum(1*(delta == 3)), by = "id"]
    tmp3[, A.1 := c(0, A[-.N]), by = "id"]
    tmp3[, L.1 := c(0, L[-.N]), by = "id"]
   
    #-- remame treatment variable to "Aobs" in tmp3 data:
    if (length(Aname)>0) {
        setnames(tmp3, Aname, "Aobs")
        #-- the treatment variable is then set to the interventional level: 
        tmp3[, (Aname) := intervention.A]
    }

    #--------------------------------
    #-- compute needed quantities for (initial) estimation and targeting

    tmp3[, at.risk.z := fit.type.z[["atrisk"]](tmp3)]
    tmp3[, at.risk.l := fit.type.l[["atrisk"]](tmp3)]
    
    tmp3[, exp.1 := exp(predict(fit.cox.1, newdata=tmp3, type="lp"))]
    tmp3[, exp.z := exp(predict(fit.cox.z, newdata=tmp3, type="lp"))]
    tmp3[, exp.l := exp(predict(fit.cox.l, newdata=tmp3, type="lp"))]

    tmp3[, P.1 := dhazard.1*exp.1]
    tmp3[, P.z := dhazard.z*exp.z]
    tmp3[, P.l := dhazard.l*exp.l]

    tmp3[, exp.0 := exp(predict(fit.cox.0, newdata=tmp3, type="lp"))]
    tmp3[, surv.0 := exp(-cumsum(exp.0*dhazard.0)), by = "id"]
    tmp3[, surv.0.1 := c(1, surv.0[-.N]), by = "id"]
    
    tmp3[, cum.hazard.z := cumsum(at.risk.z*P.z), by = "id"]
    tmp3[, cum.hazard.z.1 := c(0, cum.hazard.z[-.N]), by = "id"]

    #--------------------------------    
    #-- prediction part is over
    
    tmp3 <- tmp3[time <= tau]
    
    #--------------------------------    
    #-- to handle dependence on jump in the past:

    depend.matrix <- data.table(expand.grid(A.1 = 0:1,
                                            L.1 = 0:1))

    depend.matrix[, state := 1:.N]

    for (state.jj in depend.matrix[, unique(state)]) {

        tmp3.jj <- copy(tmp3)[, A.1 := depend.matrix[state == state.jj, A.1]][, L.1 := depend.matrix[state == state.jj, L.1]]

        tmp3.jj[, exp.1 := exp(predict(fit.cox.1, newdata=tmp3.jj, type="lp"))]
        tmp3.jj[, exp.z := exp(predict(fit.cox.z, newdata=tmp3.jj, type="lp"))]
        tmp3.jj[, exp.l := exp(predict(fit.cox.l, newdata=tmp3.jj, type="lp"))]
        
        tmp3[, (paste0("P.1.", state.jj)) := tmp3.jj[["dhazard.1"]]*tmp3.jj[["exp.1"]]]
        tmp3[, (paste0("P.z.", state.jj)) := alpha*tmp3.jj[["dhazard.z"]]*tmp3.jj[["exp.z"]]]
        tmp3[, (paste0("P.l.", state.jj)) := tmp3.jj[["dhazard.l"]]*tmp3.jj[["exp.l"]]]

        tmp3[A.1 == depend.matrix[state == state.jj, A.1] & L.1 == depend.matrix[state == state.jj, L.1],
             state := state.jj]

    }

    
    #--------------------------------    
    #-- compute clever weights:
    
    tmp3[, C := cumsum(1*(time == time.obs & delta == 0 & time %in% dt[["time"]])), by = "id"]
    tmp3[, C.1 := c(0, C[-.N]), by = "id"]
    
    if (length(Aname)>0) {
        tmp3[, probA := predict(fitA, newdata = tmp3, type = "response")]
        tmp3[, clever.weight := (Aobs == get(Aname))/((probA^(Aobs == 1)*(1-probA)^(Aobs == 0)))]
    } else {
        tmp3[, clever.weight := 1]
    }
    
    tmp3[, clever.weight := clever.weight*(C.1 == 0)/surv.0.1]

    tmp3[, clever.weight.alpha := alpha^A.1*exp(-(alpha-1)*cum.hazard.z.1)]
    tmp3[, summary(clever.weight.alpha)]

    tmp3[, final.time := max(time.obs), by = "id"]

    tmp.weight <- tmp3[time <= final.time, max(clever.weight.alpha), by = "id"]
    q100.weight <- as.numeric(quantile(tmp.weight[[2]], p = 1))
    q99.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.99))
    q975.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.975))
    q95.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.95))
    q90.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.9))
    q80.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.8))
    q70.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.7))
    q60.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.6))
    q50.weight <- as.numeric(quantile(tmp.weight[[2]], p = 0.5))

    if (truncate.weights>0 & q100.weight>truncate.weights) {
        message("NB: weights will be truncated")
        no.truncated <- tmp3[time <= final.time & clever.weight.alpha > truncate.weights, length(unique(id))]
        tmp3[time <= final.time & clever.weight.alpha > truncate.weights, clever.weight.alpha := truncate.weights]
    } else {
        if (truncate.weights>0) {
            no.truncated <- 0
        } else {
            no.truncated <- NA
        }
    }
    if (browse) browser()

    states <- copy(depend.matrix)
    setkey(states, state, A.1, L.1)

    for (iter in 1:max.iter) {

        setkey(tmp3, id, time)
        dt_list <- split(tmp3[, !(names(tmp3) %in% c("Z",
                                                     "clever.Z.D0", "clever.Z.D1",
                                                     "clever.Z.Z0", "clever.Z.Z1",
                                                     "clever.Z.L0", "clever.Z.L1")),
                              with = FALSE], by = "id", keep.by = TRUE)

        if (parameter == "target") {

            print(paste0("iter = ", iter))

            t2 <- system.time({
                dt_list <- mclapply(
                    dt_list,
                    compute_Z_and_clever_per_id_LZ,
                    states = states,
                    mc.cores = min(detectCores()-1, parallelize.Z)
                )
            })

            tmp3 <- rbindlist(dt_list)
            tmp3[, at.risk := time <= final.time]

            #-- current estimator for target parameter:  
            target.est <- mean(tmp3[, Z[1], by = "id"][[2]])

        } else {

            compute_Z_and_clever_per_id_LZ_auxiliary <- function(dt_id,
                                                                 states) {
                compute_Z_and_clever_per_id_LZ(dt_id = dt_id,
                                               states = states,
                                               parameter = "auxiliary")
            }

            t2 <- system.time({
                dt_list <- mclapply(
                    dt_list,
                    compute_Z_and_clever_per_id_LZ_auxiliary,
                    states = states,
                    mc.cores = min(detectCores()-1, parallelize.Z)
                )
            })

            tmp3 <- rbindlist(dt_list)
            tmp3[, at.risk := (time <= final.time & at.risk.z)]
            
            #-- current estimator for auxiliary parameter:  
            target.est <- mean(tmp3[, Z[1], by = "id"][[2]])

        }

        if (iter == 1) {
            g.est <- target.est
        }

        eic <- tmp3[at.risk == 1, sum(clever.weight*clever.weight.alpha*(clever.Z.D1-clever.Z.D0)*((delta == 1) - P.1)) +
                                  sum(alpha*clever.weight*clever.weight.alpha*(clever.Z.Z1-clever.Z.Z0)*((delta == 2) - at.risk.z*P.z)) +
                                  sum(clever.weight*clever.weight.alpha*(clever.Z.L1-clever.Z.L0)*((delta == 3) - at.risk.l*P.l)) +
                                  Z[1]-target.est,
                    by = "id"][[2]]

        if (iter == 1) {
            eic.init <- copy(eic)
            target.se <- sqrt(mean(eic^2/n))
        }

        print(paste0("eic equation solved at = ", abs(mean(eic))))

        if (abs(mean(eic)) <= target.se/(log(n))) break(print(paste0("finished after ", iter, " iterations")))
        
        if (target) {

            #-- otherwise update: 
            target.fun.D <- function(eps) {
                mean(tmp3[at.risk == 1, sum(clever.weight*clever.weight.alpha*(clever.Z.D1-clever.Z.D0)*((delta == 1) - P.1*exp(eps))), by = "id"][[2]])}
            target.fun.Z <- function(eps) {
                mean(tmp3[at.risk == 1, sum(alpha*clever.weight*clever.weight.alpha*(clever.Z.Z1-clever.Z.Z0)*((delta == 2) - at.risk.z*P.z*exp(eps))), by = "id"][[2]])}
            target.fun.L <- function(eps) {
                mean(tmp3[at.risk == 1, sum(clever.weight*clever.weight.alpha*(clever.Z.L1-clever.Z.L0)*((delta == 3) - at.risk.l*P.l*exp(eps))), by = "id"][[2]])}

            eps.D <- nleqslv(0.00, target.fun.D)$x
            eps.Z <- nleqslv(0.00, target.fun.Z)$x
            eps.L <- nleqslv(0.00, target.fun.L)$x

            if (verbose) print(paste0("eps.D = ", eps.D))
            if (verbose) print(paste0("eps.Z = ", eps.Z))
            if (verbose) print(paste0("eps.L = ", eps.L))

            tmp3[, P.1 := P.1*exp(eps.D)]
            tmp3[, P.z := P.z*exp(eps.Z)]
            tmp3[, P.l := P.l*exp(eps.L)]

            tmp3[, cum.hazard.z := cumsum(at.risk.z*P.z), by = "id"]
            tmp3[, cum.hazard.z.1 := c(0, cum.hazard.z[-.N]), by = "id"]

            tmp3[, clever.weight.alpha := alpha^A.1*exp(-(alpha-1)*cum.hazard.z.1)]
            
            for (state.jj in depend.matrix[, unique(state)]) {
                tmp3[, (paste0("P.1.", state.jj)) := get(paste0("P.1.", state.jj))*exp(eps.D)]
                tmp3[, (paste0("P.z.", state.jj)) := get(paste0("P.z.", state.jj))*exp(eps.Z)]
                tmp3[, (paste0("P.l.", state.jj)) := get(paste0("P.l.", state.jj))*exp(eps.L)]
            }
            
        }
    }

    if (output.eic) {
        return(list(estimate =
                        c(est = target.est, se = target.se,
                          g.est = g.est, iter = iter, 
                          eic.solved.at = abs(mean(eic)),
                          no.truncated = no.truncated,
                          q50.weight = q50.weight,
                          q60.weight = q60.weight,
                          q70.weight = q70.weight,
                          q80.weight = q80.weight,
                          q90.weight = q90.weight,
                          q95.weight = q95.weight,
                          q975.weight = q975.weight,
                          q99.weight = q99.weight,
                          q100.weight = q100.weight),
                    eic = eic.init))
    } else {
        return(c(est = target.est, se = target.se,
                 g.est = g.est, iter = iter, 
                 eic.solved.at = abs(mean(eic)),
                 no.truncated = no.truncated,
                 q50.weight = q50.weight,
                 q60.weight = q60.weight,
                 q70.weight = q70.weight,
                 q80.weight = q80.weight,
                 q90.weight = q90.weight,
                 q95.weight = q95.weight,
                 q975.weight = q975.weight,
                 q99.weight = q99.weight,
                 q100.weight = q100.weight))
    }
}
    
compute_Z_and_clever_per_id_LZ <- function(dt_id,
                                           states,
                                           P1.pref = "P.1.",
                                           Pl.pref = "P.l.",
                                           Pz.pref = "P.z.",
                                           state.idx.col = "state",
                                           parameter = "target",
                                           compute.clever = FALSE
                                           ) {

    S  <- nrow(states)
    Tn <- nrow(dt_id)

    is_target <- (parameter == "target")

    next_state <- function(a, l, da, dl) {
        states[A.1 == pmin(a + da, 1) &
               L.1 == pmin(l + dl, 1), state]
    }

    states_Z1 <- states[A.1 == 1, state]

    gamma <- data.table(
        state = states$state,
        s00 = states$state,
        s10 = mapply(next_state, states$A.1, states$L.1, 1, 0),
        s01 = mapply(next_state, states$A.1, states$L.1, 0, 1)
    )

    # hazard increment matrices: Tn x S
    P1 <- as.matrix(dt_id[, paste0(P1.pref, 1:S), with = FALSE])
    Pl <- as.matrix(dt_id[, paste0(Pl.pref, 1:S), with = FALSE])
    Pz <- as.matrix(dt_id[, paste0(Pz.pref, 1:S), with = FALSE])

    Z_by_time     <- vector("list", Tn)
    Znext_by_time <- vector("list", Tn)

    Z_row        <- numeric(Tn)
    clever_Z1    <- numeric(Tn)
    clever_Z0    <- numeric(Tn)
    clever_L1    <- numeric(Tn)
    clever_L0    <- numeric(Tn)
    clever_D1    <- numeric(Tn)
    clever_D0    <- numeric(Tn)

    # terminal condition
    if (is_target) {
        Z_T <- P1[Tn, ]
    } else {
        Z_T <- Pz[Tn, ]
        Z_T[states_Z1] <- 1
    }
    Z_by_time[[Tn]]     <- Z_T
    Znext_by_time[[Tn]] <- Z_T

    # backward recursion
    for (tt in (Tn - 1):1) {

        Z_next <- Z_by_time[[tt + 1]]

        P00 <- 1 - P1[tt, ] - Pl[tt, ] - Pz[tt, ]

        Z_t <- P00 * Z_next[gamma$s00] +
            Pz[tt, ] * Z_next[gamma$s10] +
            Pl[tt, ] * Z_next[gamma$s01]

        if (is_target) {
            Z_t <- Z_t + P1[tt, ]
        } else {
            Z_t[states_Z1] <- 1
        }

        Z_by_time[[tt]]     <- Z_t
        Znext_by_time[[tt]] <- Z_next
    }

    state_vec <- dt_id[[state.idx.col]]

    gamma_s00 <- gamma$s00
    gamma_s10 <- gamma$s10
    gamma_s01 <- gamma$s01
   
    # extract observed path + clever covariates
    for (tt in seq_len(Tn)) {

        s <- state_vec[tt]

        Znext <- Znext_by_time[[tt]]

        s_Z1    <- gamma_s10[s]
        s_Z0    <- gamma_s00[s]
        s_L1    <- gamma_s01[s]

        Z_row[tt] <- Z_by_time[[tt]][s]

        clever_D0[tt] <- Znext[s]

        if (is_target) {
            clever_D1[tt] <- 1
            clever_Z1[tt] <- Znext[s_Z1]
        } else {
            clever_D1[tt] <- 0
            clever_Z1[tt] <- 1
        }

        clever_Z0[tt] <- Znext[s_Z0] #* (1 - Pl[tt, s]) +
            #Znext[s_L1] * Pl[tt, s]

        clever_L0[tt] <- Znext[s_Z0] #* (1 - Pz[tt, s]) +
            #Znext[s_Z1] * Pz[tt, s]

        clever_L1[tt] <- Znext[s_L1]
    }

    cbind(
        dt_id,
        Z = Z_row,
        clever.Z.D0 = clever_D0,
        clever.Z.D1 = clever_D1,
        clever.Z.Z0 = clever_Z0,
        clever.Z.Z1 = clever_Z1,
        clever.Z.L0 = clever_L0,
        clever.Z.L1 = clever_L1
    )
}

######################################################################





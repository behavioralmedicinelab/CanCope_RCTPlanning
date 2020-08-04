################################################################################
####                                                                        ####
####                             Packages                                   ####
####                                                                        ####
################################################################################

library(MASS)
library(RcppEigen)
library(parallel)
library(doParallel)
library(doRNG)
library(foreach)
library(ggplot2)
library(ggpubr)
library(scales)
library(compiler)
library(data.table)
library(JWileymisc)
library(emmeans)
library(mice)


################################################################################
####                                                                        ####
####                    Simulate Data - Prior to RCT                        ####
####                                                                        ####
################################################################################

## simulate data and missingness
## this is a rough initial step just to make sure the main analysis code is
## properly setup and works
Sigma <- matrix(.3, nrow = 6, ncol = 6)
diag(Sigma) <- 1
Sigma <- Sigma * 6^2

N <- 150

percMiss <- c(.00, .12, .04, .04, .04, .10)
xmiss <- matrix(0, nrow = N, ncol = 6)
set.seed(12345)
for (i in 2:length(percMiss)) {
  n <- N - sum(xmiss[, i - 1])
  OK <- xmiss[, i - 1] == 0
  xmiss[OK, i] <- rbinom(n = n, size = 1, prob = (percMiss[i] / (n / N)))
  xmiss[!OK, i] <- 1
}

round(colSums(xmiss) / N, 2)
N - colSums(xmiss)

set.seed(12345)
mind <- mvrnorm(n = N/2,
  mu = c(58, 55, 54, 53, 50, 52), Sigma = Sigma)
lifestyle <- mvrnorm(n = N/2,
  mu = c(58, 56, 55, 55, 55, 54), Sigma = Sigma)
Xanx <- as.data.table(rbind(mind, lifestyle))
setnames(Xanx, new = paste0("Anx_", 1:6))


set.seed(23456)
mind <- mvrnorm(n = N/2,
  mu = c(59, 56, 54, 52, 50, 51), Sigma = Sigma)
lifestyle <- mvrnorm(n = N/2,
  mu = c(58, 56, 55, 55, 55, 54), Sigma = Sigma)
Xdep <- as.data.table(rbind(mind, lifestyle))
setnames(Xdep, new = paste0("Dep_", 1:6))

set.seed(34567)
mind <- mvrnorm(n = N/2,
  mu = c(54, 53, 52, 51, 50, 51), Sigma = Sigma)
lifestyle <- mvrnorm(n = N/2,
  mu = c(54, 53, 52, 51, 49, 50), Sigma = Sigma)
Xqol <- as.data.table(rbind(mind, lifestyle))
setnames(Xqol, new = paste0("QoL_", 1:6))

Sigma <- matrix(.3, nrow = 6, ncol = 6)
diag(Sigma) <- 1
Sigma <- Sigma * 14^2
set.seed(45678)
mind <- mvrnorm(n = N/2,
  mu = c(46, 40, 38, 36, 32, 33), Sigma = Sigma)
lifestyle <- mvrnorm(n = N/2,
  mu = c(46, 44, 43, 42, 41, 42), Sigma = Sigma)
Xder <- as.data.table(rbind(mind, lifestyle))
setnames(Xder, new = paste0("DERS_", 1:6))

dtrue <- cbind(
    ID = 1:N, Group = rep(c("CM", "CL"), each = N/2),
    One = 1L,
    STRATA_DERS = rbinom(N, size = 1, prob = .5),
    STRATA_DEPANX = rbinom(N, size = 1, prob = .25),
    Xanx, Xdep, Xqol, Xder)
d <- copy(dtrue)

## add unique dropout by time
for (i in 1:6) {
  d[, paste0("Anx_", i) := fifelse(xmiss[, i] == 1, NA_real_, get(paste0("Anx_", i)))]
  d[, paste0("Dep_", i) := fifelse(xmiss[, i] == 1, NA_real_, get(paste0("Dep_", i)))]
  d[, paste0("QoL_", i) := fifelse(xmiss[, i] == 1, NA_real_, get(paste0("QoL_", i)))]
  d[, paste0("DERS_", i) := fifelse(xmiss[, i] == 1, NA_real_, get(paste0("DERS_", i)))]
}

N - colSums(is.na(d))


################################################################################
####                                                                        ####
####                         Multiple Imputation                            ####
####                                                                        ####
################################################################################

if (FALSE) {
set.seed(12345)
## predictor matrix
pmat <- matrix(1L, 29, 29)
## do not use or predict ID or One (intercept constant)
pmat[, 1] <- 0L
pmat[1, ] <- 0L
pmat[, 3] <- 0L
pmat[3, ] <- 0L
## do not predict group or STRATA but can use group to predict
pmat[, 2] <- 0L
pmat[, 4] <- 0L
pmat[, 5] <- 0L
## do not predict variable from itself
diag(pmat) <- 0L

dimp <- parlmice(
  d,
  n.core = 10,
  n.imp.core = 10,
  cluster.seed = 12345,
  cl.type = "PSOCK",
  maxit = 100,
  predictorMatrix = pmat)

saveRDS(dimp, file = "dimp.RDS", compress = "xz")
} else {
  dimp <- readRDS("dimp.RDS")
}

################################################################################
####                                                                        ####
####                         Aim 1 - Analysis Code                          ####
####                                                                        ####
################################################################################


fitRCT <- function(dv, time, adj = c("none", "baseline")) {
  tmpw <- copy(d)
  fulldv <- paste0(dv, "_", time)
  if (dv %in% c("Dep", "Anx")) {
    covs <- "STRATA_DERS"
  } else if (dv %in% "DERS") {
    covs <- "STRATA_DEPANX"
  } else {
    covs <- "STRATA_DERS + STRATA_DEPANX"
  }

  if (time == 1L) {
    form <- switch(adj,
                   none = paste0(dv, "_", time, " ~ 0 + One"),
                   baseline = paste0(dv, "_", time, " ~ 0 + One + ", covs))
  } else {
    form <- switch(adj,
                   none = paste0(dv, "_", time, " ~ 1 + Group"),
                   baseline = paste0(dv, "_", time, " ~ 1 + Group + ", dv, "_1", " + ", covs))
  }

  m <- with(dimp, lm(as.formula(form)))
  pm <- pool(m)
  spm <- summary(pm, type = "all", conf.int=TRUE)

  if (time == 1L) {
    d <- NA_real_
    diff <- structure(list(term = structure(2L, .Label = c("(Intercept)",
                                                           "GroupCM", "baseline"), class = "factor"),
                           m = dimp$m, estimate = 0, std.error = 0,
                           statistic = 0, df = NA_real_, p.value = 1,
                           `2.5 %` = 0, `97.5 %` = 0, riv = 0, lambda = 0,
                           fmi = 0, ubar = 0, b = 0, t = 0, dfcom = NA_real_),
                      row.names = 2L, class = "data.frame")
    mm <- as.data.frame(summary(emmeans(m, "One", weights = "cells",
                                        df = subset(pm$pooled, term == "One")$df), inf=TRUE))
    mm <- rbind(mm, mm)
    names(mm)[1] <- "Group"
    mm[[1]] <- c("CL", "CM")
  } else {
    d <- mean(sapply(m$analyses, function(x) coef(x)[["GroupCM"]]) /
              sapply(m$analyses, sigma))
    diff <- subset(spm, term == "GroupCM")
    mm <- as.data.frame(summary(emmeans(m, "Group", weights = "cells",
                          df = subset(pm$pooled, term == "GroupCM")$df), inf=TRUE))
  }

  out <- data.table(
    DV = dv,
    Time = time,
    Adj = adj,
    Type = c("CL", "CM", "DIFF"),
    N = c(
      sum(!is.na(tmpw[Group == "CL"][[fulldv]])),
      sum(!is.na(tmpw[Group == "CM"][[fulldv]])),
      sum(!is.na(tmpw[[fulldv]]))),
    Est = c(mm$emmean, diff$estimate),
    SE = c(mm$SE, diff[["std.error"]]),
    LL = c(mm$lower.CL, diff[["2.5 %"]]),
    UL = c(mm$upper.CL, diff[["97.5 %"]]),
    ES = c(NA_real_, NA_real_, d),
    P = c(mm$p.value, diff[["p.value"]]))
  out[, Out := sprintf("%s; %d; %s; %s; %0.2f [%0.2f, %0.2f] %s%s",
                       DV, Time, Adj, Type,
                       round(Est, 2),
                       round(LL, 2),
                       round(UL, 2),
                       formatPval(P, d = 3, sd = 3, includeP=TRUE),
                       fifelse(Type == "DIFF", sprintf(" ES = %0.2f", round(ES, 2)), ""))]
  out[, Formula := form]

  return(copy(out))
}


## Makes a grid to loop through variables and time points
reg.vars <- expand.grid(var = c("Anx", "Dep", "DERS", "QoL"),
                        tp = c(1L:6L),
                        adjust = c("none", "baseline"),
                        stringsAsFactors = FALSE)

cl <- makeCluster(12)
clusterEvalQ(cl, {
  library(data.table)
  library(emmeans)
  library(mice)
  library(JWileymisc)
})
clusterExport(cl, c("reg.vars", "d", "dimp", "fitRCT"))

resall <- parLapplyLB(cl = cl, 1:nrow(reg.vars), function(i) {
    fitRCT(reg.vars$var[i], reg.vars$tp[i], adj = reg.vars$adjust[i])
})
resall <- do.call(rbind, resall)

## pd <- position_dodge(width = .3)
## ggplot(resall[Type=="DIFF"], aes(Time, Est, ymin = LL, ymax = UL, colour = Adj)) +
##   geom_rect(aes(xmin = 1.2, xmax = 4.8, ymin = -Inf, ymax = Inf),
##             colour = NA, fill = "grey90") +
##   geom_hline(yintercept = 0, linetype = 3, colour="grey50") +
##   geom_pointrange(position = pd) +
##   scale_x_continuous(
##     "Time Point", breaks = 1:6,
##     labels = c(expression(t[-1]), expression(t[1]), expression(t[2]), expression(t[3]),
##                expression(t[4]), expression(t[5]))) +
##   facet_wrap(~DV, scales = "free") +
##   theme_pubr() +
##   color_palette("jco")


#### Difficulty with Emotion Regulation (Primary Outcome)

tmpd <- resall[Type!="DIFF" & Adj == "baseline" & DV == "DERS"]
tmpd[, Ylab := factor(Time,
                      levels = 1:6,
                      labels = c("Baseline", "Post M1", "Post M2", "Post M3", "Post M4", "Follow-Up"))]
setnames(tmpd, old = "Type", new = "Group")
tmpd2 <- resall[Type=="DIFF" & Adj == "baseline" & DV == "DERS"]
tmpd2[, ES := gsub("0\\.", "\\.", sprintf("ES = %0.2f", round(abs(ES), 2)))]
tmpd2[, Max := max(tmpd$UL, na.rm=TRUE)]
tmpd2[, Min := min(tmpd$LL, na.rm=TRUE)]
tmpd2[, Max := Max + abs(Max - Min) * .02]
tmpd2[, Min := Min - abs(Max - Min) * .03]
labs <- tmpd[, sprintf("%s\n%s%d%s\n%s%d%s",
                       Ylab[1],
                       fifelse(Time==1, "CL  n: ", ""),
                       N[Group=="CL"],
                       fifelse(Time==1, "           ", ""),
                       fifelse(Time==1, "CM n: ", ""),
                       N[Group=="CM"],
                       fifelse(Time==1, "           ", "")),
             by = Time]
p.der <- ggplot(tmpd, aes(Time, Est, ymin = LL, ymax = UL)) +
  geom_rect(aes(xmin = 1.3, xmax = 4.8, ymin = -Inf, ymax = Inf),
            colour = NA, fill = "grey90") +
  ## geom_hline(yintercept = 50, linetype = 3, colour = "grey60") +
  geom_text(aes(x = Time, y = Max, label = star(P)), data = tmpd2) +
  ## geom_text(aes(x = Time, y = Max, label = formatPval(P, includeP=TRUE)), data = tmpd2) +
  geom_text(aes(x = Time, y = Min, label = ES), data = tmpd2) +
  geom_pointrange(aes(colour = Group, shape = Group), position = pd) +
  scale_x_continuous("",
                     breaks = labs$Time,
                     labels = labs$V1) +
  scale_y_continuous("Difficulty with Emotion Regulation") +
  scale_shape_manual(values=c("CL" = 17, "CM" = 19)) +
  theme_pubr() +
  color_palette("jco") +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(tmpd2$Min[1], tmpd2$Max[1]), expand = TRUE)
print(p.der)

#### Anxiety Symptoms (Secondary Outcome)

tmpd <- resall[Type!="DIFF" & Adj == "baseline" & DV == "Anx"]
tmpd[, Ylab := factor(Time,
                      levels = 1:6,
                      labels = c("Baseline", "Post M1", "Post M2", "Post M3", "Post M4", "Follow-Up"))]
setnames(tmpd, old = "Type", new = "Group")
tmpd2 <- resall[Type=="DIFF" & Adj == "baseline" & DV == "Anx"]
tmpd2[, ES := gsub("0\\.", "\\.", sprintf("ES = %0.2f", round(abs(ES), 2)))]
tmpd2[, Max := max(tmpd$UL, na.rm=TRUE)]
tmpd2[, Min := min(tmpd$LL, na.rm=TRUE)]
tmpd2[, Max := Max + abs(Max - Min) * .02]
tmpd2[, Min := Min - abs(Max - Min) * .03]
labs <- tmpd[, sprintf("%s\n%s%d%s\n%s%d%s",
                       Ylab[1],
                       fifelse(Time==1, "CL  n: ", ""),
                       N[Group=="CL"],
                       fifelse(Time==1, "           ", ""),
                       fifelse(Time==1, "CM n: ", ""),
                       N[Group=="CM"],
                       fifelse(Time==1, "           ", "")),
             by = Time]
p.anx <- ggplot(tmpd, aes(Time, Est, ymin = LL, ymax = UL)) +
  geom_rect(aes(xmin = 1.3, xmax = 4.8, ymin = -Inf, ymax = Inf),
            colour = NA, fill = "grey90") +
  geom_hline(yintercept = 50, linetype = 3, colour = "grey60") +
  geom_text(aes(x = Time, y = Max, label = star(P)), data = tmpd2) +
  ## geom_text(aes(x = Time, y = Max, label = formatPval(P, includeP=TRUE)), data = tmpd2) +
  geom_text(aes(x = Time, y = Min, label = ES), data = tmpd2) +
  geom_pointrange(aes(colour = Group, shape = Group), position = pd) +
  scale_x_continuous("",
                     breaks = labs$Time,
                     labels = labs$V1) +
  scale_y_continuous("Anxiety TScores") +
  scale_shape_manual(values=c("CL" = 17, "CM" = 19)) +
  theme_pubr() +
  color_palette("jco") +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(tmpd2$Min[1], tmpd2$Max[1]), expand = TRUE)
print(p.anx)


#### Depression Symptoms (Secondary Outcome)

tmpd <- resall[Type!="DIFF" & Adj == "baseline" & DV == "Dep"]
tmpd[, Ylab := factor(Time,
                      levels = 1:6,
                      labels = c("Baseline", "Post M1", "Post M2", "Post M3", "Post M4", "Follow-Up"))]
setnames(tmpd, old = "Type", new = "Group")
tmpd2 <- resall[Type=="DIFF" & Adj == "baseline" & DV == "Dep"]
tmpd2[, ES := gsub("0\\.", "\\.", sprintf("ES = %0.2f", round(abs(ES), 2)))]
tmpd2[, Max := max(tmpd$UL, na.rm=TRUE)]
tmpd2[, Min := min(tmpd$LL, na.rm=TRUE)]
tmpd2[, Max := Max + abs(Max - Min) * .02]
tmpd2[, Min := Min - abs(Max - Min) * .03]
labs <- tmpd[, sprintf("%s\n%s%d%s\n%s%d%s",
                       Ylab[1],
                       fifelse(Time==1, "CL  n: ", ""),
                       N[Group=="CL"],
                       fifelse(Time==1, "           ", ""),
                       fifelse(Time==1, "CM n: ", ""),
                       N[Group=="CM"],
                       fifelse(Time==1, "           ", "")),
             by = Time]
p.dep <- ggplot(tmpd, aes(Time, Est, ymin = LL, ymax = UL)) +
  geom_rect(aes(xmin = 1.3, xmax = 4.8, ymin = -Inf, ymax = Inf),
            colour = NA, fill = "grey90") +
  geom_hline(yintercept = 50, linetype = 3, colour = "grey60") +
  geom_text(aes(x = Time, y = Max, label = star(P)), data = tmpd2) +
  ## geom_text(aes(x = Time, y = Max, label = formatPval(P, includeP=TRUE)), data = tmpd2) +
  geom_text(aes(x = Time, y = Min, label = ES), data = tmpd2) +
  geom_pointrange(aes(colour = Group, shape = Group), position = pd) +
  scale_x_continuous("",
                     breaks = labs$Time,
                     labels = labs$V1) +
  scale_y_continuous("Depression TScores") +
  scale_shape_manual(values=c("CL" = 17, "CM" = 19)) +
  theme_pubr() +
  color_palette("jco") +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(tmpd2$Min[1], tmpd2$Max[1]), expand = TRUE)
print(p.dep)

#### Health Related Quality of Life (Secondary Outcome)

tmpd <- resall[Type!="DIFF" & Adj == "baseline" & DV == "QoL"]
tmpd[, Ylab := factor(Time,
                      levels = 1:6,
                      labels = c("Baseline", "Post M1", "Post M2", "Post M3", "Post M4", "Follow-Up"))]
setnames(tmpd, old = "Type", new = "Group")
tmpd2 <- resall[Type=="DIFF" & Adj == "baseline" & DV == "QoL"]
tmpd2[, ES := gsub("0\\.", "\\.", sprintf("ES = %0.2f", round(abs(ES), 2)))]
tmpd2[, Max := max(tmpd$UL, na.rm=TRUE)]
tmpd2[, Min := min(tmpd$LL, na.rm=TRUE)]
tmpd2[, Max := Max + abs(Max - Min) * .02]
tmpd2[, Min := Min - abs(Max - Min) * .03]
labs <- tmpd[, sprintf("%s\n%s%d%s\n%s%d%s",
                       Ylab[1],
                       fifelse(Time==1, "CL  n: ", ""),
                       N[Group=="CL"],
                       fifelse(Time==1, "           ", ""),
                       fifelse(Time==1, "CM n: ", ""),
                       N[Group=="CM"],
                       fifelse(Time==1, "           ", "")),
             by = Time]
p.qol <- ggplot(tmpd, aes(Time, Est, ymin = LL, ymax = UL)) +
  geom_rect(aes(xmin = 1.3, xmax = 4.8, ymin = -Inf, ymax = Inf),
            colour = NA, fill = "grey90") +
  geom_hline(yintercept = 50, linetype = 3, colour = "grey60") +
  geom_text(aes(x = Time, y = Max, label = star(P)), data = tmpd2) +
  ## geom_text(aes(x = Time, y = Max, label = formatPval(P, includeP=TRUE)), data = tmpd2) +
  geom_text(aes(x = Time, y = Min, label = ES), data = tmpd2) +
  geom_pointrange(aes(colour = Group, shape = Group), position = pd) +
  scale_x_continuous("",
                     breaks = labs$Time,
                     labels = labs$V1) +
  scale_y_continuous("QoL TScores") +
  scale_shape_manual(values=c("CL" = 17, "CM" = 19)) +
  theme_pubr() +
  color_palette("jco") +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(tmpd2$Min[1], tmpd2$Max[1]), expand = TRUE)
print(p.qol)

#### Put all figures together

ggarrange(
  p.der, p.anx, p.dep, p.qol,
  ncol = 2, nrow = 2,
  align = "hv", labels = c("A", "B", "C", "D"),
  legend = "bottom", common.legend = TRUE)

################################################################################
####                                                                        ####
####                             Packages                                   ####
####                                                                        ####
################################################################################

library(MASS)
library(ggplot2)
library(ggpubr)
library(scales)
library(data.table)
library(MplusAutomation)
library(mice)

################################################################################
####                                                                        ####
####                    Simulate Data - Prior to RCT                        ####
####                                                                        ####
################################################################################

## simulate data and missingness
## this is a rough initial step just to make sure the main analysis code is
## properly setup and works
Sigma <- matrix(.3, nrow = 5, ncol = 5)
diag(Sigma) <- 1
Sigma <- Sigma * 6^2

N <- 150

percMiss <- c(.00, .12, .04, .04, .04)
xmiss <- matrix(0, nrow = N, ncol = 5)
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
  mu = c(58, 55, 54, 53, 50), Sigma = Sigma)
lifestyle <- mvrnorm(n = N/2,
  mu = c(58, 56, 55, 55, 55), Sigma = Sigma)
Xbes <- as.data.table(rbind(mind, lifestyle))
setnames(Xbes, new = paste0("BES_", 1:5))


set.seed(23456)
mind <- mvrnorm(n = N/2,
  mu = c(59, 56, 54, 52, 50), Sigma = Sigma)
lifestyle <- mvrnorm(n = N/2,
  mu = c(58, 56, 55, 55, 55), Sigma = Sigma)
Xsmq <- as.data.table(rbind(mind, lifestyle))
setnames(Xsmq, new = paste0("SMQ_", 1:5))

dtrue <- cbind(
    ID = 1:N, Group = rep(c("CM", "CL"), each = N/2),
    One = 1L,
    STRATA_DERS = rbinom(N, size = 1, prob = .5),
    STRATA_DEPANX = rbinom(N, size = 1, prob = .25),
    Xbes, Xsmq)
d <- copy(dtrue)

## add unique dropout by time
for (i in 1:5) {
  d[, paste0("BES_", i) := fifelse(xmiss[, i] == 1, NA_real_, get(paste0("BES_", i)))]
  d[, paste0("SMQ_", i) := fifelse(xmiss[, i] == 1, NA_real_, get(paste0("SMQ_", i)))]
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
pmat <- matrix(1L, ncol(d), ncol(d))
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

saveRDS(dimp, file = "dimp_aim3.RDS", compress = "xz")
} else {
  dimp <- readRDS("dimp_aim3.RDS")
}


################################################################################
####                                                                        ####
####                         Aim 3 - Analysis Code                          ####
####                                                                        ####
################################################################################


dimpl <- complete(dimp, action = "long", include = TRUE)
dimpl$Group <- as.integer(dimpl$Group == "CM")
dimpl <- lapply(1:100, function(i) {
  subset(dimpl, .imp == i)
})

m <- mplusObject(
  TITLE = "Hypothesis 3A - Module 1 Mechanisms",
  ANALYSIS = "ESTIMATOR = ML;",
  MODEL = "
    BES_5 ON BES_4 Group STRATA_DEPANX STRATA_DERS (d1 - d4);
    BES_4 ON BES_3 Group STRATA_DEPANX STRATA_DERS (c1 - c4);
    BES_3 ON BES_2 Group STRATA_DEPANX STRATA_DERS (b1 - b4);
    BES_2 ON BES_1 Group STRATA_DEPANX STRATA_DERS (a1 - a4);
    [BES_1* BES_2* BES_3* BES_4* BES_5*] (i1 - i5)
    BES_1* BES_2* BES_3* BES_4* BES_5* (e1 - e5);
",
MODELCONSTRAINT = "
 NEW(avg diff2 diff3 diff4);
 !avg = (a2/sqrt(e2)) - (( (b2/sqrt(e3)) + (c2/sqrt(e4)) + (d2/sqrt(e5)))/3);
 !diff2 = (a2/sqrt(e2)) - (b2/sqrt(e3));
 !diff3 = (a2/sqrt(e2)) - (c2/sqrt(e4));
 !diff4 = (a2/sqrt(e2)) - (d2/sqrt(e5));
 avg = (a2) - (( (b2) + (c2) + (d2))/3);
 diff2 = (a2) - (b2);
 diff3 = (a2) - (c2);
 diff4 = (a2) - (d2);
",
OUTPUT = "STDY; CINTERVAL;",
usevariables = c("Group", "STRATA_DEPANX", "STRATA_DERS", paste0("BES_", 1:5)),
rdata = dimpl,
imputed = TRUE)
mfit <- mplusModeler(m, dataout = "sap_bes_implist.dat", run = TRUE)

mfit$results$summaries
coef(mfit)
coef(mfit, type = "stdy")


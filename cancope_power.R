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

doit <- cmpfun(function(n, r, d) {
  stopifnot(identical(n %/% 2, n / 2))
  g1 <- g2 <- n / 2
  X <- cbind(
    Int = 1,
    x = rnorm(n = n, mean = 0, sd = 1),
    g = rep(0:1, c(g1, g2)))
  B <- c(0, r, d)
  vhat <- var(X %*% B)
  y <- rnorm(n = n, mean = X %*% B, sd = sqrt(1 - vhat))
  m <- coef(summary(fastLm(X, y)))
  m["g", "Pr(>|t|)"]
})

set.seed(1234)
doit(100, r = .1, d = .2)

grid <- expand.grid(
  n = seq(80, 160, 2),
  r = seq(.3, .6, .1),
  d = c(.5, .6))

cl <- makeCluster(24)
clusterEvalQ(cl, {
  library(RcppEigen)
})
clusterExport(cl, c("doit", "grid"))

registerDoParallel(cl)

set.seed(1234)
res <- foreach(i = 1:nrow(grid), .combine = c) %dorng% list(replicate(10000, doit(n = grid$n[i], r = grid$r[i], d = grid$d[i])))

grid$alpha05 <- unlist(lapply(res, function(x) mean(x < .05)))

grid$r <- paste0("correlation over time: ", grid$r)
grid$d <- paste0("d = ", grid$d)

ggplot(subset(grid, n >= 0), aes(n, alpha05, linetype = r, colour = r)) +
  geom_hline(yintercept = .8) +
  geom_vline(xintercept = 110, linetype = 2) +
  stat_smooth(se=FALSE) +
  facet_wrap(~ d) +
  scale_colour_viridis_d() +
  scale_y_continuous("Power", breaks = c(.7, .8, .9, 1), labels = percent) +
  scale_x_continuous("Total Sample Size (N)",
                     breaks = c(80, 100, 110, 120, 140, 160),
                     labels = c(80, 100, "\n110", 120, 140, 160)) +
  theme_pubr()

### Models for binary outcome - common ###

setwd("~/Aim1")

library(dplyr)
library(reshape2)
library(ggplot2)
library(tvmediation)
library(logbin)
library(brglm)

set.seed(112322)

load("~/Aim1/aim1_dat.rda")

ID <- unique(aim1_dat[,c("id","treat")])
pop_split <- split(ID, ID$treat)

draw_sample <- function(n, data){
  data[sample((1:nrow(data)), n, replace=FALSE),]
}

R = 2000
ss = 500

mat_id <- matrix(nrow = R, ncol = ss)

mat_a0 <- matrix(nrow = R, ncol = 20)
mat_a1 <- matrix(nrow = R, ncol = 20)
mat_lm <- matrix(nrow = R, ncol = 20)
mat_x1m <- matrix(nrow = R, ncol = 20)
# mat_x2m <- matrix(nrow = R, ncol = 20)
mat_m <- matrix(nrow = R, ncol = 19)

mat_b0 <- matrix(nrow = R, ncol = 20)
mat_g <- matrix(nrow = R, ncol = 20)
mat_b1 <- matrix(nrow = R, ncol = 20)
mat_lo <- matrix(nrow = R, ncol = 20)
mat_x1o <- matrix(nrow = R, ncol = 20)
# mat_x2o <- matrix(nrow = R, ncol = 20)
mat_o <- matrix(nrow = R, ncol = 19)

mean.y <- as.data.frame(table(aim1_dat$time_pt, aim1_dat$Yc)) %>%
  filter(Var2 == 1)
start.p <- mean.y$Freq/nrow(ID)

for (r in 1:R) {
  print(paste("Simulation_",r))
  
  samp <- do.call("rbind", mapply(draw_sample, 
                                  n = c(ss/2,ss/2),
                                  data = pop_split, SIMPLIFY = FALSE))
  
  id <- as.vector(samp$id)
  mat_id[r,] <- id
  
  aim1 <- aim1_dat[aim1_dat$id %in% id,] %>%
    arrange(id, time_pt)
  
  t.seq <- sort(unique(aim1$time_pt))
  
  trt <- as.numeric(unique(aim1[ ,c("id","treat")])[,2])
  
  x1 <- as.numeric(unique(aim1[ ,c("id","x1")])[,2])
  # x2 <- as.numeric(unique(aim1[ ,c("id","x2")])[,2])
  
  mediator <- LongToWide(aim1$id, aim1$time_pt, aim1$M)
  outcome <- LongToWide(aim1$id, aim1$time_pt, aim1$Yc)
  tvcon <- LongToWide(aim1$id, aim1$time_pt, aim1$L)
  
  n <- length(trt)
  nm <- nrow(outcome)
  
  a0 <- vector()
  a1 <- vector()
  lm <- vector()
  x_1m <- vector()
  # x_2m <- vector()
  m <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      fit1 <- glm(factor(mediator[i,]) ~ factor(trt) + factor(tvcon[i,]) + x1,
                  family=binomial(link = "logit"), na.action=na.omit)
    } else {
      fit1 <- glm(factor(mediator[i,]) ~ factor(trt) + factor(tvcon[i,]) + x1 + factor(mediator[i-1,]),
                  family=binomial(link = "logit"), na.action=na.omit)
      m <- append(m, fit1$coefficients[[5]])
    }
    a0 <- append(a0, fit1$coefficients[[1]])
    a1 <- append(a1, fit1$coefficients[[2]])
    lm <- append(lm, fit1$coefficients[[3]])
    x_1m <- append(x_1m, fit1$coefficients[[4]])
    # x_2m <- append(x_2m, fit1$coefficients[[5]])
  }
  
  mat_a0[r,] <- a0
  mat_a1[r,] <- a1
  mat_lm[r,] <- lm
  mat_x1m[r,] <- x_1m
  # mat_x2m[r,] <- x_2m
  mat_m[r,] <- m
  
  b0 <- vector()
  g <- vector()
  b1 <- vector()
  lo <- vector()
  x_1o <- vector()
  # x_2o <- vector()
  o <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1,
      #                         family = poisson(link = "log"),
      #                         na.action=na.omit))
      fit2 <- logbin(factor(outcome[i,]) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1, 
                     start = c(log(start.p[i]), rep(0,4)), maxit = 25000,
                     method = "em",
                     na.action=na.omit)
    } else {
      # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,],
      #                         family = poisson(link = "log"),
      #                         na.action=na.omit))
      fit2 <- logbin(factor(outcome[i,]) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + factor(outcome[i-1,]), 
                     start = c(log(start.p[i]), rep(0,5)), maxit = 25000,
                     method = "em",
                     na.action=na.omit)
      o <- append(o, fit2$coefficients[[6]])
    }
    b0 <- append(b0, fit2$coefficients[[1]])
    g <- append(g, fit2$coefficients[[2]])
    b1 <- append(b1, fit2$coefficients[[3]])
    lo <- append(lo, fit2$coefficients[[4]])
    x_1o <- append(x_1o, fit2$coefficients[[5]])
    # x_2o <- append(x_2o, fit2$coefficients[[6]])
  }
  
  mat_b0[r,] <- b0
  mat_g[r,] <- g
  mat_b1[r,] <- b1
  mat_lo[r,] <- lo
  mat_x1o[r,] <- x_1o
  # mat_x2o[r,] <- x_2o
  mat_o[r,] <- o
}

results_continuous <- vector("list")

results_continuous$mat_id <- mat_id

results_continuous$mat_a0 <-  mat_a0
results_continuous$mat_a1 <- mat_a1
results_continuous$mat_lm <- mat_lm
results_continuous$mat_x1m <- mat_x1m
# results_continuous$mat_x2m <- mat_x2m
results_continuous$mat_m <- mat_m

results_continuous$mat_b0 <- mat_b0
results_continuous$mat_g <- mat_g
results_continuous$mat_b1 <- mat_b1
results_continuous$mat_lo <- mat_lo
results_continuous$mat_x1o <- mat_x1o
# results_continuous$mat_x2o <- mat_x2o
results_continuous$mat_o <- mat_o

save(results_continuous, file = paste("results_binary_500_modelc.rda", sep = ""))

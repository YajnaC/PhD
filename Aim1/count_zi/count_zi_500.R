### Models for count outcome - zero-inflated with log link for zero component using the Nelder-Mead method###

setwd("~/Aim1")

library(dplyr)
library(reshape2)
library(ggplot2)
library(tvmediation)
library(logbin)
library(brglm)
library(pscl)

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

mat_b0_count <- matrix(nrow = R, ncol = 20)
mat_g_count <- matrix(nrow = R, ncol = 20)
mat_b1_count <- matrix(nrow = R, ncol = 20)
mat_lo_count <- matrix(nrow = R, ncol = 20)
mat_x1o_count <- matrix(nrow = R, ncol = 20)
# mat_x2o_count <- matrix(nrow = R, ncol = 20)
mat_o_count <- matrix(nrow = R, ncol = 19)

mat_b0_zero <- matrix(nrow = R, ncol = 20)
mat_g_zero <- matrix(nrow = R, ncol = 20)
mat_b1_zero <- matrix(nrow = R, ncol = 20)
mat_lo_zero <- matrix(nrow = R, ncol = 20)
mat_x1o_zero <- matrix(nrow = R, ncol = 20)
# mat_x2o_zero <- matrix(nrow = R, ncol = 20)
mat_o_zero <- matrix(nrow = R, ncol = 19)

mean.y <- as.data.frame(table(aim1_dat$time_pt, aim1_dat$Cigc)) %>%
  filter(Var2 == 0)
# logit_pr <- (nrow(ID)-mean.y$Freq)/nrow(ID)
logit_pr <- (mean.y$Freq)/nrow(ID)
# start.p <- logit_pr/(1-logit_pr)
start.p <- logit_pr

for (r in 1:R) {
  print(paste("Simulation_",r))
  
  samp <- do.call("rbind", mapply(draw_sample, 
                                  n = c(ss/2,ss/2),
                                  data = pop_split, SIMPLIFY = FALSE))
  
  id <- as.vector(samp$id)
  mat_id[r,] <- id
  
  aim1 <- aim1_dat[aim1_dat$id %in% id,] %>%
    arrange(id, time_pt)
  
  # for(i in 1:nrow(aim1)){
  #   if(aim1$Cigc[i] == 0)
  #     aim1$SI[i] <- 0
  #   else
  #     aim1$SI[i] <- 1
  # }
  
  t.seq <- sort(unique(aim1$time_pt))
  
  trt <- as.numeric(unique(aim1[ ,c("id","treat")])[,2])
  
  x1 <- as.numeric(unique(aim1[ ,c("id","x1")])[,2])
  # x2 <- as.numeric(unique(aim1[ ,c("id","x2")])[,2])
  
  mediator <- LongToWide(aim1$id, aim1$time_pt, aim1$M)
  outcome <- LongToWide(aim1$id, aim1$time_pt, aim1$Cigc)
  tvcon <- LongToWide(aim1$id, aim1$time_pt, aim1$L)
  # si <- LongToWide(aim1$id, aim1$time_pt, aim1$SI)
  
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
  
  b0_count <- vector()
  g_count <- vector()
  b1_count <- vector()
  lo_count <- vector()
  x_1o_count <- vector()
  # x_2o_count <- vector()
  o_count <- vector()
  
  b0_zero <- vector()
  g_zero <- vector()
  b1_zero <- vector()
  lo_zero <- vector()
  x_1o_zero <- vector()
  # x_2o_zero <- vector()
  o_zero <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      start.coeff <- glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1,
                         family = poisson(link = "log"),
                         na.action=na.omit)
      start.c.p <- exp(coef(start.coeff)[1])
      # x <- model.matrix(start.coeff)
      
      start.coeff1 <- logbin(factor(outcome[i,]==0) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1,
                             # family = binomial(link = "log"),
                             start = c(log(start.p[i]), rep(0,4)), maxit = 25000,
                             method = "em",
                             na.action=na.omit)
      start.z.p <- 1 - exp(coef(start.coeff1)[1])
      
      # z <- model.matrix(start.coeff1)
      # y <- outcome[i,]
      # 
      # nll <- function(par) {
      #   lambda <- exp(x %*% par[1:5])
      #   ziprob <- plogis(z %*% par[6:10])
      #   zidens <- ziprob * (y == 0) + (1 - ziprob) * dpois(y, lambda = lambda)
      #   -sum(log(zidens))
      # }
      # start.c.p <- optim(c(coef(start.coeff), coef(start.coeff1)), nll, method = "BFGS")
      
      fit2 <- zeroinfl(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 | factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1,
                       na.action=na.omit, method = "Nelder-Mead",  link = "log",
                       start = list(count = c(log(start.c.p), rep(0, 4)), 
                                    zero = c(log(start.z.p), rep(0, 4))))
    } else {
      # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,] ,
      #                         family = poisson(link = "log"),
      #                         na.action=na.omit))
      # start.c.p <- exp(start.coeff)/(1+exp(start.coeff))
      
      start.coeff <- glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,],
                         family = poisson(link = "log"),
                         na.action=na.omit)
      start.c.p <- exp(coef(start.coeff)[1])
      # x <- model.matrix(start.coeff)
      
      start.coeff1 <- logbin(factor(outcome[i,]==0) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,],
                             # family = binomial(link = "log"),
                             start = c(log(start.p[i]), rep(0,5)), maxit = 25000,
                             method = "em",
                             na.action=na.omit)
      start.z.p <- 1 - exp(coef(start.coeff1)[1])
      
      fit2 <- zeroinfl(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,] | factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,],
                       na.action=na.omit, method = "Nelder-Mead", link = "log",
                       start = list(count = c(log(start.c.p), rep(0, 5)), 
                                    zero = c(log(start.z.p), rep(0, 5))))
      
      o_count <- append(o_count, fit2$coefficients$count[[6]])
      o_zero <- append(o_zero, fit2$coefficients$zero[[6]])
    }
    b0_count <- append(b0_count, fit2$coefficients$count[[1]])
    g_count <- append(g_count, fit2$coefficients$count[[2]])
    b1_count <- append(b1_count, fit2$coefficients$count[[3]])
    lo_count <- append(lo_count, fit2$coefficients$count[[4]])
    x_1o_count <- append(x_1o_count, fit2$coefficients$count[[5]])
    # x_2o_count <- append(x_2o_count, fit2$coefficients$count[[6]])
    
    b0_zero <- append(b0_zero, fit2$coefficients$zero[[1]])
    g_zero <- append(g_zero, fit2$coefficients$zero[[2]])
    b1_zero <- append(b1_zero, fit2$coefficients$zero[[3]])
    lo_zero <- append(lo_zero, fit2$coefficients$zero[[4]])
    x_1o_zero <- append(x_1o_zero, fit2$coefficients$zero[[5]])
    # x_2o_zero <- append(x_2o_zero, fit2$coefficients$zero[[6]])
  }
  
  mat_b0_count[r,] <- b0_count
  mat_g_count[r,] <- g_count
  mat_b1_count[r,] <- b1_count
  mat_lo_count[r,] <- lo_count
  mat_x1o_count[r,] <- x_1o_count
  # mat_x2o_count[r,] <- x_2o_count
  mat_o_count[r,] <- o_count
  
  mat_b0_zero[r,] <- b0_zero
  mat_g_zero[r,] <- g_zero
  mat_b1_zero[r,] <- b1_zero
  mat_lo_zero[r,] <- lo_zero
  mat_x1o_zero[r,] <- x_1o_zero
  # mat_x2o_zero[r,] <- x_2o_zero
  mat_o_zero[r,] <- o_zero
  
}

results_continuous <- vector("list")

results_continuous$mat_id <- mat_id

results_continuous$mat_a0 <-  mat_a0
results_continuous$mat_a1 <- mat_a1
results_continuous$mat_lm <- mat_lm
results_continuous$mat_x1m <- mat_x1m
# results_continuous$mat_x2m <- mat_x2m
results_continuous$mat_m <- mat_m

results_continuous$mat_b0_count <- mat_b0_count
results_continuous$mat_g_count <- mat_g_count
results_continuous$mat_b1_count <- mat_b1_count
results_continuous$mat_lo_count <- mat_lo_count
results_continuous$mat_x1o_count <- mat_x1o_count
# results_continuous$mat_x2o_count <- mat_x2o_count
results_continuous$mat_o_count <- mat_o_count

results_continuous$mat_b0_zero <- mat_b0_zero
results_continuous$mat_g_zero <- mat_g_zero
results_continuous$mat_b1_zero <- mat_b1_zero
results_continuous$mat_lo_zero <- mat_lo_zero
results_continuous$mat_x1o_zero <- mat_x1o_zero
# results_continuous$mat_x2o_zero <- mat_x2o_zero
results_continuous$mat_o_zero <- mat_o_zero

save(results_continuous, file = paste("results_count_zi_500_v2.rda", sep = ""))

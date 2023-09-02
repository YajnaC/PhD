### Although the file name says binary_c_bootstrap.R this is in fact for the count zero inflated model ###

library(dplyr)
library(reshape2)
library(ggplot2)
library(tvmediation)
library(logbin)
library(brglm)
library(pscl)

start.time <- Sys.time()

setwd("~/Aim1/count_zi")

results <- get(load("~/Aim1/count_zi/results_count_zi_500_v2.rda"))

load("~/Aim1/aim1_dat.rda")

### True coefficients of association ###

j <- seq(from = 1, to = 20, length.out = 20)

alpha0 <- -0.4
alpha1 <- exp((-0.25 -(2/j)))

lm <- -1.2
lo <- 0.36

x1m <- 0.05
x1o <- -0.02

m <- 1.25
o <- 0.05

beta0 <- -0.1 + 0.5
gamma_A <- -0.1*exp(1/j)
beta1_c <- -(1/exp(1/(j^3)))

pr1_beta0 <- -0.1 - 1.17
pr1_trt <- 0.2
pr1_m <- 0.9
pr1_l <- -0.36
pr1_x1 <- 0.02
pr1_o <- -0.05

true_coeff <- vector("list",30)
true_coeff[[1]] <- alpha0
true_coeff[[2]] <- alpha1
true_coeff[[3]] <- lm
true_coeff[[4]] <- x1m
true_coeff[[5]] <- m
true_coeff[[6]] <- beta0
true_coeff[[7]] <- exp(gamma_A)
true_coeff[[8]] <- beta1_c
true_coeff[[9]] <- lo
true_coeff[[10]] <- x1o
true_coeff[[11]] <- o
true_coeff[[12]] <- (exp(alpha1)/(1+exp(alpha1)))*exp(beta1_c)

for(i in 1:20){
  if(i == 1){
    ### Overall Indirect Effect ###
    ## At baseline L(t+1)=1; C=2.25
    orNIE_l1 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + lm + x1m*2.25)))
    
    ## At baseline L(t+1)=0; C=2.25
    orNIE_l0 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[13]][i] <- orNIE_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[14]][i] <- orNIE_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[15]][i] <- orNIE_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[16]][i] <- orNIE_l0
    
    
    ### Indirect Effect through susceptibility latent variable ###
    ## At baseline L(t+1)=1; C=2.25
    rrIES_l1 <- ((1 + exp(pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + lm + x1m*2.25)))
    
    ## At baseline L(t+1)=0; C=2.25
    rrIES_l0 <- ((1 + exp(pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[17]][i] <- rrIES_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[18]][i] <- rrIES_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[19]][i] <- rrIES_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[20]][i] <- rrIES_l0
    
    ### Indirect Effect not through susceptibility latent variable ###
    ## At baseline L(t+1)=1; C=2.25
    rrIENS_l1 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25)))
    
    ## At baseline L(t+1)=0; C=2.25
    rrIENS_l0 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[21]][i] <- rrIENS_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[22]][i] <- rrIENS_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[23]][i] <- rrIENS_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[24]][i] <- rrIENS_l0
    
  }else{
    ### Overall Indirect Effect ###
    ## At t+1 M(t)=1; L(t+1)=1; C=2.25
    orNIE_m1l1 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(alpha0 + lm + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + lm + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=2.25
    orNIE_m0l1 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + lm + x1m*2.25)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=2.25
    orNIE_m1l0 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(alpha0 + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=2.25
    orNIE_m0l0 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[13]][i] <- orNIE_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[14]][i] <- orNIE_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[15]][i] <- orNIE_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[16]][i] <- orNIE_m0l0
    
    
    ### Indirect Effect through susceptibility latent variable ###
    ## At t+1 M(t)=1; L(t+1)=1; C=2.25
    rrIES_m1l1 <- ((1 + exp(pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + lm + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=2.25
    rrIES_m0l1 <- ((1 + exp(pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + lm + x1m*2.25)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=2.25
    rrIES_m1l0 <- ((1 + exp(pr1_m + alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=2.25
    rrIES_m0l0 <- ((1 + exp(pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + pr1_m + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[17]][i] <- rrIES_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[18]][i] <- rrIES_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[19]][i] <- rrIES_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[20]][i] <- rrIES_m0l0
    
    ### Indirect Effect not through susceptibility latent variable ###
    ## At t+1 M(t)=1; L(t+1)=1; C=2.25
    rrIENS_m1l1 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(alpha0 + lm + x1m*2.25 + m)))/((1 + exp(pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=2.25
    rrIENS_m0l1 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(pr1_m + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=2.25
    rrIENS_m1l0 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(alpha0 + x1m*2.25 + m)))/((1 + exp(pr1_m + alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=2.25
    rrIENS_m0l0 <- ((1 + exp(beta1_c[i] + pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(pr1_m + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[21]][i] <- rrIENS_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[22]][i] <- rrIENS_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[23]][i] <- rrIENS_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[24]][i] <- rrIENS_m0l0
  }
}

true_coeff[[25]] <- pr1_beta0
true_coeff[[26]] <- pr1_trt
true_coeff[[27]] <- pr1_m
true_coeff[[28]] <- pr1_l
true_coeff[[29]] <- pr1_x1
true_coeff[[30]] <- pr1_o


CI <- vector("list", 30)
# CI[[1]] <- vector("list",2000)
# CI[[2]] <- vector("list",2000)
# CI[[3]] <- vector("list",2000)
# CI[[4]] <- vector("list",2000)
# CI[[5]] <- vector("list",2000)
# CI[[6]] <- vector("list",2000)
# CI[[7]] <- vector("list",2000)
# CI[[8]] <- vector("list",2000)
# CI[[9]] <- vector("list",2000)
# CI[[10]] <- vector("list",2000)
# CI[[11]] <- vector("list",2000)
# CI[[12]] <- vector("list",2000)
# CI[[13]] <- vector("list",2000)
# CI[[14]] <- vector("list",2000)
# CI[[15]] <- vector("list",2000)
# CI[[16]] <- vector("list",2000)
# CI[[17]] <- vector("list",2000)
# CI[[18]] <- vector("list",2000)
# CI[[19]] <- vector("list",2000)
# CI[[20]] <- vector("list",2000)
# CI[[21]] <- vector("list",2000)
# CI[[22]] <- vector("list",2000)
# CI[[23]] <- vector("list",2000)
# CI[[24]] <- vector("list",2000)
# CI[[25]] <- vector("list",2000)
# CI[[26]] <- vector("list",2000)
# CI[[27]] <- vector("list",2000)
# CI[[28]] <- vector("list",2000)
# CI[[29]] <- vector("list",2000)
# CI[[30]] <- vector("list",2000)


CI[[1]] <- vector("list",1)
CI[[2]] <- vector("list",1)
CI[[3]] <- vector("list",1)
CI[[4]] <- vector("list",1)
CI[[5]] <- vector("list",1)
CI[[6]] <- vector("list",1)
CI[[7]] <- vector("list",1)
CI[[8]] <- vector("list",1)
CI[[9]] <- vector("list",1)
CI[[10]] <- vector("list",1)
CI[[11]] <- vector("list",1)
CI[[12]] <- vector("list",1)
CI[[13]] <- vector("list",1)
CI[[14]] <- vector("list",1)
CI[[15]] <- vector("list",1)
CI[[16]] <- vector("list",1)
CI[[17]] <- vector("list",1)
CI[[18]] <- vector("list",1)
CI[[19]] <- vector("list",1)
CI[[20]] <- vector("list",1)
CI[[21]] <- vector("list",1)
CI[[22]] <- vector("list",1)
CI[[23]] <- vector("list",1)
CI[[24]] <- vector("list",1)
CI[[25]] <- vector("list",1)
CI[[26]] <- vector("list",1)
CI[[27]] <- vector("list",1)
CI[[28]] <- vector("list",1)
CI[[29]] <- vector("list",1)
CI[[30]] <- vector("list",1)


# for(k in 1:nrow(results$mat_id)) {
# k = 1
  print(paste("Starting bootstrap for Sample_",k))
  
  cov_avg <- vector()
  
  # samp_id <- results[["mat_id"]][k,]
  id <- results$mat_id[k,]
  coeff_data <- aim1_dat[aim1_dat$id %in% id,] %>%
    arrange(id, time_pt)
  
  cov_avg <- append(cov_avg, mean(coeff_data$x1))
  
  R = 300
  
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
  
  mat_b0_zero <- matrix(nrow = R, ncol = 20)
  mat_g_zero <- matrix(nrow = R, ncol = 20)
  mat_b1_zero <- matrix(nrow = R, ncol = 20)
  mat_lo_zero <- matrix(nrow = R, ncol = 20)
  mat_x1o_zero <- matrix(nrow = R, ncol = 20)
  # mat_x2o_zero <- matrix(nrow = R, ncol = 20)
  mat_o_zero <- matrix(nrow = R, ncol = 19)
  
  # mean.y <- as.data.frame(table(coeff_data$time_pt, coeff_data$Yc)) %>%
  #   filter(Var2 == 1)
  # start.p <- mean.y$Freq/length(id)
  
  mean.y <- as.data.frame(table(coeff_data$time_pt, coeff_data$Cigc)) %>%
    filter(Var2 == 0)
  # logit_pr <- (nrow(ID)-mean.y$Freq)/nrow(ID)
  logit_pr <- (mean.y$Freq)/length(id)
  # start.p <- logit_pr/(1-logit_pr)
  start.p <- logit_pr
  
  for (r in 1:R) {
    # print(paste("CI_Simulation_",r))
    
    t.seq <- sort(unique(coeff_data$time_pt))
    
    trt <- as.numeric(unique(coeff_data[ ,c("id","treat")])[,2])
    
    x1 <- as.numeric(unique(coeff_data[ ,c("id","x1")])[,2])
    # x2 <- as.numeric(unique(coeff_data[ ,c("id","x2")])[,2])
    
    mediator <- LongToWide(coeff_data$id, coeff_data$time_pt, coeff_data$M)
    outcome <- LongToWide(coeff_data$id, coeff_data$time_pt, coeff_data$Cigc)
    tvcon <- LongToWide(coeff_data$id, coeff_data$time_pt, coeff_data$L)
    
    n <- length(trt)
    nm <- nrow(outcome)
    
    index1 <- sample(1:n, size=n, replace=TRUE)
    
    a0 <- vector()
    a1 <- vector()
    lm <- vector()
    x_1m <- vector()
    # x_2m <- vector()
    m <- vector()
    
    for(i in 1:nm){
      if(i == 1){
        fit1 <- glm(factor(mediator[i,index1]) ~ factor(trt[index1]) + factor(tvcon[i,index1]) + x1[index1],
                                                  family=binomial(link = "logit"))
      } else {
        fit1 <- glm(factor(mediator[i,index1]) ~ factor(trt[index1]) + factor(tvcon[i,index1]) + x1[index1] + factor(mediator[i-1,index1]), 
                                                  family=binomial(link = "logit"))
        m <- append(m, fit1$coefficients[[5]])
      }
      a0 <- append(a0, fit1$coefficients[[1]])
      a1 <- append(a1, fit1$coefficients[[2]])
      lm <- append(lm, fit1$coefficients[[3]])
      x_1m <- append(x_1m, fit1$coefficients[[4]])
      # x_2m <- append(x_2m, fit1$coefficients[[5]])
    }
    
    mat_a0[r,] <- a0
    # mat_a1[r,] <- pnorm(predict(loess(a1 ~ t.seq[1:nm], span = 0.5, degree = 1)))
    mat_a1[r,] <- predict(loess(a1 ~ t.seq[1:nm], span = 0.5, degree = 1))
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
        start.coeff <- glm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1],
                           family = poisson(link = "log"),
                           na.action=na.omit)
        start.c.p <- exp(coef(start.coeff)[1])
        # x <- model.matrix(start.coeff)
        
        start.coeff1 <- logbin(factor(outcome[i,index1]==0) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1],
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
        
        
        fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1],
                         na.action=na.omit, method = "Nelder-Mead",  link = "log",
                         start = list(count = c(log(start.c.p), rep(0, 4)), 
                                      zero = c(log(start.z.p), rep(0, 4))))
      } else {
        # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,] ,
        #                         family = poisson(link = "log"),
        #                         na.action=na.omit))
        # start.c.p <- exp(start.coeff)/(1+exp(start.coeff))
        
        start.coeff <- glm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1],
                           family = poisson(link = "log"),
                           na.action=na.omit)
        start.c.p <- exp(coef(start.coeff)[1])
        # x <- model.matrix(start.coeff)
        
        start.coeff1 <- logbin(factor(outcome[i,index1]==0) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1],
                               # family = binomial(link = "log"),
                               start = c(log(start.p[i]), rep(0,5)), maxit = 25000,
                               method = "em",
                               na.action=na.omit)
        start.z.p <- 1 - exp(coef(start.coeff1)[1])
        
        fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1],
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
    
    mat_b0[r,] <- b0_count
    mat_g[r,] <- predict(loess(g_count ~ t.seq[1:nm], span = 0.5, degree = 1))
    mat_b1[r,] <- predict(loess(b1_count ~ t.seq[1:nm], span = 0.5, degree = 1))
    mat_lo[r,] <- lo_count
    mat_x1o[r,] <- x_1o_count
    # mat_x2o[r,] <- x_2o_count
    mat_o[r,] <- o_count
    
    mat_b0_zero[r,] <- b0_zero
    mat_g_zero[r,] <- g_zero
    mat_b1_zero[r,] <- b1_zero
    mat_lo_zero[r,] <- lo_zero
    mat_x1o_zero[r,] <- x_1o_zero
    # mat_x2o_zero[r,] <- x_2o_zero
    mat_o_zero[r,] <- o_zero
    
  }
  
  effect_interest <- vector("list",30)
  effect_interest[[1]] <- mat_a0
  effect_interest[[2]] <- mat_a1
  effect_interest[[3]] <- mat_lm
  effect_interest[[4]] <- mat_x1m
  effect_interest[[5]] <- mat_m
  effect_interest[[6]] <- mat_b0
  effect_interest[[7]] <- exp(mat_g)
  effect_interest[[8]] <- mat_b1
  effect_interest[[9]] <- mat_lo
  effect_interest[[10]] <- mat_x1o
  effect_interest[[11]] <- mat_o
  
  effect_interest[[12]] <- (exp(mat_a1)/(1+exp(mat_a1)))*exp(mat_b1)
  
  effect_interest[[13]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[14]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[15]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[16]] <- matrix(nrow = R, ncol = 20)
  
  effect_interest[[17]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[18]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[19]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[20]] <- matrix(nrow = R, ncol = 20)
  
  effect_interest[[21]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[22]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[23]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[24]] <- matrix(nrow = R, ncol = 20)
  
  effect_interest[[25]] <- mat_b0_zero
  effect_interest[[26]] <- mat_g_zero
  effect_interest[[27]] <- mat_b1_zero
  effect_interest[[28]] <- mat_lo_zero
  effect_interest[[29]] <- mat_x1o_zero
  effect_interest[[30]] <- mat_o_zero
  
  for (q in 1:nrow(mat_a0)) {
    IE11 <- vector()
    IE01 <- vector()
    IE10 <- vector()
    IE00 <- vector()
    
    IES11 <- vector()
    IES01 <- vector()
    IES10 <- vector()
    IES00 <- vector()
    
    IENS11 <- vector()
    IENS01 <- vector()
    IENS10 <- vector()
    IENS00 <- vector()
    
    for(i in 1:20){
      if(i == 1){
        ### Overall Indirect Effect ###
        ## At baseline L(t+1)=1; C=cov_avg
        orNIE_l1 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## At baseline L(t+1)=0; C=cov_avg
        orNIE_l0 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## IE at M(t)=1; L(t+1)=1; C=cov_avg
        IE11 <- append(IE11, orNIE_l1)
        
        ## IE at M(t)=0; L(t+1)=1; C=cov_avg
        IE01 <- append(IE01, orNIE_l1)
        
        ## IE at M(t)=1; L(t+1)=0; C=cov_avg
        IE10 <- append(IE10, orNIE_l0)
        
        ## IE at M(t)=0; L(t+1)=0; C=cov_avg
        IE00 <- append(IE00, orNIE_l0)
        
        
        ### Indirect Effect through susceptibility latent variable ###
        ## At baseline L(t+1)=1; C=cov_avg
        rrIES_l1 <- ((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## At baseline L(t+1)=0; C=cov_avg
        rrIES_l0 <- ((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## IE at M(t)=1; L(t+1)=1; C=cov_avg
        IES11 <- append(IES11, rrIES_l1)
        
        ## IE at M(t)=0; L(t+1)=1; C=cov_avg
        IES01 <- append(IES01, rrIES_l1)
        
        ## IE at M(t)=1; L(t+1)=0; C=cov_avg
        IES10 <- append(IES10, rrIES_l0)
        
        ## IE at M(t)=0; L(t+1)=0; C=cov_avg
        IES00 <- append(IES00, rrIES_l0)
        
        
        ### Indirect Effect not through susceptibility latent variable ###
        ## At baseline L(t+1)=1; C=cov_avg
        rrIENS_l1 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## At baseline L(t+1)=0; C=cov_avg
        rrIENS_l0 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## IE at M(t)=1; L(t+1)=1; C=cov_avg
        IENS11 <- append(IENS11, rrIENS_l1)
        
        ## IE at M(t)=0; L(t+1)=1; C=cov_avg
        IENS01 <- append(IENS01, rrIENS_l1)
        
        ## IE at M(t)=1; L(t+1)=0; C=cov_avg
        IENS10 <- append(IENS10, rrIENS_l0)
        
        ## IE at M(t)=0; L(t+1)=0; C=cov_avg
        IENS00 <- append(IENS00, rrIENS_l0)
      
      }else{
        ### Overall Indirect Effect ###
        ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
        orNIE_m1l1 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))
        
        ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
        orNIE_m0l1 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
        orNIE_m1l0 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))
        
        ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
        orNIE_m0l0 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## IE at M(t)=1; L(t+1)=1; C=cov_avg
        IE11 <- append(IE11, orNIE_m1l1)
        
        ## IE at M(t)=0; L(t+1)=1; C=cov_avg
        IE01 <- append(IE01, orNIE_m0l1)
        
        ## IE at M(t)=1; L(t+1)=0; C=cov_avg
        IE10 <- append(IE10, orNIE_m1l0)
        
        ## IE at M(t)=0; L(t+1)=0; C=cov_avg
        IE00 <- append(IE00, orNIE_m0l0)
        
        
        ### Indirect Effect through susceptibility latent variable ###
        ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
        rrIES_m1l1 <- ((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))
        
        ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
        rrIES_m0l1 <- ((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
        rrIES_m1l0 <- ((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))
        
        ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
        rrIES_m0l0 <- ((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## IE at M(t)=1; L(t+1)=1; C=cov_avg
        IES11 <- append(IES11, rrIES_m1l1)
        
        ## IE at M(t)=0; L(t+1)=1; C=cov_avg
        IES01 <- append(IES01, rrIES_m0l1)
        
        ## IE at M(t)=1; L(t+1)=0; C=cov_avg
        IES10 <- append(IES10, rrIES_m1l0)
        
        ## IE at M(t)=0; L(t+1)=0; C=cov_avg
        IES00 <- append(IES00, rrIES_m0l0)
        
        
        ### Indirect Effect not through susceptibility latent variable ###
        ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
        rrIENS_m1l1 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))/((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))
        
        ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
        rrIENS_m0l1 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
        rrIENS_m1l0 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))/((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])))
        
        ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
        rrIENS_m0l0 <- ((1 + exp(mat_b1[q,i] + mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))/((1 + exp(mat_b1_zero[q,i] + mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))*(1 + exp(mat_b1[q,i] + mat_a0[q,i] + mat_x1m[q,i]*cov_avg)))
        
        ## IE at M(t)=1; L(t+1)=1; C=cov_avg
        IENS11 <- append(IENS11, rrIENS_m1l1)
        
        ## IE at M(t)=0; L(t+1)=1; C=cov_avg
        IENS01 <- append(IENS01, rrIENS_m0l1)
        
        ## IE at M(t)=1; L(t+1)=0; C=cov_avg
        IENS10 <- append(IENS10, rrIENS_m1l0)
        
        ## IE at M(t)=0; L(t+1)=0; C=cov_avg
        IENS00 <- append(IENS00, rrIENS_m0l0)
        
      }
    }
    effect_interest[[13]][q,] <- IE11
    effect_interest[[14]][q,] <- IE01
    effect_interest[[15]][q,] <- IE10
    effect_interest[[16]][q,] <- IE00
    
    effect_interest[[17]][q,] <- IES11
    effect_interest[[18]][q,] <- IES01
    effect_interest[[19]][q,] <- IES10
    effect_interest[[20]][q,] <- IES00
    
    effect_interest[[21]][q,] <- IENS11
    effect_interest[[22]][q,] <- IENS01
    effect_interest[[23]][q,] <- IENS10
    effect_interest[[24]][q,] <- IENS00
    
  }
  
  for(z in 1:length(effect_interest)){
    mat_interest <- effect_interest[[z]]
    
    # if(z == 5 || z == 11){
    if(z %in% c(5,11,30)){
      quantiles <- matrix(NA, nrow=nm-1, ncol=2)

      lower <- 0.025
      upper <- 1 - lower

      for(i in 1:nm-1){
        quantiles[i,1] <- quantile(mat_interest[,i], c(lower), na.rm=TRUE)
        quantiles[i,2] <- quantile(mat_interest[,i], c(upper), na.rm=TRUE)
      }
      
    }else{
      quantiles <- matrix(NA, nrow=nm, ncol=2)
      
      lower <- 0.025
      upper <- 1 - lower
      
      for(i in 1:nm){
        quantiles[i,1] <- quantile(mat_interest[,i], c(lower), na.rm=TRUE)
        quantiles[i,2] <- quantile(mat_interest[,i], c(upper), na.rm=TRUE)
      }
      
    }
    CI[[z]][[k]] <- quantiles
  }
  
  print(paste("End of bootstrap for Sample_",k))
  save(CI, file = paste("CI_count_zi_500_",k,".rda", sep = ""))
  
# }

end.time <- Sys.time()
total.time <- end.time - start.time
print(sprintf("Bootstrap complete. Elapsed time = %.3f secs", 
              as.numeric(total.time, units = "secs")))

# save(CI, file = paste("CI_count_nzi_500.rda", sep = ""))


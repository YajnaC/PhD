library(tvmediation)
library(dplyr)

start.time <- Sys.time()

setwd("~/Aim1/continuous")

results <- get(load("~/Aim1/continuous/results_cont_500_dat.rda"))

load("~/Aim1/aim1_dat.rda")

j <- seq(from = 1, to = 20, length.out = 20)

alpha0 <- -0.4
alpha1 <- exp((-0.25 -(2/j)))

lm <- -1.2
lo <- 0.36

x1m <- -0.05
x1o <- -0.02

m <- 1.25
o <- 0.5

beta0_cont <- 2.3
gamma_A <- -0.1*exp(1/j)
beta1_c <- -(1/exp(1/(j^3)))

true_coeff <- vector("list",16)
true_coeff[[1]] <- alpha0
true_coeff[[2]] <- alpha1
true_coeff[[3]] <- lm
true_coeff[[4]] <- x1m
true_coeff[[5]] <- m
true_coeff[[6]] <- beta0_cont
true_coeff[[7]] <- gamma_A
true_coeff[[8]] <- beta1_c
true_coeff[[9]] <- lo
true_coeff[[10]] <- x1o
true_coeff[[11]] <- o
true_coeff[[12]] <- exp(alpha1)*beta1_c
# true_coeff[[12]] <- (exp(alpha1)/(1+exp(alpha1)))*beta1_c

for(i in 1:20){
  if(i == 1){
    ## At baseline L(t+1)=1; C=2.25
    prob1_l1 <- exp(alpha0 + alpha1[i] + lm + x1m*2.25)/(1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))
    prob0_l1 <- exp(alpha0 + lm + x1m*2.25)/(1 + exp(alpha0 + lm + x1m*2.25))
    
    ## At baseline L(t+1)=0; C=2.25
    prob1_l0 <- exp(alpha0 + alpha1[i] + x1m*2.25)/(1 + exp(alpha0 + alpha1[i] + x1m*2.25))
    prob0_l0 <- exp(alpha0 + x1m*2.25)/(1 + exp(alpha0 + x1m*2.25))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[13]][i] <- beta1_c[i]*(prob1_l1 - prob0_l1)
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[14]][i] <- beta1_c[i]*(prob1_l1 - prob0_l1)
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[15]][i] <- beta1_c[i]*(prob1_l0 - prob0_l0)
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[16]][i] <- beta1_c[i]*(prob1_l0 - prob0_l0)
    
  }else{
    ## At t+1 M(t)=1; L(t+1)=1; C=2.25
    prob1_m1l1 <- exp(alpha0 + alpha1[i] + lm + x1m*2.25 + m)/(1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25 + m))
    prob0_m1l1 <- exp(alpha0 + lm + x1m*2.25 + m)/(1 + exp(alpha0 + lm + x1m*2.25 + m))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=2.25
    prob1_m0l1 <- exp(alpha0 + alpha1[i] + lm + x1m*2.25)/(1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))
    prob0_m0l1 <- exp(alpha0 + lm + x1m*2.25)/(1 + exp(alpha0 + lm + x1m*2.25))
    
    ## At t+1 M(t)1; L(t+1)=0; C=2.25
    prob1_m1l0 <- exp(alpha0 + alpha1[i] + x1m*2.25 + m)/(1 + exp(alpha0 + alpha1[i] + x1m*2.25 + m))
    prob0_m1l0 <- exp(alpha0 + x1m*2.25 + m)/(1 + exp(alpha0 + x1m*2.25 + m))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=2.25
    prob1_m0l0 <- exp(alpha0 + alpha1[i] + x1m*2.25)/(1 + exp(alpha0 + alpha1[i] + x1m*2.25))
    prob0_m0l0 <- exp(alpha0 + x1m*2.25)/(1 + exp(alpha0 + x1m*2.25))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[13]][i] <- beta1_c[i]*(prob1_m1l1 - prob0_m1l1)
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[14]][i] <- beta1_c[i]*(prob1_m0l1 - prob0_m0l1)
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[15]][i] <- beta1_c[i]*(prob1_m1l0 - prob0_m1l0)
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[16]][i] <- beta1_c[i]*(prob1_m0l0 - prob0_m0l0)
  }
}

CI <- vector("list", 16)
CI[[1]] <- vector("list",2000)
CI[[2]] <- vector("list",2000)
CI[[3]] <- vector("list",2000)
CI[[4]] <- vector("list",2000)
CI[[5]] <- vector("list",2000)
CI[[6]] <- vector("list",2000)
CI[[7]] <- vector("list",2000)
CI[[8]] <- vector("list",2000)
CI[[9]] <- vector("list",2000)
CI[[10]] <- vector("list",2000)
CI[[11]] <- vector("list",2000)
CI[[12]] <- vector("list",2000)
CI[[13]] <- vector("list",2000)
CI[[14]] <- vector("list",2000)
CI[[15]] <- vector("list",2000)
CI[[16]] <- vector("list",2000)

for (k in 1:nrow(results$mat_id)) {
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
  
  for (r in 1:R) {
    # print(paste("CI_Simulation_",r))
    
    t.seq <- sort(unique(coeff_data$time_pt))
    
    trt <- as.numeric(unique(coeff_data[ ,c("id","treat")])[,2])
    
    x1 <- as.numeric(unique(coeff_data[ ,c("id","x1")])[,2])
    # x2 <- as.numeric(unique(coeff_data[ ,c("id","x2")])[,2])
    
    mediator <- LongToWide(coeff_data$id, coeff_data$time_pt, coeff_data$M)
    outcome <- LongToWide(coeff_data$id, coeff_data$time_pt, coeff_data$Ycrav)
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
    
    b0 <- vector()
    g <- vector()
    b1 <- vector()
    lo <- vector()
    x_1o <- vector()
    x_2o <- vector()
    o <- vector()
    
    for(i in 1:nm){
      if(i == 1){
        fit2 <- lm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1])
      } else {
        fit2 <- lm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1])
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
    mat_g[r,] <- predict(loess(g ~ t.seq[1:nm], span = 0.5, degree = 1))
    mat_b1[r,] <- predict(loess(b1 ~ t.seq[1:nm], span = 0.5, degree = 1))
    mat_lo[r,] <- lo
    mat_x1o[r,] <- x_1o
    # mat_x2o[r,] <- x_2o
    mat_o[r,] <- o
  }
  
  effect_interest <- vector("list",16)
  effect_interest[[1]] <- mat_a0
  effect_interest[[2]] <- mat_a1
  effect_interest[[3]] <- mat_lm
  effect_interest[[4]] <- mat_x1m
  effect_interest[[5]] <- mat_m
  effect_interest[[6]] <- mat_b0
  effect_interest[[7]] <- mat_g
  effect_interest[[8]] <- mat_b1
  effect_interest[[9]] <- mat_lo
  effect_interest[[10]] <- mat_x1o
  effect_interest[[11]] <- mat_o

  effect_interest[[12]] <- exp(mat_a1)*mat_b1
  
  effect_interest[[13]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[14]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[15]] <- matrix(nrow = R, ncol = 20)
  effect_interest[[16]] <- matrix(nrow = R, ncol = 20)
  
  
for (q in 1:nrow(mat_a0)) {
  IE11 <- vector()
  IE01 <- vector()
  IE10 <- vector()
  IE00 <- vector()
  
  for(i in 1:20){
    if(i == 1){
      ## At baseline L(t+1)=1; C=cov_avg
      prob1_l1 <- exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))
      prob0_l1 <- exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))
      
      ## At baseline L(t+1)=0; C=cov_avg
      prob1_l0 <- exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))
      prob0_l0 <- exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg))
      
      ## IE at M(t)=1; L(t+1)=1; C=cov_avg
      IE11 <- append(IE11, (mat_b1[q,i]*(prob1_l1 - prob0_l1)))
      
      ## IE at M(t)=0; L(t+1)=1; C=cov_avg
      IE01 <- append(IE01, (mat_b1[q,i]*(prob1_l1 - prob0_l1)))
      
      ## IE at M(t)=1; L(t+1)=0; C=cov_avg
      IE10 <- append(IE10, (mat_b1[q,i]*(prob1_l0 - prob0_l0)))
      
      ## IE at M(t)=0; L(t+1)=0; C=cov_avg
      IE00 <- append(IE00, (mat_b1[q,i]*(prob1_l0 - prob0_l0)))
      
    }else{
      ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
      prob1_m1l1 <- exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])/(1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))
      prob0_m1l1 <- exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])/(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      prob1_m0l1 <- exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))
      prob0_m0l1 <- exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_lm[q,i] + mat_x1m[q,i]*cov_avg))
      
      ## At t+1 M(t)1; L(t+1)=0; C=cov_avg
      prob1_m1l0 <- exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])/(1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))
      prob0_m1l0 <- exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1])/(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg + mat_m[q,i-1]))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      prob1_m0l0 <- exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_a1[q,i] + mat_x1m[q,i]*cov_avg))
      prob0_m0l0 <- exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg)/(1 + exp(mat_a0[q,i] + mat_x1m[q,i]*cov_avg))
      
      ## IE at M(t)=1; L(t+1)=1; C=cov_avg
      IE11 <- append(IE11, (mat_b1[q,i]*(prob1_m1l1 - prob0_m1l1)))
      
      ## IE at M(t)=0; L(t+1)=1; C=cov_avg
      IE01 <- append(IE01, (mat_b1[q,i]*(prob1_m0l1 - prob0_m0l1)))
      
      ## IE at M(t)=1; L(t+1)=0; C=cov_avg
      IE10 <- append(IE10, (mat_b1[q,i]*(prob1_m1l0 - prob0_m1l0)))
      
      ## IE at M(t)=0; L(t+1)=0; C=cov_avg
      IE00 <- append(IE00, (mat_b1[q,i]*(prob1_m0l0 - prob0_m0l0)))
    }
  }
  effect_interest[[13]][q,] <- IE11
  effect_interest[[14]][q,] <- IE01
  effect_interest[[15]][q,] <- IE10
  effect_interest[[16]][q,] <- IE00
}

for(z in 1:length(effect_interest)){
    mat_interest <- effect_interest[[z]]

    if(z == 5 || z == 11){
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
}

end.time <- Sys.time()
total.time <- end.time - start.time
print(sprintf("Bootstrap complete. Elapsed time = %.3f secs", 
              as.numeric(total.time, units = "secs")))

save(CI, file = paste("CI_cont_500_final.rda", sep = ""))


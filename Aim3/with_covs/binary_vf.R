library(dplyr)
library(reshape2)
library(ggplot2)
library(tvmediation)
library(logbin)
library(brglm)

### Complete data analysis ###
setwd("~/Aim3")
load("~/Aim3/original.data.imp.rda")

set.seed(1234)
index1 <- sample(1:5, 1)
original.data <- original.data.imp[[index1]]

# Varenicline vs Nicotine Patch Only #
data <- original.data %>%
  filter(ConditionID != 4) %>%
  mutate(treat = ifelse(ConditionID == 2, 0, 1)) %>%
  arrange(SubjectID, DaysFromTQD)

# cNRT vs Nicotine Patch Only #
# data <- original.data %>%
#   filter(ConditionID != 3) %>%
#   mutate(treat = ifelse(ConditionID == 2, 0, 1)) %>%
#   arrange(SubjectID, DaysFromTQD)

t.seq <- sort(unique(data$DaysFromTQD))

trt <- as.numeric(unique(data[ ,c("SubjectID","ConditionID")])[,2])

x1 <- as.numeric(unique(data[ ,c("SubjectID","sdi_score")])[,2])
x2 <- as.numeric(unique(data[ ,c("SubjectID","FTND_Score")])[,2])
# x3 <- as.numeric(unique(data[ ,c("SubjectID","AgeAtConsent")])[,2])
x3 <- as.numeric(unique(data[ ,c("SubjectID","Gender")])[,2])

mean.y <- as.data.frame(table(data$DaysFromTQD, data$Smoked)) %>%
  filter(Var2 == 1)
start.p <- mean.y$Freq/length(trt)

### Full Adherence - Mediator ###
mediator <- LongToWide(data$SubjectID, data$DaysFromTQD, data$adh_full)

### Partial Adherence - Mediator ###
# mediator <- LongToWide(data$SubjectID, data$DaysFromTQD, data$adh_part)

tvcon <- LongToWide(data$SubjectID, data$DaysFromTQD, data$StressEventToday)

### Binary Outcome - Daily Smoking Status ###
outcome <- LongToWide(data$SubjectID, data$DaysFromTQD, data$Smoked)

n <- length(trt)
nm <- nrow(outcome)

### All covariates ###

a0 <- vector()
a1 <- vector()
lm <- vector()
x_1m <- vector()
x_2m <- vector()
x_3m <- vector()
m <- vector()

for(i in 1:nm){
  if(i == 1){
    fit1 <- glm(mediator[i,] ~ factor(trt) + factor(tvcon[i,]) + x1 + x2 + as.factor(x3),
                family=binomial(link = "logit"))
  } else {
    fit1 <- glm(mediator[i,] ~ factor(trt) + factor(tvcon[i,]) + x1 + x2 + as.factor(x3) + 
                  factor(mediator[i-1,]), 
                family=binomial(link = "logit"))
    m <- append(m, fit1$coefficients[[7]])
  }
  a0 <- append(a0, fit1$coefficients[[1]])
  a1 <- append(a1, fit1$coefficients[[2]])
  lm <- append(lm, fit1$coefficients[[3]])
  x_1m <- append(x_1m, fit1$coefficients[[4]])
  x_2m <- append(x_2m, fit1$coefficients[[5]])
  x_3m <- append(x_2m, fit1$coefficients[[6]])
}

b0 <- vector()
g <- vector()
b1 <- vector()
lo <- vector()
x_1o <- vector()
x_2o <- vector()
x_3o <- vector()
o <- vector()

for(i in 1:nm){
  if(i == 1){
    # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1,
    #                         family = poisson(link = "log"),
    #                         na.action=na.omit))
    fit2 <- logbin(factor(outcome[i,]) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                   + x1 + x2 + as.factor(x3), 
                   start = c(log(start.p[i]), rep(0,6)), maxit = 25000,
                   method = "em",
                   na.action=na.omit)
  } else {
    # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,],
    #                         family = poisson(link = "log"),
    #                         na.action=na.omit))
    fit2 <- logbin(factor(outcome[i,]) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                   + x1 + x2 + as.factor(x3) + factor(outcome[i-1,]), 
                   start = c(log(start.p[i]), rep(0,7)), maxit = 25000,
                   method = "em",
                   na.action=na.omit)
    o <- append(o, fit2$coefficients[[8]])
  }
  b0 <- append(b0, fit2$coefficients[[1]])
  g <- append(g, fit2$coefficients[[2]])
  b1 <- append(b1, fit2$coefficients[[3]])
  lo <- append(lo, fit2$coefficients[[4]])
  x_1o <- append(x_1o, fit2$coefficients[[5]])
  x_2o <- append(x_2o, fit2$coefficients[[6]])
  x_3o <- append(x_3o, fit2$coefficients[[7]])
}

results_continuous <- vector("list")

results_continuous$mat_a0 <-  a0
results_continuous$mat_a1 <- a1
results_continuous$mat_lm <- lm
results_continuous$mat_x1m <- x_1m
results_continuous$mat_x2m <- x_2m
results_continuous$mat_x3m <- x_3m
results_continuous$mat_m <- m

results_continuous$mat_b0 <- b0
results_continuous$mat_g <- g
results_continuous$mat_b1 <- b1
results_continuous$mat_lo <- lo
results_continuous$mat_x1o <- x_1o
results_continuous$mat_x2o <- x_2o
results_continuous$mat_x3o <- x_3o
results_continuous$mat_o <- o

save(results_continuous, file = paste("results_varen_medf_binary.rda", sep = ""))
# save(results_continuous, file = paste("results_cNRT_medf_binary.rda", sep = ""))

# save(results_continuous, file = paste("results_varen_medp_binary.rda", sep = ""))
# save(results_continuous, file = paste("results_cNRT_medp_binary.rda", sep = ""))

### Compute (in)direct effects ###
results <- get(load("~/Aim3/results_varen_medf_binary.rda"))
# results <- get(load("~/Aim3/results_cNRT_medf_binary.rda"))
# results <- get(load("~/Aim3/results_varen_medp_binary.rda"))
# results <- get(load("~/Aim3/results_cNRT_medp_binary.rda"))


cov_x1 <- mean(x1, na.rm = TRUE)
cov_x2 <- mean(x2, na.rm = TRUE)

###  Read each coefficient from the results ###
for (i in 1:(length(results))) {
  mat1 <- results[[i]]
  assign(paste0("est_", i), mat1)
}

est_9 <- exp(est_9)

est_16 <- vector(length = length(est_2))
est_17 <- vector(length = length(est_2))
est_18 <- vector(length = length(est_2))
est_19 <- vector(length = length(est_2))
est_20 <- vector(length = length(est_2))

# Gender = 1 #

for (i in 1:length(est_2)) {
  est_16[i] <- (exp(est_2[i])/(1+exp(est_2[i])))*exp(est_10[i])
  
  if(i == 1){
    ## At baseline L(t+1)=1; C=avg cov
    orNIE_l1 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At baseline L(t+1)=0; C=avg cov
    orNIE_l0 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=avg cov
    est_17[i] <- orNIE_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=avg cov
    est_18[i] <- orNIE_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=avg cov
    est_19[i] <- orNIE_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=avg cov
    est_20[i] <- orNIE_l0
    
  }else{
    ## At t+1 M(t)=1; L(t+1)=1; C=avg cov
    orNIE_m1l1 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=avg cov
    orNIE_m0l1 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=avg cov
    orNIE_m1l0 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=avg cov
    orNIE_m0l0 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=avg cov
    est_17[i] <- orNIE_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=avg cov
    est_18[i] <- orNIE_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=avg cov
    est_19[i] <- orNIE_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=avg cov
    est_20[i] <- orNIE_m0l0
  }
}

effect_9 <- est_9
effect_16 <- est_16
effect_17 <- est_17
effect_18 <- est_18
effect_19 <- est_19
effect_20 <- est_20

est_9 <- predict(loess(effect_9 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_16 <- predict(loess(effect_16 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_17 <- predict(loess(effect_17 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_18 <- predict(loess(effect_18 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_19 <- predict(loess(effect_19 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_20 <- predict(loess(effect_20 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_smooth <- vector("list", 20)

for (i in 1:20) {
  est_smooth[[i]] <-  get(paste0("est_", i))
}

## Bootstrap 95% CIs ##

CI <- vector("list", 6)
CI[[1]] <- vector()
CI[[2]] <- vector()
CI[[3]] <- vector()
CI[[4]] <- vector()
CI[[5]] <- vector()
CI[[6]] <- vector()

R = 1000

mat_g <- matrix(nrow = R, ncol = length(t.seq))
mat_ab <- matrix(nrow = R, ncol = length(t.seq))
mat_IE11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE00 <- matrix(nrow = R, ncol = length(t.seq))

for (r in 1:R) {
  print(paste("CI_Simulation_",r))
  
  n <- length(trt)
  nm <- nrow(outcome)
  
  index1 <- sample(1:n, size=n, replace=TRUE)
  
  a0 <- vector()
  a1 <- vector()
  lm <- vector()
  x_1m <- vector()
  x_2m <- vector()
  x_3m <- vector()
  m <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      fit1 <- glm(factor(mediator[i,index1]) ~ factor(trt[index1]) + factor(tvcon[i,index1]) 
                  + x1[index1] + x2[index1] + x3[index1],
                  family=binomial(link = "logit"))
    } else {
      fit1 <- glm(factor(mediator[i,index1]) ~ factor(trt[index1]) + factor(tvcon[i,index1]) +  
                    x1[index1] + x2[index1] + x3[index1] + factor(mediator[i-1,index1]), 
                  family=binomial(link = "logit"))
      m <- append(m, fit1$coefficients[[7]])
    }
    a0 <- append(a0, fit1$coefficients[[1]])
    a1 <- append(a1, fit1$coefficients[[2]])
    lm <- append(lm, fit1$coefficients[[3]])
    x_1m <- append(x_1m, fit1$coefficients[[4]])
    x_2m <- append(x_2m, fit1$coefficients[[5]])
    x_3m <- append(x_2m, fit1$coefficients[[6]])
  }
  
  b0 <- vector()
  g <- vector()
  b1 <- vector()
  lo <- vector()
  x_1o <- vector()
  x_2o <- vector()
  x_3o <- vector()
  o <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      fit2 <- logbin(factor(outcome[i,index1]) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                     + x1[index1] + x2[index1] + x3[index1], 
                     start = c(log(start.p[i]), rep(c(0), 6)),
                     method = "em", maxit = 25000,
                     na.action=na.omit)
    } else {
      fit2 <- logbin(factor(outcome[i,index1]) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                     + x1[index1] + x2[index1] + x3[index1] + factor(outcome[i-1,index1]),
                     start = c(log(start.p[i]), rep(c(0), 7)),
                     method = "em", maxit = 25000,
                     na.action=na.omit)
      o <- append(o, fit2$coefficients[[8]])
    }
    b0 <- append(b0, fit2$coefficients[[1]])
    g <- append(g, fit2$coefficients[[2]])
    b1 <- append(b1, fit2$coefficients[[3]])
    lo <- append(lo, fit2$coefficients[[4]])
    x_1o <- append(x_1o, fit2$coefficients[[5]])
    x_2o <- append(x_2o, fit2$coefficients[[6]])
    x_3o <- append(x_2o, fit2$coefficients[[7]])
  }
  
  mat_g[r,] <- predict(loess(exp(g) ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_ab[r,] <- predict(loess(((exp(a1)/(1+exp(a1)))*exp(b1)) ~ t.seq[1:nm], span = 0.5, degree = 1))
  
  IE11 <- vector()
  IE01 <- vector()
  IE10 <- vector()
  IE00 <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      ## At baseline L(t+1)=1; C=cov_avg
      orNIE_l1 <- ((1 + exp(b1[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At baseline L(t+1)=0; C=cov_avg
      orNIE_l0 <- ((1 + exp(b1[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## IE at M(t)=1; L(t+1)=1; C=cov_avg
      IE11 <- append(IE11, orNIE_l1)
      
      ## IE at M(t)=0; L(t+1)=1; C=cov_avg
      IE01 <- append(IE01, orNIE_l1)
      
      ## IE at M(t)=1; L(t+1)=0; C=cov_avg
      IE10 <- append(IE10, orNIE_l0)
      
      ## IE at M(t)=0; L(t+1)=0; C=cov_avg
      IE00 <- append(IE00, orNIE_l0)
      
    }else{
      ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
      orNIE_m1l1 <- ((1 + exp(b1[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      orNIE_m0l1 <- ((1 + exp(b1[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      orNIE_m1l0 <- ((1 + exp(b1[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      orNIE_m0l0 <- ((1 + exp(b1[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## IE at M(t)=1; L(t+1)=1; C=cov_avg
      IE11 <- append(IE11, orNIE_m1l1)
      
      ## IE at M(t)=0; L(t+1)=1; C=cov_avg
      IE01 <- append(IE01, orNIE_m0l1)
      
      ## IE at M(t)=1; L(t+1)=0; C=cov_avg
      IE10 <- append(IE10, orNIE_m1l0)
      
      ## IE at M(t)=0; L(t+1)=0; C=cov_avg
      IE00 <- append(IE00, orNIE_m0l0)
    }
  }
  
  mat_IE11[r,] <- predict(loess(IE11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE01[r,] <- predict(loess(IE01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE10[r,] <- predict(loess(IE10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE00[r,] <- predict(loess(IE00 ~ t.seq[1:nm], span = 0.5, degree = 1))
}

effect_interest <- vector("list",6)
effect_interest[[1]] <- mat_g
effect_interest[[2]] <- mat_ab
effect_interest[[3]] <- mat_IE11
effect_interest[[4]] <- mat_IE01
effect_interest[[5]] <- mat_IE10
effect_interest[[6]] <- mat_IE00

for(z in 1:length(effect_interest)){
  mat_interest <- effect_interest[[z]]
  
  # if(z == 5 || z == 11){
  #   quantiles <- matrix(NA, nrow=nm-1, ncol=2)
  #   
  #   lower <- 0.025
  #   upper <- 1 - lower
  #   
  #   for(i in 1:nm-1){
  #     quantiles[i,1] <- quantile(mat_interest[,i], c(lower), na.rm=TRUE)
  #     quantiles[i,2] <- quantile(mat_interest[,i], c(upper), na.rm=TRUE)
  #   }
  #   
  # }else{
  quantiles <- matrix(NA, nrow=nm, ncol=2)
  
  lower <- 0.025
  upper <- 1 - lower
  
  for(i in 1:nm){
    quantiles[i,1] <- quantile(mat_interest[,i], c(lower), na.rm=TRUE)
    quantiles[i,2] <- quantile(mat_interest[,i], c(upper), na.rm=TRUE)
  }
  
  # }
  CI[[z]] <- quantiles
}

save(CI, file = paste("CI__varen_medf_binary_male.rda"))
# save(CI, file = paste("CI__cNRT_medf_binary_male.rda"))
# 
# save(CI, file = paste("CI__varen_medp_binary_male.rda"))
# save(CI, file = paste("CI__cNRT_medp_binary_male.rda"))


# Gender = 0 #

for (i in 1:length(est_2)) {
  est_16[i] <- (exp(est_2[i])/(1+exp(est_2[i])))*exp(est_10[i])
  
  if(i == 1){
    ## At baseline L(t+1)=1; C=avg cov
    orNIE_l1 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At baseline L(t+1)=0; C=avg cov
    orNIE_l0 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=avg cov
    est_17[i] <- orNIE_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=avg cov
    est_18[i] <- orNIE_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=avg cov
    est_19[i] <- orNIE_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=avg cov
    est_20[i] <- orNIE_l0
    
  }else{
    ## At t+1 M(t)=1; L(t+1)=1; C=avg cov
    orNIE_m1l1 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=avg cov
    orNIE_m0l1 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=avg cov
    orNIE_m1l0 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=avg cov
    orNIE_m0l0 <- ((1 + exp(est_10[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=avg cov
    est_17[i] <- orNIE_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=avg cov
    est_18[i] <- orNIE_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=avg cov
    est_19[i] <- orNIE_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=avg cov
    est_20[i] <- orNIE_m0l0
  }
}

effect_9 <- est_9
effect_16 <- est_16
effect_17 <- est_17
effect_18 <- est_18
effect_19 <- est_19
effect_20 <- est_20

est_9 <- predict(loess(effect_9 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_16 <- predict(loess(effect_16 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_17 <- predict(loess(effect_17 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_18 <- predict(loess(effect_18 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_19 <- predict(loess(effect_19 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_20 <- predict(loess(effect_20 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_smooth <- vector("list", 20)

for (i in 1:20) {
  est_smooth[[i]] <-  get(paste0("est_", i))
}

## Bootstrap 95% CIs ##

CI <- vector("list", 6)
CI[[1]] <- vector()
CI[[2]] <- vector()
CI[[3]] <- vector()
CI[[4]] <- vector()
CI[[5]] <- vector()
CI[[6]] <- vector()

R = 1000

mat_g <- matrix(nrow = R, ncol = length(t.seq))
mat_ab <- matrix(nrow = R, ncol = length(t.seq))
mat_IE11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE00 <- matrix(nrow = R, ncol = length(t.seq))

for (r in 1:R) {
  print(paste("CI_Simulation_",r))
  
  n <- length(trt)
  nm <- nrow(outcome)
  
  index1 <- sample(1:n, size=n, replace=TRUE)
  
  a0 <- vector()
  a1 <- vector()
  lm <- vector()
  x_1m <- vector()
  x_2m <- vector()
  x_3m <- vector()
  m <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      fit1 <- glm(factor(mediator[i,index1]) ~ factor(trt[index1]) + factor(tvcon[i,index1]) 
                  + x1[index1] + x2[index1] + x3[index1],
                  family=binomial(link = "logit"))
    } else {
      fit1 <- glm(factor(mediator[i,index1]) ~ factor(trt[index1]) + factor(tvcon[i,index1]) +  
                    x1[index1] + x2[index1] + x3[index1] + factor(mediator[i-1,index1]), 
                  family=binomial(link = "logit"))
      m <- append(m, fit1$coefficients[[7]])
    }
    a0 <- append(a0, fit1$coefficients[[1]])
    a1 <- append(a1, fit1$coefficients[[2]])
    lm <- append(lm, fit1$coefficients[[3]])
    x_1m <- append(x_1m, fit1$coefficients[[4]])
    x_2m <- append(x_2m, fit1$coefficients[[5]])
    x_3m <- append(x_2m, fit1$coefficients[[6]])
  }
  
  b0 <- vector()
  g <- vector()
  b1 <- vector()
  lo <- vector()
  x_1o <- vector()
  x_2o <- vector()
  x_3o <- vector()
  o <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      fit2 <- logbin(factor(outcome[i,index1]) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                     + x1[index1] + x2[index1] + x3[index1], 
                     start = c(log(start.p[i]), rep(c(0), 6)),
                     method = "em", maxit = 25000,
                     na.action=na.omit)
    } else {
      fit2 <- logbin(factor(outcome[i,index1]) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                     + x1[index1] + x2[index1] + x3[index1] + factor(outcome[i-1,index1]),
                     start = c(log(start.p[i]), rep(c(0), 7)),
                     method = "em", maxit = 25000,
                     na.action=na.omit)
      o <- append(o, fit2$coefficients[[8]])
    }
    b0 <- append(b0, fit2$coefficients[[1]])
    g <- append(g, fit2$coefficients[[2]])
    b1 <- append(b1, fit2$coefficients[[3]])
    lo <- append(lo, fit2$coefficients[[4]])
    x_1o <- append(x_1o, fit2$coefficients[[5]])
    x_2o <- append(x_2o, fit2$coefficients[[6]])
    x_3o <- append(x_2o, fit2$coefficients[[7]])
  }
  
  mat_g[r,] <- predict(loess(exp(g) ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_ab[r,] <- predict(loess(((exp(a1)/(1+exp(a1)))*exp(b1)) ~ t.seq[1:nm], span = 0.5, degree = 1))
  
  IE11 <- vector()
  IE01 <- vector()
  IE10 <- vector()
  IE00 <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      ## At baseline L(t+1)=1; C=cov_avg
      orNIE_l1 <- ((1 + exp(b1[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At baseline L(t+1)=0; C=cov_avg
      orNIE_l0 <- ((1 + exp(b1[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## IE at M(t)=1; L(t+1)=1; C=cov_avg
      IE11 <- append(IE11, orNIE_l1)
      
      ## IE at M(t)=0; L(t+1)=1; C=cov_avg
      IE01 <- append(IE01, orNIE_l1)
      
      ## IE at M(t)=1; L(t+1)=0; C=cov_avg
      IE10 <- append(IE10, orNIE_l0)
      
      ## IE at M(t)=0; L(t+1)=0; C=cov_avg
      IE00 <- append(IE00, orNIE_l0)
      
    }else{
      ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
      orNIE_m1l1 <- ((1 + exp(b1[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      orNIE_m0l1 <- ((1 + exp(b1[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      orNIE_m1l0 <- ((1 + exp(b1[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      orNIE_m0l0 <- ((1 + exp(b1[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## IE at M(t)=1; L(t+1)=1; C=cov_avg
      IE11 <- append(IE11, orNIE_m1l1)
      
      ## IE at M(t)=0; L(t+1)=1; C=cov_avg
      IE01 <- append(IE01, orNIE_m0l1)
      
      ## IE at M(t)=1; L(t+1)=0; C=cov_avg
      IE10 <- append(IE10, orNIE_m1l0)
      
      ## IE at M(t)=0; L(t+1)=0; C=cov_avg
      IE00 <- append(IE00, orNIE_m0l0)
    }
  }
  
  mat_IE11[r,] <- predict(loess(IE11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE01[r,] <- predict(loess(IE01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE10[r,] <- predict(loess(IE10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE00[r,] <- predict(loess(IE00 ~ t.seq[1:nm], span = 0.5, degree = 1))
}

effect_interest <- vector("list",6)
effect_interest[[1]] <- mat_g
effect_interest[[2]] <- mat_ab
effect_interest[[3]] <- mat_IE11
effect_interest[[4]] <- mat_IE01
effect_interest[[5]] <- mat_IE10
effect_interest[[6]] <- mat_IE00

for(z in 1:length(effect_interest)){
  mat_interest <- effect_interest[[z]]
  
  # if(z == 5 || z == 11){
  #   quantiles <- matrix(NA, nrow=nm-1, ncol=2)
  #   
  #   lower <- 0.025
  #   upper <- 1 - lower
  #   
  #   for(i in 1:nm-1){
  #     quantiles[i,1] <- quantile(mat_interest[,i], c(lower), na.rm=TRUE)
  #     quantiles[i,2] <- quantile(mat_interest[,i], c(upper), na.rm=TRUE)
  #   }
  #   
  # }else{
  quantiles <- matrix(NA, nrow=nm, ncol=2)
  
  lower <- 0.025
  upper <- 1 - lower
  
  for(i in 1:nm){
    quantiles[i,1] <- quantile(mat_interest[,i], c(lower), na.rm=TRUE)
    quantiles[i,2] <- quantile(mat_interest[,i], c(upper), na.rm=TRUE)
  }
  
  # }
  CI[[z]] <- quantiles
}

save(CI, file = paste("CI__varen_medf_binary_female.rda"))
# save(CI, file = paste("CI__cNRT_medf_binary_female.rda"))
# 
# save(CI, file = paste("CI__varen_medp_binary_female.rda"))
# save(CI, file = paste("CI__cNRT_medp_binary_female.rda"))


# ### Plot the (In)direct effects curves ###
# effect_interest <- vector("list",6)
# effect_interest[[1]] <- est_smooth[[7]]
# effect_interest[[2]] <- est_smooth[[12]]
# 
# effect_interest[[3]] <- est_smooth[[13]]
# effect_interest[[4]] <- est_smooth[[14]]
# effect_interest[[5]] <- est_smooth[[15]]
# effect_interest[[6]] <- est_smooth[[16]]
# 
# t.est <- t.seq
# 
# dat_g <- as.data.frame(cbind(effect_interest[[1]], CI[[1]], t.est))
# dat_ab <- as.data.frame(cbind(effect_interest[[2]], CI[[2]], t.est))
# dat_IE11 <- as.data.frame(cbind(effect_interest[[3]], CI[[3]]))
# dat_IE01 <- as.data.frame(cbind(effect_interest[[4]], CI[[4]]))
# dat_IE10 <- as.data.frame(cbind(effect_interest[[5]], CI[[5]]))
# dat_IE00 <- as.data.frame(cbind(effect_interest[[6]], CI[[6]]))
# dat_POF <- as.data.frame(cbind(dat_IE11, dat_IE01, dat_IE10, dat_IE00, t.est))
# 
# var_names <- c("IE11", "LL_IE11", "UL_IE11",
#                "IE01", "LL_IE01", "UL_IE01",
#                "IE10", "LL_IE10", "UL_IE10",
#                "IE00", "LL_IE00", "UL_IE00", "t.est")
# 
# colnames(dat_POF) <- var_names
# 
# ggplot(data = dat_g, aes(t.est, V1)) +
#   geom_line(color = "blueviolet", linewidth = 0.8) +
#   geom_ribbon(data=dat_g, aes(ymin=V2, ymax=V3), fill="lightcyan3", alpha=0.5) +
#   geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
#   labs(title = "Plot of the time-varying direct effect - γ(t)",
#        x = "Time Sequence",
#        y = "γ(t) Curve") +
#   scale_x_continuous(breaks = t.est)
# 
# ggplot(data = dat_ab, aes(t.est, V1)) +
#   geom_line(color = "darkgreen", linewidth = 0.75) +
#   geom_ribbon(data=dat_ab, aes(ymin=V2, ymax=V3), fill="azure4", alpha=0.5) +
#   geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
#   labs(title = "Plot of the time-varying mediation effect - POC",
#        x = "Time Sequence",
#        y = "IE POC Curve") +
#   scale_x_continuous(breaks = t.est)
# 
# plot_POF1 <- ggplot(data = dat_POF, aes(t.est, IE11)) +
#   geom_line(color = "red", linewidth = 0.75) +
#   geom_ribbon(data=dat_POF, aes(ymin=LL_IE11, ymax=UL_IE11), fill="lightpink", alpha=0.5) +
#   
#   geom_line(aes(t.est, IE01), color = "blue", linewidth = 0.75) +
#   geom_ribbon(data=dat_POF, aes(ymin=LL_IE01, ymax=UL_IE01), fill="lightblue", alpha=0.5) + 
#   
#   geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
#   
#   labs(title = "Plot of the time-varying mediation effect - POF",
#        x = "Time Sequence",
#        y = "Indirect Effect Curves") +
#   scale_x_continuous(breaks = t.est)
# 
# plot_POF2 <- ggplot(data = dat_POF, aes(t.est, IE10)) +
#   geom_line(color = "darkred", linewidth = 0.75) +
#   geom_ribbon(data=dat_POF, aes(ymin=LL_IE10, ymax=UL_IE10), fill="yellow", alpha=0.5) + 
#   
#   
#   geom_line(aes(t.est, IE00), color = "darkgreen", linewidth = 0.75) +
#   geom_ribbon(data=dat_POF, aes(ymin=LL_IE00, ymax=UL_IE00), fill="lightgreen", alpha=0.5) + 
#   
#   geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
#   
#   labs(title = "Plot of the time-varying mediation effect - POF",
#        x = "Time Sequence",
#        y = "Indirect Effect Curves") +
#   scale_x_continuous(breaks = t.est)



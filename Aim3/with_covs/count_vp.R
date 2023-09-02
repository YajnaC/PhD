library(dplyr)
library(reshape2)
library(ggplot2)
library(tvmediation)
library(logbin)
library(brglm)
library(pscl)

### Complete data analysis ###
setwd("~/Aim3/with_covs")
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

mean.y <- as.data.frame(table(data$DaysFromTQD, data$CigCount_Today)) %>%
  filter(Var2 == 0)
# logit_pr <- (nrow(ID)-mean.y$Freq)/nrow(ID)
logit_pr <- (mean.y$Freq)/length(trt)
# start.p <- logit_pr/(1-logit_pr)
start.p <- logit_pr

### Full Adherence - Mediator ###
# mediator <- LongToWide(data$SubjectID, data$DaysFromTQD, data$adh_full)

### Partial Adherence - Mediator ###
mediator <- LongToWide(data$SubjectID, data$DaysFromTQD, data$adh_part)

tvcon <- LongToWide(data$SubjectID, data$DaysFromTQD, data$StressEventToday)

### Count Outcome - Daily Number of Cigarettes Smoked ###
outcome <- LongToWide(data$SubjectID, data$DaysFromTQD, data$CigCount_Today)

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

b0_count <- vector()
g_count <- vector()
b1_count <- vector()
lo_count <- vector()
x_1o_count <- vector()
x_2o_count <- vector()
x_3o_count <- vector()
# o_count <- vector()

b0_zero <- vector()
g_zero <- vector()
b1_zero <- vector()
lo_zero <- vector()
x_1o_zero <- vector()
x_2o_zero <- vector()
x_3o_zero <- vector()
o_zero <- vector()

for(i in 1:nm){
  if(i == 1){
    start.coeff <- glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                       + x1 + x2 + as.factor(x3),
                       family = poisson(link = "log"),
                       na.action=na.omit)
    start.c.p <- exp(coef(start.coeff)[1])
    # x <- model.matrix(start.coeff)
    
    # start.coeff1 <- logbin(factor(outcome[i,]==0) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1,
    #                        # family = binomial(link = "log"),
    #                        start = c(log(start.p[i]), rep(0,4)), maxit = 25000,
    #                        method = "em",
    #                        na.action=na.omit)
    # start.z.p <- 1 - exp(coef(start.coeff1)[1])
    
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
    
    # fit2 <- zeroinfl(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 | factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1,
    #                  na.action=na.omit, method = "Nelder-Mead",  link = "log",
    #                  start = list(count = c(log(start.c.p), rep(0, 4)), 
    #                               zero = c(log(start.z.p), rep(0, 4))))
    
    fit2 <- zeroinfl(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                     + x1 + x2 + as.factor(x3) | factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                     + x1 + x2 + as.factor(x3),
                     na.action=na.omit, method = "Nelder-Mead",  link = "log",
                     start = list(count = c(log(start.c.p), rep(0, 6)), 
                                  zero = c(log(start.p[i]), rep(0, 6))))
  } else {
    # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,] ,
    #                         family = poisson(link = "log"),
    #                         na.action=na.omit))
    # start.c.p <- exp(start.coeff)/(1+exp(start.coeff))
    
    start.coeff <- glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                       + x1 + x2 + as.factor(x3),
                       family = poisson(link = "log"),
                       na.action=na.omit)
    start.c.p <- exp(coef(start.coeff)[1])
    # x <- model.matrix(start.coeff)
    
    # start.coeff1 <- logbin(factor(outcome[i,]==0) ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,],
    #                        # family = binomial(link = "log"),
    #                        start = c(log(start.p[i]), rep(0,5)), maxit = 25000,
    #                        method = "em",
    #                        na.action=na.omit)
    # start.z.p <- 1 - exp(coef(start.coeff1)[1])
    
    # fit2 <- zeroinfl(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,] | factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,],
    #                  na.action=na.omit, method = "Nelder-Mead", link = "log",
    #                  start = list(count = c(log(start.c.p), rep(0, 5)), 
    #                               zero = c(log(start.z.p), rep(0, 5))))
    
    fit2 <- zeroinfl(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                     + x1 + x2 + as.factor(x3) | factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) 
                     + x1 + x2 + as.factor(x3) + outcome[i-1,],
                     na.action=na.omit, method = "Nelder-Mead", link = "log",
                     start = list(count = c(log(start.c.p), rep(0, 6)), 
                                  zero = c(log(start.p[i]), rep(0, 7))))
    
    # o_count <- append(o_count, fit2$coefficients$count[[8]])
    o_zero <- append(o_zero, fit2$coefficients$zero[[8]])
  }
  b0_count <- append(b0_count, fit2$coefficients$count[[1]])
  g_count <- append(g_count, fit2$coefficients$count[[2]])
  b1_count <- append(b1_count, fit2$coefficients$count[[3]])
  lo_count <- append(lo_count, fit2$coefficients$count[[4]])
  x_1o_count <- append(x_1o_count, fit2$coefficients$count[[5]])
  x_2o_count <- append(x_2o_count, fit2$coefficients$count[[6]])
  x_3o_count <- append(x_3o_count, fit2$coefficients$count[[7]])
  
  b0_zero <- append(b0_zero, fit2$coefficients$zero[[1]])
  g_zero <- append(g_zero, fit2$coefficients$zero[[2]])
  b1_zero <- append(b1_zero, fit2$coefficients$zero[[3]])
  lo_zero <- append(lo_zero, fit2$coefficients$zero[[4]])
  x_1o_zero <- append(x_1o_zero, fit2$coefficients$zero[[5]])
  x_2o_zero <- append(x_2o_zero, fit2$coefficients$zero[[6]])
  x_3o_zero <- append(x_3o_zero, fit2$coefficients$zero[[7]])
}

results_continuous <- vector("list")

results_continuous$mat_a0 <-  a0
results_continuous$mat_a1 <- a1
results_continuous$mat_lm <- lm
results_continuous$mat_x1m <- x_1m
results_continuous$mat_x2m <- x_2m
results_continuous$mat_x3m <- x_3m
results_continuous$mat_m <- m

results_continuous$mat_b0_count <- b0_count
results_continuous$mat_g_count <- g_count
results_continuous$mat_b1_count <- b1_count
results_continuous$mat_lo_count <- lo_count
results_continuous$mat_x1o_count <- x_1o_count
results_continuous$mat_x2o_count <- x_2o_count
results_continuous$mat_x3o_count <- x_3o_count
# results_continuous$mat_o_count <- o_count

results_continuous$mat_b0_zero <- b0_zero
results_continuous$mat_g_zero <- g_zero
results_continuous$mat_b1_zero <- b1_zero
results_continuous$mat_lo_zero <- lo_zero
results_continuous$mat_x1o_zero <- x_1o_zero
results_continuous$mat_x2o_zero <- x_2o_zero
results_continuous$mat_x3o_zero <- x_3o_zero
results_continuous$mat_o_zero <- o_zero

# save(results_continuous, file = paste("results_varen_medf_count.rda", sep = ""))
# save(results_continuous, file = paste("results_cNRT_medf_count.rda", sep = ""))

save(results_continuous, file = paste("results_varen_medp_count.rda", sep = ""))
# save(results_continuous, file = paste("results_cNRT_medp_count.rda", sep = ""))

### Compute (in)direct effects ###
# results <- get(load("~/Aim3/with_covs/results_varen_medf_count.rda"))
# results <- get(load("~/Aim3/with_covs/results_cNRT_medf_count.rda"))
# 
results <- get(load("~/Aim3/with_covs/results_varen_medp_count.rda"))
# results <- get(load("~/Aim3/with_covs/results_cNRT_medp_count.rda"))

cov_x1 <- mean(x1, na.rm = TRUE)
cov_x2 <- mean(x2, na.rm = TRUE)

###  Read each coefficient from the results ###
for (i in 1:(length(results))) {
  mat1 <- results[[i]]
  assign(paste0("est_", i), mat1)
}

# Gender = 1 #

est_9 <- exp(est_9)

est_23 <- vector(length = length(est_2)) # this is just to keep the numbers same

est_24 <- vector(length = length(est_2))
est_25 <- vector(length = length(est_2))
est_26 <- vector(length = length(est_2))
est_27 <- vector(length = length(est_2))
est_28 <- vector(length = length(est_2))
est_29 <- vector(length = length(est_2))
est_30 <- vector(length = length(est_2))
est_31 <- vector(length = length(est_2))
est_32 <- vector(length = length(est_2))
est_33 <- vector(length = length(est_2))
est_34 <- vector(length = length(est_2))
est_35 <- vector(length = length(est_2))
est_36 <- vector(length = length(est_2))

for (i in 1:length(est_2)) {
  est_24[i] <- (exp(est_2[i])/(1+exp(est_2[i])))*exp(est_10[i])
  
  if(i == 1){
    ### Overall Indirect Effect ###
    ## At baseline L(t+1)=1; C=cov_avg
    orNIE_l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At baseline L(t+1)=0; C=cov_avg
    orNIE_l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_25[i] <- orNIE_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_26[i] <- orNIE_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_27[i] <- orNIE_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_28[i] <- orNIE_l0
    
    
    ### Indirect Effect through susceptibility latent variable ###
    ## At baseline L(t+1)=1; C=cov_avg
    rrIES_l1 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At baseline L(t+1)=0; C=cov_avg
    rrIES_l0 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_29[i] <- rrIES_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_30[i] <- rrIES_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_31[i] <- rrIES_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_32[i] <- rrIES_l0
    
    ### Indirect Effect not through susceptibility latent variable ###
    ## At baseline L(t+1)=1; C=cov_avg
    rrIENS_l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At baseline L(t+1)=0; C=cov_avg
    rrIENS_l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_33[i] <- rrIENS_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_34[i] <- rrIENS_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_35[i] <- rrIENS_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_36[i] <- rrIENS_l0
    
  }else{
    ### Overall Indirect Effect ###
    ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
    orNIE_m1l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
    orNIE_m0l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
    orNIE_m1l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
    orNIE_m0l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_25[i] <- orNIE_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_26[i] <- orNIE_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_27[i] <- orNIE_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_28[i] <- orNIE_m0l0
    
    
    ### Indirect Effect through susceptibility latent variable ###
    ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
    rrIES_m1l1 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
    rrIES_m0l1 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
    rrIES_m1l0 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
    rrIES_m0l0 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_29[i] <- rrIES_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_30[i] <- rrIES_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_31[i] <- rrIES_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_32[i] <- rrIES_m0l0
    
    ### Indirect Effect not through susceptibility latent variable ###
    ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
    rrIENS_m1l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
    rrIENS_m0l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
    rrIENS_m1l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i] + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
    rrIENS_m0l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_6[i])))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_33[i] <- rrIENS_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_34[i] <- rrIENS_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_35[i] <- rrIENS_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_36[i] <- rrIENS_m0l0
  }
}

effect_9 <- est_9
effect_24 <- est_24

effect_25 <- est_25
effect_26 <- est_26
effect_27 <- est_27
effect_28 <- est_28

effect_29 <- est_29
effect_30 <- est_30
effect_31 <- est_31
effect_32 <- est_32

effect_33 <- est_33
effect_34 <- est_34
effect_35 <- est_35
effect_36 <- est_36

est_9 <- predict(loess(effect_9 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_24 <- predict(loess(effect_24 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_25 <- predict(loess(effect_25 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_26 <- predict(loess(effect_26 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_27 <- predict(loess(effect_27 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_28 <- predict(loess(effect_28 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_29 <- predict(loess(effect_29 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_30 <- predict(loess(effect_30 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_31 <- predict(loess(effect_31 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_32 <- predict(loess(effect_32 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_33 <- predict(loess(effect_33 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_34 <- predict(loess(effect_34 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_35 <- predict(loess(effect_35 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_36 <- predict(loess(effect_36 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_smooth <- vector("list", 36)

for (i in 1:36) {
  est_smooth[[i]] <-  get(paste0("est_", i))
}

## Bootstrap 95% CIs ##

CI <- vector("list", 14)
CI[[1]] <- vector()
CI[[2]] <- vector()
CI[[3]] <- vector()
CI[[4]] <- vector()
CI[[5]] <- vector()
CI[[6]] <- vector()
CI[[7]] <- vector()
CI[[8]] <- vector()
CI[[9]] <- vector()
CI[[10]] <- vector()
CI[[11]] <- vector()
CI[[12]] <- vector()
CI[[13]] <- vector()
CI[[14]] <- vector()

R = 1000

mat_g <- matrix(nrow = R, ncol = length(t.seq))
mat_ab <- matrix(nrow = R, ncol = length(t.seq))

mat_IE11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE00 <- matrix(nrow = R, ncol = length(t.seq))

mat_IES11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IES01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IES10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IES00 <- matrix(nrow = R, ncol = length(t.seq))

mat_IENS11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IENS01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IENS10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IENS00 <- matrix(nrow = R, ncol = length(t.seq))


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
  
  b0_count <- vector()
  g_count <- vector()
  b1_count <- vector()
  lo_count <- vector()
  x_1o_count <- vector()
  x_2o_count <- vector()
  x_3o_count <- vector()
  o_count <- vector()
  
  b0_zero <- vector()
  g_zero <- vector()
  b1_zero <- vector()
  lo_zero <- vector()
  x_1o_zero <- vector()
  x_2o_zero <- vector()
  x_3o_zero <- vector()
  o_zero <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      start.coeff <- glm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                         + x1[index1] + x2[index1] + x3[index1],
                         family = poisson(link = "log"),
                         na.action=na.omit)
      start.c.p <- exp(coef(start.coeff)[1])
      # x <- model.matrix(start.coeff)
      
      # start.coeff1 <- logbin(factor(outcome[i,index1]==0) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1],
      #                        # family = binomial(link = "log"),
      #                        start = c(log(start.p[i]), rep(0,4)), maxit = 25000,
      #                        method = "em",
      #                        na.action=na.omit)
      # start.z.p <- 1 - exp(coef(start.coeff1)[1])
      
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
      
      
      # fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1],
      #                  na.action=na.omit, method = "Nelder-Mead",  link = "log",
      #                  start = list(count = c(log(start.c.p), rep(0, 4)), 
      #                               zero = c(log(start.z.p), rep(0, 4))))
      
      fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1],
                       na.action=na.omit, method = "Nelder-Mead",  link = "log",
                       start = list(count = c(log(start.c.p), rep(0, 6)), 
                                    zero = c(log(start.p[i]), rep(0, 6))))
    } else {
      # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,] ,
      #                         family = poisson(link = "log"),
      #                         na.action=na.omit))
      # start.c.p <- exp(start.coeff)/(1+exp(start.coeff))
      
      start.coeff <- glm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                         + x1[index1] + x2[index1] + x3[index1],
                         family = poisson(link = "log"),
                         na.action=na.omit)
      start.c.p <- exp(coef(start.coeff)[1])
      # x <- model.matrix(start.coeff)
      
      # start.coeff1 <- logbin(factor(outcome[i,index1]==0) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1],
      #                        # family = binomial(link = "log"),
      #                        start = c(log(start.p[i]), rep(0,5)), maxit = 25000,
      #                        method = "em",
      #                        na.action=na.omit)
      # start.z.p <- 1 - exp(coef(start.coeff1)[1])
      
      # fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1],
      #                  na.action=na.omit, method = "Nelder-Mead", link = "log",
      #                  start = list(count = c(log(start.c.p), rep(0, 5)), 
      #                               zero = c(log(start.z.p), rep(0, 5))))
      
      fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1] + outcome[i-1,index1],
                       na.action=na.omit, method = "Nelder-Mead", link = "log",
                       start = list(count = c(log(start.c.p), rep(0, 6)), 
                                    zero = c(log(start.p[i]), rep(0, 7))))
      
      # o_count <- append(o_count, fit2$coefficients$count[[8]])
      o_zero <- append(o_zero, fit2$coefficients$zero[[8]])
    }
    b0_count <- append(b0_count, fit2$coefficients$count[[1]])
    g_count <- append(g_count, fit2$coefficients$count[[2]])
    b1_count <- append(b1_count, fit2$coefficients$count[[3]])
    lo_count <- append(lo_count, fit2$coefficients$count[[4]])
    x_1o_count <- append(x_1o_count, fit2$coefficients$count[[5]])
    x_2o_count <- append(x_2o_count, fit2$coefficients$count[[6]])
    x_3o_count <- append(x_3o_count, fit2$coefficients$count[[7]])
    
    
    b0_zero <- append(b0_zero, fit2$coefficients$zero[[1]])
    g_zero <- append(g_zero, fit2$coefficients$zero[[2]])
    b1_zero <- append(b1_zero, fit2$coefficients$zero[[3]])
    lo_zero <- append(lo_zero, fit2$coefficients$zero[[4]])
    x_1o_zero <- append(x_1o_zero, fit2$coefficients$zero[[5]])
    x_2o_zero <- append(x_2o_zero, fit2$coefficients$zero[[6]])
    x_3o_zero <- append(x_3o_zero, fit2$coefficients$zero[[7]])
  }
  
  # mat_b0[r,] <- b0_count
  mat_g[r,] <- predict(loess(exp(g_count) ~ t.seq[1:nm], span = 0.5, degree = 1))
  # mat_b1[r,] <- predict(loess(b1_count ~ t.seq[1:nm], span = 0.5, degree = 1))
  # mat_lo[r,] <- lo_count
  # mat_x1o[r,] <- x_1o_count
  # # mat_x2o[r,] <- x_2o_count
  # mat_o[r,] <- o_count
  # 
  # mat_b0_zero[r,] <- b0_zero
  # mat_g_zero[r,] <- g_zero
  # mat_b1_zero[r,] <- b1_zero
  # mat_lo_zero[r,] <- lo_zero
  # mat_x1o_zero[r,] <- x_1o_zero
  # # mat_x2o_zero[r,] <- x_2o_zero
  # mat_o_zero[r,] <- o_zero
  mat_ab[r,] <- predict(loess(((exp(a1)/(1+exp(a1)))*exp(b1_count)) ~ t.seq[1:nm], span = 0.5, degree = 1))
  
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
  
  for(i in 1:nm){
    if(i == 1){
      ### Overall Indirect Effect ###
      ## At baseline L(t+1)=1; C=cov_avg
      orNIE_l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At baseline L(t+1)=0; C=cov_avg
      orNIE_l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
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
      rrIES_l1 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At baseline L(t+1)=0; C=cov_avg
      rrIES_l0 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
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
      rrIENS_l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At baseline L(t+1)=0; C=cov_avg
      rrIENS_l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
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
      orNIE_m1l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      orNIE_m0l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      orNIE_m1l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      orNIE_m0l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
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
      rrIES_m1l1 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      rrIES_m0l1 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      rrIES_m1l0 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      rrIES_m0l0 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
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
      rrIENS_m1l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      rrIENS_m0l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      rrIENS_m1l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i] + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      rrIENS_m0l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + x_3m[i])))
      
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
  
  mat_IE11[r,] <- predict(loess(IE11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE01[r,] <- predict(loess(IE01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE10[r,] <- predict(loess(IE10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE00[r,] <- predict(loess(IE00 ~ t.seq[1:nm], span = 0.5, degree = 1))
  
  mat_IES11[r,] <- predict(loess(IES11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IES01[r,] <- predict(loess(IES01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IES10[r,] <- predict(loess(IES10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IES00[r,] <- predict(loess(IES00 ~ t.seq[1:nm], span = 0.5, degree = 1))
  
  mat_IENS11[r,] <- predict(loess(IENS11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IENS01[r,] <- predict(loess(IENS01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IENS10[r,] <- predict(loess(IENS10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IENS00[r,] <- predict(loess(IENS00 ~ t.seq[1:nm], span = 0.5, degree = 1))
}

effect_interest <- vector("list",14)

effect_interest[[1]] <- mat_g
effect_interest[[2]] <- mat_ab
effect_interest[[3]] <- mat_IE11
effect_interest[[4]] <- mat_IE01
effect_interest[[5]] <- mat_IE10
effect_interest[[6]] <- mat_IE00

effect_interest[[7]] <- mat_IES11
effect_interest[[8]] <- mat_IES01
effect_interest[[9]] <- mat_IES10
effect_interest[[10]] <- mat_IES00

effect_interest[[11]] <- mat_IENS11
effect_interest[[12]] <- mat_IENS01
effect_interest[[13]] <- mat_IENS10
effect_interest[[14]] <- mat_IENS00

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

# save(CI, file = paste("CI__varen_medf_count_male.rda"))
# save(CI, file = paste("CI__cNRT_medf_count_male.rda"))
# 
save(CI, file = paste("CI__varen_medp_count_male.rda"))
# save(CI, file = paste("CI__cNRT_medp_count_male.rda"))

# Gender = 0 #

est_9 <- exp(est_9)

est_23 <- vector(length = length(est_2))

est_24 <- vector(length = length(est_2))
est_25 <- vector(length = length(est_2))
est_26 <- vector(length = length(est_2))
est_27 <- vector(length = length(est_2))
est_28 <- vector(length = length(est_2))
est_29 <- vector(length = length(est_2))
est_30 <- vector(length = length(est_2))
est_31 <- vector(length = length(est_2))
est_32 <- vector(length = length(est_2))
est_33 <- vector(length = length(est_2))
est_34 <- vector(length = length(est_2))
est_35 <- vector(length = length(est_2))
est_36 <- vector(length = length(est_2))

for (i in 1:length(est_2)) {
  est_24[i] <- (exp(est_2[i])/(1+exp(est_2[i])))*exp(est_10[i])
  
  if(i == 1){
    ### Overall Indirect Effect ###
    ## At baseline L(t+1)=1; C=cov_avg
    orNIE_l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At baseline L(t+1)=0; C=cov_avg
    orNIE_l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_25[i] <- orNIE_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_26[i] <- orNIE_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_27[i] <- orNIE_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_28[i] <- orNIE_l0
    
    
    ### Indirect Effect through susceptibility latent variable ###
    ## At baseline L(t+1)=1; C=cov_avg
    rrIES_l1 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At baseline L(t+1)=0; C=cov_avg
    rrIES_l0 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_29[i] <- rrIES_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_30[i] <- rrIES_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_31[i] <- rrIES_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_32[i] <- rrIES_l0
    
    ### Indirect Effect not through susceptibility latent variable ###
    ## At baseline L(t+1)=1; C=cov_avg
    rrIENS_l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At baseline L(t+1)=0; C=cov_avg
    rrIENS_l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_33[i] <- rrIENS_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_34[i] <- rrIENS_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_35[i] <- rrIENS_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_36[i] <- rrIENS_l0
    
  }else{
    ### Overall Indirect Effect ###
    ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
    orNIE_m1l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
    orNIE_m0l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
    orNIE_m1l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
    orNIE_m0l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_25[i] <- orNIE_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_26[i] <- orNIE_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_27[i] <- orNIE_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_28[i] <- orNIE_m0l0
    
    
    ### Indirect Effect through susceptibility latent variable ###
    ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
    rrIES_m1l1 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
    rrIES_m0l1 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
    rrIES_m1l0 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
    rrIES_m0l0 <- ((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_17[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_29[i] <- rrIES_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_30[i] <- rrIES_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_31[i] <- rrIES_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_32[i] <- rrIES_m0l0
    
    ### Indirect Effect not through susceptibility latent variable ###
    ## At t+1 M(t)=1; L(t+1)=1; C=cov_avg
    rrIENS_m1l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
    rrIENS_m0l1 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_3[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
    rrIENS_m1l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1]))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2 + est_7[i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
    rrIENS_m0l0 <- ((1 + exp(est_10[i] + est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))/((1 + exp(est_17[i] + est_1[i] + est_2[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2))*(1 + exp(est_10[i] + est_1[i] + est_4[i]*cov_x1 + est_5[i]*cov_x2)))
    
    ## IE at M(t)=1; L(t+1)=1; C=cov_avg
    est_33[i] <- rrIENS_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=cov_avg
    est_34[i] <- rrIENS_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=cov_avg
    est_35[i] <- rrIENS_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=cov_avg
    est_36[i] <- rrIENS_m0l0
  }
}

effect_9 <- est_9
effect_24 <- est_24

effect_25 <- est_25
effect_26 <- est_26
effect_27 <- est_27
effect_28 <- est_28

effect_29 <- est_29
effect_30 <- est_30
effect_31 <- est_31
effect_32 <- est_32

effect_33 <- est_33
effect_34 <- est_34
effect_35 <- est_35
effect_36 <- est_36

est_9 <- predict(loess(effect_9 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_24 <- predict(loess(effect_24 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_25 <- predict(loess(effect_25 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_26 <- predict(loess(effect_26 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_27 <- predict(loess(effect_27 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_28 <- predict(loess(effect_28 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_29 <- predict(loess(effect_29 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_30 <- predict(loess(effect_30 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_31 <- predict(loess(effect_31 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_32 <- predict(loess(effect_32 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_33 <- predict(loess(effect_33 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_34 <- predict(loess(effect_34 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_35 <- predict(loess(effect_35 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
est_36 <- predict(loess(effect_36 ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))

est_smooth <- vector("list", 36)

for (i in 1:36) {
  est_smooth[[i]] <-  get(paste0("est_", i))
}

## Bootstrap 95% CIs ##

CI <- vector("list", 14)
CI[[1]] <- vector()
CI[[2]] <- vector()
CI[[3]] <- vector()
CI[[4]] <- vector()
CI[[5]] <- vector()
CI[[6]] <- vector()
CI[[7]] <- vector()
CI[[8]] <- vector()
CI[[9]] <- vector()
CI[[10]] <- vector()
CI[[11]] <- vector()
CI[[12]] <- vector()
CI[[13]] <- vector()
CI[[14]] <- vector()

R = 1000

mat_g <- matrix(nrow = R, ncol = length(t.seq))
mat_ab <- matrix(nrow = R, ncol = length(t.seq))

mat_IE11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IE00 <- matrix(nrow = R, ncol = length(t.seq))

mat_IES11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IES01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IES10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IES00 <- matrix(nrow = R, ncol = length(t.seq))

mat_IENS11 <- matrix(nrow = R, ncol = length(t.seq))
mat_IENS01 <- matrix(nrow = R, ncol = length(t.seq))
mat_IENS10 <- matrix(nrow = R, ncol = length(t.seq))
mat_IENS00 <- matrix(nrow = R, ncol = length(t.seq))


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
  
  b0_count <- vector()
  g_count <- vector()
  b1_count <- vector()
  lo_count <- vector()
  x_1o_count <- vector()
  x_2o_count <- vector()
  x_3o_count <- vector()
  # o_count <- vector()
  
  b0_zero <- vector()
  g_zero <- vector()
  b1_zero <- vector()
  lo_zero <- vector()
  x_1o_zero <- vector()
  x_2o_zero <- vector()
  x_3o_zero <- vector()
  o_zero <- vector()
  
  for(i in 1:nm){
    if(i == 1){
      start.coeff <- glm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                         + x1[index1] + x2[index1] + x3[index1],
                         family = poisson(link = "log"),
                         na.action=na.omit)
      start.c.p <- exp(coef(start.coeff)[1])
      # x <- model.matrix(start.coeff)
      
      # start.coeff1 <- logbin(factor(outcome[i,index1]==0) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1],
      #                        # family = binomial(link = "log"),
      #                        start = c(log(start.p[i]), rep(0,4)), maxit = 25000,
      #                        method = "em",
      #                        na.action=na.omit)
      # start.z.p <- 1 - exp(coef(start.coeff1)[1])
      
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
      
      
      # fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1],
      #                  na.action=na.omit, method = "Nelder-Mead",  link = "log",
      #                  start = list(count = c(log(start.c.p), rep(0, 4)), 
      #                               zero = c(log(start.z.p), rep(0, 4))))
      
      fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1],
                       na.action=na.omit, method = "Nelder-Mead",  link = "log",
                       start = list(count = c(log(start.c.p), rep(0, 6)), 
                                    zero = c(log(start.p[i]), rep(0, 6))))
    } else {
      # start.coeff <- coef(glm(outcome[i,] ~ factor(trt) + factor(mediator[i,]) + factor(tvcon[i,]) + x1 + outcome[i-1,] ,
      #                         family = poisson(link = "log"),
      #                         na.action=na.omit))
      # start.c.p <- exp(start.coeff)/(1+exp(start.coeff))
      
      start.coeff <- glm(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                         + x1[index1] + x2[index1] + x3[index1],
                         family = poisson(link = "log"),
                         na.action=na.omit)
      start.c.p <- exp(coef(start.coeff)[1])
      # x <- model.matrix(start.coeff)
      
      # start.coeff1 <- logbin(factor(outcome[i,index1]==0) ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1],
      #                        # family = binomial(link = "log"),
      #                        start = c(log(start.p[i]), rep(0,5)), maxit = 25000,
      #                        method = "em",
      #                        na.action=na.omit)
      # start.z.p <- 1 - exp(coef(start.coeff1)[1])
      
      # fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) + x1[index1] + outcome[i-1,index1],
      #                  na.action=na.omit, method = "Nelder-Mead", link = "log",
      #                  start = list(count = c(log(start.c.p), rep(0, 5)), 
      #                               zero = c(log(start.z.p), rep(0, 5))))
      
      fit2 <- zeroinfl(outcome[i,index1] ~ factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1] | factor(trt[index1]) + factor(mediator[i,index1]) + factor(tvcon[i,index1]) 
                       + x1[index1] + x2[index1] + x3[index1] + outcome[i-1,index1],
                       na.action=na.omit, method = "Nelder-Mead", link = "log",
                       start = list(count = c(log(start.c.p), rep(0, 6)), 
                                    zero = c(log(start.p[i]), rep(0, 7))))
      
      # o_count <- append(o_count, fit2$coefficients$count[[8]])
      o_zero <- append(o_zero, fit2$coefficients$zero[[8]])
    }
    b0_count <- append(b0_count, fit2$coefficients$count[[1]])
    g_count <- append(g_count, fit2$coefficients$count[[2]])
    b1_count <- append(b1_count, fit2$coefficients$count[[3]])
    lo_count <- append(lo_count, fit2$coefficients$count[[4]])
    x_1o_count <- append(x_1o_count, fit2$coefficients$count[[5]])
    x_2o_count <- append(x_2o_count, fit2$coefficients$count[[6]])
    x_3o_count <- append(x_3o_count, fit2$coefficients$count[[7]])
    
    
    b0_zero <- append(b0_zero, fit2$coefficients$zero[[1]])
    g_zero <- append(g_zero, fit2$coefficients$zero[[2]])
    b1_zero <- append(b1_zero, fit2$coefficients$zero[[3]])
    lo_zero <- append(lo_zero, fit2$coefficients$zero[[4]])
    x_1o_zero <- append(x_1o_zero, fit2$coefficients$zero[[5]])
    x_2o_zero <- append(x_2o_zero, fit2$coefficients$zero[[6]])
    x_3o_zero <- append(x_3o_zero, fit2$coefficients$zero[[7]])
  }
  
  # mat_b0[r,] <- b0_count
  mat_g[r,] <- predict(loess(exp(g_count) ~ t.seq[1:nm], span = 0.5, degree = 1))
  # mat_b1[r,] <- predict(loess(b1_count ~ t.seq[1:nm], span = 0.5, degree = 1))
  # mat_lo[r,] <- lo_count
  # mat_x1o[r,] <- x_1o_count
  # # mat_x2o[r,] <- x_2o_count
  # mat_o[r,] <- o_count
  # 
  # mat_b0_zero[r,] <- b0_zero
  # mat_g_zero[r,] <- g_zero
  # mat_b1_zero[r,] <- b1_zero
  # mat_lo_zero[r,] <- lo_zero
  # mat_x1o_zero[r,] <- x_1o_zero
  # # mat_x2o_zero[r,] <- x_2o_zero
  # mat_o_zero[r,] <- o_zero
  mat_ab[r,] <- predict(loess(((exp(a1)/(1+exp(a1)))*exp(b1_count)) ~ t.seq[1:nm], span = 0.5, degree = 1))
  
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
  
  for(i in 1:nm){
    if(i == 1){
      ### Overall Indirect Effect ###
      ## At baseline L(t+1)=1; C=cov_avg
      orNIE_l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
      ## At baseline L(t+1)=0; C=cov_avg
      orNIE_l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
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
      rrIES_l1 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
      ## At baseline L(t+1)=0; C=cov_avg
      rrIES_l0 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
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
      rrIENS_l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
      ## At baseline L(t+1)=0; C=cov_avg
      rrIENS_l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
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
      orNIE_m1l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      orNIE_m0l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      orNIE_m1l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      orNIE_m0l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
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
      rrIES_m1l1 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      rrIES_m0l1 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      rrIES_m1l0 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      rrIES_m0l0 <- ((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + b1_zero[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
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
      rrIENS_m1l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=1; C=cov_avg
      rrIENS_m0l1 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + lm[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
      ## At t+1 M(t)=1; L(t+1)=0; C=cov_avg
      rrIENS_m1l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1]))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2 + m[i-1])))
      
      ## At t+1 M(t)=0; L(t+1)=0; C=cov_avg
      rrIENS_m0l0 <- ((1 + exp(b1_count[i] + b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))/((1 + exp(b1_zero[i] + a0[i] + a1[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2))*(1 + exp(b1_count[i] + a0[i] + x_1m[i]*cov_x1 + x_2m[i]*cov_x2)))
      
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
  
  mat_IE11[r,] <- predict(loess(IE11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE01[r,] <- predict(loess(IE01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE10[r,] <- predict(loess(IE10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IE00[r,] <- predict(loess(IE00 ~ t.seq[1:nm], span = 0.5, degree = 1))
  
  mat_IES11[r,] <- predict(loess(IES11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IES01[r,] <- predict(loess(IES01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IES10[r,] <- predict(loess(IES10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IES00[r,] <- predict(loess(IES00 ~ t.seq[1:nm], span = 0.5, degree = 1))
  
  mat_IENS11[r,] <- predict(loess(IENS11 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IENS01[r,] <- predict(loess(IENS01 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IENS10[r,] <- predict(loess(IENS10 ~ t.seq[1:nm], span = 0.5, degree = 1))
  mat_IENS00[r,] <- predict(loess(IENS00 ~ t.seq[1:nm], span = 0.5, degree = 1))
}

effect_interest <- vector("list",14)

effect_interest[[1]] <- mat_g
effect_interest[[2]] <- mat_ab
effect_interest[[3]] <- mat_IE11
effect_interest[[4]] <- mat_IE01
effect_interest[[5]] <- mat_IE10
effect_interest[[6]] <- mat_IE00

effect_interest[[7]] <- mat_IES11
effect_interest[[8]] <- mat_IES01
effect_interest[[9]] <- mat_IES10
effect_interest[[10]] <- mat_IES00

effect_interest[[11]] <- mat_IENS11
effect_interest[[12]] <- mat_IENS01
effect_interest[[13]] <- mat_IENS10
effect_interest[[14]] <- mat_IENS00

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

# save(CI, file = paste("CI__varen_medf_count_female.rda"))
# save(CI, file = paste("CI__cNRT_medf_count_female.rda"))
# 
save(CI, file = paste("CI__varen_medp_count_female.rda"))
# save(CI, file = paste("CI__cNRT_medp_count_female.rda"))

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
# 
# 

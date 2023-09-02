### Bias analysis for binary outcome -  rare ###

setwd("~/Aim1/binary_r_var1")

library("RColorBrewer")
library(forcats)
library(reshape2)

# results <- get(load("~/Aim1/binary_r_var1/results_binary_200_modelr.rda"))
results <- get(load("~/Aim1/binary_r_var1/results_binary_500_modelr.rda"))
# results <- get(load("~/Aim1/binary_r_var1/results_binary_700_modelr.rda"))

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
o <- 0.9

beta0 <- -0.1 -2
gamma_A <- -0.1*exp(1/j)
beta1_r <- -(1/exp(1/j))-1.9

true_coeff <- vector("list",16)
true_coeff[[1]] <- alpha0
true_coeff[[2]] <- alpha1
true_coeff[[3]] <- lm
true_coeff[[4]] <- x1m
true_coeff[[5]] <- m
true_coeff[[6]] <- beta0
true_coeff[[7]] <- exp(gamma_A)
true_coeff[[8]] <- beta1_r
true_coeff[[9]] <- lo
true_coeff[[10]] <- x1o
true_coeff[[11]] <- o
true_coeff[[12]] <- (exp(alpha1)/(1+exp(alpha1)))*exp(beta1_r)

for(i in 1:20){
  if(i == 1){
    ## At baseline L(t+1)=1; C=2.25
    orNIE_l1 <- ((1 + exp(beta1_r[i] + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_r[i] + alpha0 + lm + x1m*2.25)))
    
    ## At baseline L(t+1)=0; C=2.25
    orNIE_l0 <- ((1 + exp(beta1_r[i] + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_r[i] + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[13]][i] <- orNIE_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[14]][i] <- orNIE_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[15]][i] <- orNIE_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[16]][i] <- orNIE_l0
    
  }else{
    ## At t+1 M(t)=1; L(t+1)=1; C=2.25
    orNIE_m1l1 <- ((1 + exp(beta1_r[i] + alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(alpha0 + lm + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(beta1_r[i] + alpha0 + lm + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=2.25
    orNIE_m0l1 <- ((1 + exp(beta1_r[i] + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_r[i] + alpha0 + lm + x1m*2.25)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=2.25
    orNIE_m1l0 <- ((1 + exp(beta1_r[i] + alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(alpha0 + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(beta1_r[i] + alpha0 + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=2.25
    orNIE_m0l0 <- ((1 + exp(beta1_r[i] + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_r[i] + alpha0 + x1m*2.25)))
    
    ## IE at M(t)=1; L(t+1)=1; C=2.25
    true_coeff[[13]][i] <- orNIE_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=2.25
    true_coeff[[14]][i] <- orNIE_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=2.25
    true_coeff[[15]][i] <- orNIE_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=2.25
    true_coeff[[16]][i] <- orNIE_m0l0
  }
}

R = 2000
t.seq <- j

cov_avg <- vector()

for(i in 1:nrow(results$mat_id)){
  id <- results$mat_id[i,]
  data <- aim1_dat[aim1_dat$id %in% id,] %>%
    arrange(id, time_pt)
  cov_avg <- append(cov_avg, mean(data$x1))
}

### Smooth the coefficient estimates where applicable###
for (i in 2:(length(results))) {
  mat1 <- results[[i]]
  mat2 <- matrix(nrow = R, ncol = ncol(mat1))
  
  ### Although the following lines of code are the same for both the mediator and the outcome model,
  ### since the required link transformation takes place for the final expression, 
  ### keep the code execution separate for the mediator and the outcome
  
  # if(i < 7){
  #   if(length(true_coeff[[i-1]]) > 1){
  #     for (j in 1:nrow(mat1)) {
  #       t1 <- mat1[j,]
  #       nm <- ncol(mat1)
  #       t2 <- loess(t1 ~ t.seq[1:nm], span = 0.5, degree = 1)
  #       mat2[j,] <- t2$fitted
  #     }
  #   }else if(length(true_coeff[[i-1]]) == 1){
  #     mat2 <- mat1
  #   }
  # }else{
  #   if(length(true_coeff[[i-1]]) > 1){
  #     for (j in 1:nrow(mat1)) {
  #       t1 <- mat1[j,]
  #       nm <- ncol(mat1)
  #       t2 <- loess(t1 ~ t.seq[1:nm], span = 0.5, degree = 1)
  #       mat2[j,] <- t2$fitted
  #     }
  #   }else if(length(true_coeff[[i-1]]) == 1){
  #     mat2 <- mat1
  #   }
  # }
  assign(paste0("est_", i-1), mat1)
}

est_12 <- matrix(nrow = R, ncol = ncol(est_2))
est_13 <- matrix(nrow = R, ncol = ncol(est_2))
est_14 <- matrix(nrow = R, ncol = ncol(est_2))
est_15 <- matrix(nrow = R, ncol = ncol(est_2))
est_16 <- matrix(nrow = R, ncol = ncol(est_2))

for (i in 1:ncol(est_2)) {
  
  est_12[,i] <- (exp(est_2[,i])/(1+exp(est_2[,i])))*exp(est_8[,i])
  
  if(i == 1){
    ## At baseline L(t+1)=1; C=avg cov
    orNIE_l1 <- ((1 + exp(est_8[,i] + est_1[,i] + est_2[,i] + est_3[,i] + est_4[,i]*cov_avg))*(1 + exp(est_1[,i] + est_3[,i] + est_4[,i]*cov_avg)))/((1 + exp(est_1[,i] + est_2[,i] + est_3[,i] + est_4[,i]*cov_avg))*(1 + exp(est_8[,i] + est_1[,i] + est_3[,i] + est_4[,i]*cov_avg)))
    
    ## At baseline L(t+1)=0; C=avg cov
    orNIE_l0 <- ((1 + exp(est_8[,i] + est_1[,i] + est_2[,i] + est_4[,i]*cov_avg))*(1 + exp(est_1[,i] + est_4[,i]*cov_avg)))/((1 + exp(est_1[,i] + est_2[,i] + est_4[,i]*cov_avg))*(1 + exp(est_8[,i] + est_1[,i] + est_4[,i]*cov_avg)))
    
    ## IE at M(t)=1; L(t+1)=1; C=avg cov
    est_13[,i] <- orNIE_l1
    
    ## IE at M(t)=0; L(t+1)=1; C=avg cov
    est_14[,i] <- orNIE_l1
    
    ## IE at M(t)=1; L(t+1)=0; C=avg cov
    est_15[,i] <- orNIE_l0
    
    ## IE at M(t)=0; L(t+1)=0; C=avg cov
    est_16[,i] <- orNIE_l0
    
  }else{
    ## At t+1 M(t)=1; L(t+1)=1; C=avg cov
    orNIE_m1l1 <- ((1 + exp(est_8[,i] + est_1[,i] + est_2[,i] + est_3[,i] + est_4[,i]*cov_avg + est_5[,i-1]))*(1 + exp(est_1[,i] + est_3[,i] + est_4[,i]*cov_avg + est_5[,i-1])))/((1 + exp(est_1[,i] + est_2[,i] + est_3[,i] + est_4[,i]*cov_avg + est_5[,i-1]))*(1 + exp(est_8[,i] + est_1[,i] + est_3[,i] + est_4[,i]*cov_avg + est_5[,i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=avg cov
    orNIE_m0l1 <- ((1 + exp(est_8[,i] + est_1[,i] + est_2[,i] + est_3[,i] + est_4[,i]*cov_avg))*(1 + exp(est_1[,i] + est_3[,i] + est_4[,i]*cov_avg)))/((1 + exp(est_1[,i] + est_2[,i] + est_3[,i] + est_4[,i]*cov_avg))*(1 + exp(est_8[,i] + est_1[,i] + est_3[,i] + est_4[,i]*cov_avg)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=avg cov
    orNIE_m1l0 <- ((1 + exp(est_8[,i] + est_1[,i] + est_2[,i] + est_4[,i]*cov_avg + est_5[,i-1]))*(1 + exp(est_1[,i] + est_4[,i]*cov_avg + est_5[,i-1])))/((1 + exp(est_1[,i] + est_2[,i] + est_4[,i]*cov_avg + est_5[,i-1]))*(1 + exp(est_8[,i] + est_1[,i] + est_4[,i]*cov_avg + est_5[,i-1])))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=avg cov
    orNIE_m0l0 <- ((1 + exp(est_8[,i] + est_1[,i] + est_2[,i] + est_4[,i]*cov_avg))*(1 + exp(est_1[,i] + est_4[,i]*cov_avg)))/((1 + exp(est_1[,i] + est_2[,i] + est_4[,i]*cov_avg))*(1 + exp(est_8[,i] + est_1[,i] + est_4[,i]*cov_avg)))
    
    ## IE at M(t)=1; L(t+1)=1; C=avg cov
    est_13[,i] <- orNIE_m1l1
    
    ## IE at M(t)=0; L(t+1)=1; C=avg cov
    est_14[,i] <- orNIE_m0l1
    
    ## IE at M(t)=1; L(t+1)=0; C=avg cov
    est_15[,i] <- orNIE_m1l0
    
    ## IE at M(t)=0; L(t+1)=0; C=avg cov
    est_16[,i] <- orNIE_m0l0
  }
}

effect_12 <- est_12
effect_13 <- est_13
effect_14 <- est_14
effect_15 <- est_15
effect_16 <- est_16


for(i in 1:nrow(effect_12)){
  est_12[i,] <- predict(loess(effect_12[i,] ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
  est_13[i,] <- predict(loess(effect_13[i,] ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
  est_14[i,] <- predict(loess(effect_14[i,] ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
  est_15[i,] <- predict(loess(effect_15[i,] ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
  est_16[i,] <- predict(loess(effect_16[i,] ~ t.seq[1:length(t.seq)], span = 0.5, degree = 1))
}

for (i in 2:(length(results))) {
  mat1 <- results[[i]]
  mat2 <- matrix(nrow = R, ncol = ncol(mat1))
  
  ### Although the following lines of code are the same for both the mediator and the outcome model,
  ### since the required link transformation takes place for the final expression, 
  ### keep the code execution separate for the mediator and the outcome
  
  if(i < 7){
    if(length(true_coeff[[i-1]]) > 1){
      for (j in 1:nrow(mat1)) {
        t1 <- mat1[j,]
        nm <- ncol(mat1)
        t2 <- loess(t1 ~ t.seq[1:nm], span = 0.5, degree = 1)
        mat2[j,] <- t2$fitted
      }
    }else if(length(true_coeff[[i-1]]) == 1){
      mat2 <- mat1
    }
  }else{
    if(length(true_coeff[[i-1]]) > 1){
      for (j in 1:nrow(mat1)) {
        t1 <- mat1[j,]
        nm <- ncol(mat1)
        t2 <- loess(t1 ~ t.seq[1:nm], span = 0.5, degree = 1)
        mat2[j,] <- t2$fitted
      }
    }else if(length(true_coeff[[i-1]]) == 1){
      mat2 <- mat1
    }
  }
  assign(paste0("est_", i-1), mat2)
}

est_7 <- exp(est_7)

est_smooth <- vector("list", 16)

for (i in 1:(length(true_coeff))) {
  est_smooth[[i]] <-  get(paste0("est_", i))
}

### Bias analysis ###
columns_tv_made = c("a1_m", "g_m", "b1_m", "ab_m",
                    "IEm1l1_m", "IEm0l1_m", "IEm1l0_m", "IEm0l0_m")
columns_tv_wase = c("a1_w", "g_w", "b1_w", "ab_w", 
                    "IEm1l1_w", "IEm0l1_w", "IEm1l0_w", "IEm0l0_w")
columns_tv_mase_nts = c("a1_nts", "g_nts", "b1_nts", "ab_nts",
                        "IEm1l1_nts", "IEm0l1_nts", "IEm1l0_nts", "IEm0l0_nts")
columns_tv_mase_nsts = c("a1_nsts", "g_nsts", "b1_nsts", "ab_nsts",
                         "IEm1l1_nsts", "IEm0l1_nsts", "IEm1l0_nsts", "IEm0l0_nsts")

columns = c("a0", "lm", "x1m", "m",
            "b0", "lo", "x1o", "o")

bias_tv_made <- as.data.frame(matrix(nrow = R, ncol = 8))
bias_tv_wase <- as.data.frame(matrix(nrow = R, ncol = 8))
bias_tv_mase_nts <- as.data.frame(matrix(nrow = R, ncol = 8))
bias_tv_mase_nsts <- as.data.frame(matrix(nrow = R, ncol = 8))

bias <- as.data.frame(matrix(nrow = R, ncol = 8))

colnames(bias_tv_made) <- columns_tv_made
colnames(bias_tv_wase) <- columns_tv_wase
colnames(bias_tv_mase_nts) <- columns_tv_mase_nts
colnames(bias_tv_mase_nsts) <- columns_tv_mase_nsts

colnames(bias) <- columns

tv=0
f=0

for (i in 1:length(true_coeff)){
  mat3 <- est_smooth[[i]]
  mat4 <- true_coeff[[i]]
  nm <- ncol(mat3)
  
  MADE <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  WASE <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  MASE_nsts <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  MASE_nts <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  RMSE <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  
  if(length(true_coeff[[i]]) > 1){
    tv <- tv+1
    
    denom_nsts <- vector()
    for (k in 2:length(t.seq)) {
      denom_nsts <- append(denom_nsts, abs(mat4[k] - mat4[k-1]))
    }
    
    denom_nts <- vector()
    for (k in 1:length(t.seq)) {
      denom_nts <- append(denom_nts, abs(mat4[k] - mean(mat4)))
    }
    
    for(j in 1:R){
      MADE[j,] <- abs(mat4 - mat3[j,])/(max(mat4) - min(mat4))
      WASE[j,] <- ((mat4 - mat3[j,])^2)/((max(mat4) - min(mat4))^2)
      
      # MADE[j,] <- abs((mat4 - mat3[j,])/mat4)
      
      MASE_nsts[j,] <- abs(mat4 - mat3[j,])/mean(denom_nsts)
      MASE_nts[j,] <- abs(mat4 - mat3[j,])/mean(denom_nts)
      
    }
    
    t_made <- as.data.frame(rowMeans(MADE)/nm)
    t_wase <- as.data.frame(rowMeans(WASE)/nm)
    t_mase_nsts <- as.data.frame(rowMeans(MASE_nsts))
    t_mase_nts <- as.data.frame(rowMeans(MASE_nts))
    
    bias_tv_made[,tv] <- t_made
    bias_tv_wase[,tv] <- t_wase
    bias_tv_mase_nsts[,tv] <- t_mase_nsts
    bias_tv_mase_nts[,tv] <- t_mase_nts
    
  }else if(length(true_coeff[[i]]) == 1){
    f <- f+1
    for(j in 1:nrow(mat3)){
      RMSE[j,] <- ((mat4 - mat3[j,]))^2
    }
    # if(nm < 20){
    #   t <- c(NA, sqrt(rowMeans(RMSE)))
    # }else if(nm == 20){
    #   t <- as.data.frame(sqrt(rowMeans(RMSE)))
    # }
    t <- as.data.frame(sqrt(rowMeans(RMSE)))
    bias[,f] <- t
  }
}

### Box plot for MADE, WASE and RMSE ###

bias_tv_made <- bias_tv_made %>%
  dplyr::select(-c(a1_m, b1_m)) 
  # bias_tv_made <- bias_tv_made[-1625,]  

test1 <- melt(bias_tv_made)
# test1$var <- substring(test1$variable, 1, 2)
test1$type <- "MADE"
# test1$var <- factor(test1$variable, 
#                     levels = c("a1_m", "g_m", "b1_m", "ab_m", "IEm1l1_m", "IEm0l1_m", "IEm1l0_m", "IEm0l0_m"),
#                     labels = c("α1(t)", "γ(t) - DE", "β1(t)", "α1(t)×β1(t)", 
#                                "POF IE11", "POF IE01", "POF IE10", "POF IE00"))

test1$var <- factor(test1$variable, 
                    levels = c("g_m", "ab_m", "IEm1l1_m", "IEm0l1_m", "IEm1l0_m", "IEm0l0_m"),
                    labels = c("DE", "POC IE", 
                               "POF IE11", "POF IE01", "POF IE10", "POF IE00"))

# bias_tv_wase <- bias_tv_wase[-1625,]  

bias_tv_wase <- bias_tv_wase %>%
  dplyr::select(-c(a1_w, b1_w))

test2 <- melt(bias_tv_wase)
# test2$var <- substring(test2$variable, 1, 2)
test2$type <- "WASE"
# test2$var <- factor(test2$variable, 
#                     levels = c("a1_w", "g_w", "b1_w", "ab_w", "IEm1l1_w", "IEm0l1_w", "IEm1l0_w", "IEm0l0_w"),
#                     labels = c("α1(t)", "γ(t) - DE", "β1(t)", "α1(t)×β1(t)", 
#                                "POF IE11", "POF IE01", "POF IE10", "POF IE00"))

test2$var <- factor(test2$variable, 
                    levels = c("g_w", "ab_w", "IEm1l1_w", "IEm0l1_w", "IEm1l0_w", "IEm0l0_w"),
                    labels = c("DE", "POC IE", 
                               "POF IE11", "POF IE01", "POF IE10", "POF IE00"))

tv_bias <- rbind(test1, test2)
tv_bias1 <- tv_bias %>%
  filter(variable %in% c("g_m", "g_w", "ab_m", "ab_w"))
tv_bias2 <- tv_bias %>%
  filter(variable %in% c("IEm1l1_m", "IEm0l1_m", "IEm1l0_m", "IEm0l0_m",
                         "IEm1l1_w", "IEm0l1_w", "IEm1l0_w", "IEm0l0_w"))

p1 <- ggplot(tv_bias1, aes(x=var, y=value, fill=type)) + 
  geom_boxplot() +
  facet_grid(type ~ var, scale="free") +
  labs(# title = "Bias estimates for the time-varying coefficients/effects of interest",
    x = "Coeffcient/Effect of Interest"
    # y = "Bias") +
  ) + 
  guides(fill=guide_legend(title="Type of Bias Measure")) +
  theme(axis.title.y = element_blank())

p2 <- ggplot(tv_bias2, aes(x=var, y=value, fill=type)) + 
  geom_boxplot() +
  facet_grid(type ~ var, scale="free") +
  labs(# title = "Bias estimates for the time-varying coefficients/effects of interest",
    x = "Coeffcient/Effect of Interest"
    # y = "Bias") +
  ) + 
  guides(fill=guide_legend(title="Type of Bias Measure")) +
  theme(axis.title.y = element_blank())

tiv_bias <- melt(bias)
tiv_bias$var <- factor(tiv_bias$variable, 
                       levels = c("a0", "lm", "x1m", "m", "b0", "lo", "x1o", "o"),
                       labels = c("α0", "α2", "α3'", "α4", "β0", "β2", "β3''", "β4"))

p3 <- ggplot(tiv_bias, aes(x=var, y=value)) + 
  geom_boxplot(color="black", fill="#69b3a2", alpha=0.2) +
  labs(title = "Bias estimates for the time-invariant coefficients of interest",
       x = "Coeffcient of Interest"
       # y = "Bias") +
  ) +
  theme(axis.title.y = element_blank())


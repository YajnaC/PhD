### Bias analysis for count outcome - not zero-inflated ###

library("RColorBrewer")
library(forcats)
results <- get(load("~/Aim1/results_binary_500_modelr.rda"))

### True coefficients of association ###

j <- seq(from = 1, to = 20, length.out = 20)

alpha0 <- -0.4
alpha1 <- exp((-0.25 -(2/j)))

lm <- -1.2
lo <- 0.36

x1m <- -0.05
x1o <- -0.02

x2m <- 0.3
x2o <- 0.06

m <- 1.25
o <- 0.05

beta0 <- -0.1 + 0.5
gamma_A <- -0.1*exp(1/j)
beta1_r <- -(1/exp(1/j))-1.9

true_coeff <- vector("list",13)
true_coeff[[1]] <- alpha0
true_coeff[[2]] <- alpha1
true_coeff[[3]] <- lm
true_coeff[[4]] <- x1m
true_coeff[[5]] <- x2m
true_coeff[[6]] <- m
true_coeff[[7]] <- beta0
true_coeff[[8]] <- gamma_A
true_coeff[[9]] <- beta1_r
true_coeff[[10]] <- lo
true_coeff[[11]] <- x1o
true_coeff[[12]] <- x2o
true_coeff[[13]] <- o

R = 2000
t.seq <- j

### Smooth the coefficient estimates where applicable###
for (i in 2:(length(results))) {
  mat1 <- results[[i]]
  mat2 <- matrix(nrow = R, ncol = ncol(mat1))
  
  ### Although the following lines of code are the same for both the mediator and the outcome model,
  ### since the required transformation takes place for the final expression, keep the code execution
  ### separate for the mediator and the outcome
  
  if(i < 8){
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

est_smooth <- vector("list", 13)

for (i in 1:(length(results)-1)) {
  est_smooth[[i]] <-  get(paste0("est_", i))
}

### Bias analysis ###
columns_tv_made = c("a1_m", "g_m", "b1_m")
columns_tv_wase = c("a1_w", "g_w", "b1_w")

columns = c("a0", "lm", "x1m", "x2m", "m",
            "b0", "lo", "x1o", "x2o", "o")

bias_tv_made <- as.data.frame(matrix(nrow = R, ncol = 3))
bias_tv_wase <- as.data.frame(matrix(nrow = R, ncol = 3))

bias <- as.data.frame(matrix(nrow = 20, ncol = 10))

colnames(bias_tv_made) <- columns_tv_made
colnames(bias_tv_wase) <- columns_tv_wase
colnames(bias) <- columns
tv=0
f=0
for (i in 1:length(true_coeff)){
  mat3 <- est_smooth[[i]]
  mat4 <- true_coeff[[i]]
  nm <- ncol(mat3)
  MADE <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  WASE <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  RMSE <- matrix(nrow = nrow(mat3), ncol = ncol(mat3))
  
  if(length(true_coeff[[i]]) > 1){
    tv <- tv+1
    for(j in 1:R){
      MADE[j,] <- abs(mat4 - mat3[j,])/(max(mat4) - min(mat4))
      WASE[j,] <- ((mat4 - mat3[j,])^2)/((max(mat4) - min(mat4))^2)
    }
    
    t_made <- as.data.frame(rowMeans(MADE)/nm)
    t_wase <- as.data.frame(rowMeans(WASE)/nm)
    
    bias_tv_made[,tv] <- t_made
    bias_tv_wase[,tv] <- t_wase
  }else if(length(true_coeff[[i]]) == 1){
    f <- f+1
    for(j in 1:nrow(mat3)){
      RMSE[j,] <- ((mat4 - mat3[j,]))^2
    }
    if(nm < 20){
      t <- c(NA, sqrt(colMeans(RMSE)/R))
    }
    else if(nm == 20){
      t <- as.data.frame(sqrt(colMeans(RMSE)/R))
    }
    bias[,f] <- t
  }
}

### Box plot for MADE/WASE ###
test1 <- melt(bias_tv_made)
test1$var <- substring(test1$variable, 1, 2)
test1$type <- "MADE"

test2 <- melt(bias_tv_wase)
test2$var <- substring(test2$variable, 1, 2)
test2$type <- "WASE"

tv_bias <- rbind(test1, test2)

p1 <- ggplot(tv_bias, aes(x=var, y=value, fill=type)) + 
  geom_boxplot() +
  facet_wrap(~var, scale="free")


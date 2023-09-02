### Combine all CI results and find coverage ###

setwd("~/Aim1/count_nzi")

CI_full <- vector("list", 16)
CI_full[[1]] <- vector("list",2000)
CI_full[[2]] <- vector("list",2000)
CI_full[[3]] <- vector("list",2000)
CI_full[[4]] <- vector("list",2000)
CI_full[[5]] <- vector("list",2000)
CI_full[[6]] <- vector("list",2000)
CI_full[[7]] <- vector("list",2000)
CI_full[[8]] <- vector("list",2000)
CI_full[[9]] <- vector("list",2000)
CI_full[[10]] <- vector("list",2000)
CI_full[[11]] <- vector("list",2000)
CI_full[[12]] <- vector("list",2000)
CI_full[[13]] <- vector("list",2000)
CI_full[[14]] <- vector("list",2000)
CI_full[[15]] <- vector("list",2000)
CI_full[[16]] <- vector("list",2000)

for (k in 1:2000) {
  CI_part <- get(load(paste("~/Aim1/count_nzi/CI_count_nzi_500_",k,".rda", sep = "")))
  for(z in 1:16){
    CI_full[[z]][[k]] <- CI_part[[z]][[k]] 
  }
}

save(CI_full, file = paste("CI_count_nzi_500.rda", sep = ""))

CI <- get(load("~/Aim1/continuous/CI_count_nzi_500.rda"))

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

CI_cov <- vector("list", 16)

for(z in 1:length(true_coeff)){
  for (k in 1:length(CI[[z]])) {
    mat_interest <- CI[[z]][[k]]
    
    coverage <- ifelse(mat_interest[,1] < true_coeff[[z]] & 
                         mat_interest[,2] > true_coeff[[z]], 1, 0)
    # if(length(true_coeff[[z]]) > 1){
    #   coverage <- ifelse(mat_interest[,1] < true_coeff[[z]] & 
    #                        mat_interest[,2] > true_coeff[[z]], 1, 0)  
    # 
    # }else{
    #   coverage <- ifelse(mat_interest[,1] < true_coeff[[z]] & 
    #                        mat_interest[,2] > true_coeff[[z]], 1, 0)  
    # }
    CI_cov[[z]][[k]] <- coverage
  }
}

save(CI_cov, file = paste("CIcoverage_cont_500.rda", sep = ""))



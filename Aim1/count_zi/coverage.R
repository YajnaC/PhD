setwd("~/Aim1/count_zi")

load("~/Aim1/aim1_dat.rda")

CI_full <- vector("list", 30)
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
CI_full[[17]] <- vector("list",2000)
CI_full[[18]] <- vector("list",2000)
CI_full[[19]] <- vector("list",2000)
CI_full[[20]] <- vector("list",2000)
CI_full[[21]] <- vector("list",2000)
CI_full[[22]] <- vector("list",2000)
CI_full[[23]] <- vector("list",2000)
CI_full[[24]] <- vector("list",2000)
CI_full[[25]] <- vector("list",2000)
CI_full[[26]] <- vector("list",2000)
CI_full[[27]] <- vector("list",2000)
CI_full[[28]] <- vector("list",2000)
CI_full[[29]] <- vector("list",2000)
CI_full[[30]] <- vector("list",2000)

for (k in 1:2000) {
  CI_part <- get(load(paste("~/Aim1/count_zi/CI_count_zi_500_",k,".rda", sep = "")))
  for(z in 1:30){
    CI_full[[z]][[k]] <- CI_part[[z]][[k]] 
  }
}

save(CI_full, file = paste("CI_count_zi_500_full.rda", sep = ""))

CI <- get(load("~/Aim1/count_zi/CI_count_zi_500_full.rda"))

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

CI_cov <- vector("list", 30)

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

save(CI_cov, file = paste("CIcoverage_count_zi_500.rda", sep = ""))

cov_perc <- matrix(nrow = length(j), ncol = length(CI_cov))

for (i in 1:length(CI_cov)) {
  temp <- as.data.frame(CI_cov[[i]])
  temp$cov <- (rowSums(temp)/2000)*100
  
  if(nrow(temp) == 20){
    cov_perc[,i] <- temp$cov
  }else{
    cov_perc[,i] <- append(NA, temp$cov)
  }
}

var_names <- c("α0", "α1(t)", "α2", "α3'", "α4", 
               "β0", "γ(t) - DE", "β1(t)", "β2", "β3'", "β4",
               "α1(t)×β1(t)", "POF OIE11", "POF OIE01", "POF OIE10", "POF OIE00",
               "POF IES11", "POF IES01", "POF IES10", "POF IES00",
               "POF IENS11", "POF IENS01", "POF IENS10", "POF IENS00",
               "β0'", "γ(t)'", "β1(t)'", "β2'", "β3''", "β4'")

colnames(cov_perc) <- var_names
cov_perc_select <- cov_perc[,c(7,12:24)]

write.csv(cov_perc_select, paste("coverage_countzi.csv"), row.names=FALSE)

# melted_cov_perc <- melt(cov_perc_select) 
# 
# ggplot(melted_cov_perc, aes(fill=Var2, y=value, x=Var1)) + 
#   geom_bar(position="dodge", stat="identity")

med_ab_made <- median(bias_tv_made$ab_m)
med_diff <- vector()
for(i in 1:nrow(bias_tv_made)){
  dif_med <- abs(bias_tv_made$ab_m[i]-med_ab_made)
  med_diff <- append(med_diff,dif_med)
}

med_diff <- as.data.frame(med_diff)
id_med = 903
### Sample ID 903 rendered a product of coefficients estimate closest to the median ###

t.seq <- j

effect_interest <- vector("list",14)
effect_interest[[1]] <- est_smooth[[7]]
effect_interest[[2]] <- est_smooth[[18]]

effect_interest[[3]] <- est_smooth[[19]]
effect_interest[[4]] <- est_smooth[[20]]
effect_interest[[5]] <- est_smooth[[21]]
effect_interest[[6]] <- est_smooth[[22]]

effect_interest[[7]] <- est_smooth[[23]]
effect_interest[[8]] <- est_smooth[[24]]
effect_interest[[9]] <- est_smooth[[25]]
effect_interest[[10]] <- est_smooth[[26]]

effect_interest[[11]] <- est_smooth[[27]]
effect_interest[[12]] <- est_smooth[[28]]
effect_interest[[13]] <- est_smooth[[29]]
effect_interest[[14]] <- est_smooth[[30]]

t.est <- seq(from = 1, to = 20, length.out = 20)

dat_g <- as.data.frame(cbind(effect_interest[[1]][id_med,], CI[[7]][[id_med]], t.est))
dat_ab <- as.data.frame(cbind(effect_interest[[2]][id_med,], CI[[12]][[id_med]], t.est))

dat_OIE11 <- as.data.frame(cbind(effect_interest[[3]][id_med,], CI[[13]][[id_med]]))
dat_OIE01 <- as.data.frame(cbind(effect_interest[[4]][id_med,], CI[[14]][[id_med]]))
dat_OIE10 <- as.data.frame(cbind(effect_interest[[5]][id_med,], CI[[15]][[id_med]]))
dat_OIE00 <- as.data.frame(cbind(effect_interest[[6]][id_med,], CI[[16]][[id_med]]))

dat_IES11 <- as.data.frame(cbind(effect_interest[[7]][id_med,], CI[[17]][[id_med]]))
dat_IES01 <- as.data.frame(cbind(effect_interest[[8]][id_med,], CI[[18]][[id_med]]))
dat_IES10 <- as.data.frame(cbind(effect_interest[[9]][id_med,], CI[[19]][[id_med]]))
dat_IES00 <- as.data.frame(cbind(effect_interest[[10]][id_med,], CI[[20]][[id_med]]))

dat_IENS11 <- as.data.frame(cbind(effect_interest[[11]][id_med,], CI[[21]][[id_med]]))
dat_IENS01 <- as.data.frame(cbind(effect_interest[[12]][id_med,], CI[[22]][[id_med]]))
dat_IENS10 <- as.data.frame(cbind(effect_interest[[13]][id_med,], CI[[23]][[id_med]]))
dat_IENS00 <- as.data.frame(cbind(effect_interest[[14]][id_med,], CI[[24]][[id_med]]))


dat_POF <- as.data.frame(cbind(dat_OIE11, dat_OIE01, dat_OIE10, dat_OIE00,
                               dat_IES11, dat_IES01, dat_IES10, dat_IES00,
                               dat_IENS11, dat_IENS01, dat_IENS10, dat_IENS00, t.est))

var_names <- c("OIE11", "LL_OIE11", "UL_OIE11",
               "OIE01", "LL_OIE01", "UL_OIE01",
               "OIE10", "LL_OIE10", "UL_OIE10",
               "OIE00", "LL_OIE00", "UL_OIE00",
               "IES11", "LL_IES11", "UL_IES11",
               "IES01", "LL_IES01", "UL_IES01",
               "IES10", "LL_IES10", "UL_IES10",
               "IES00", "LL_IES00", "UL_IES00",
               "IENS11", "LL_IENS11", "UL_IENS11",
               "IENS01", "LL_IENS01", "UL_IENS01",
               "IENS10", "LL_IENS10", "UL_IENS10",
               "IENS00", "LL_IENS00", "UL_IENS00",
               "t.est")

colnames(dat_POF) <- var_names

ggplot(data = dat_g, aes(t.est, V1)) +
  geom_line(color = "blueviolet", linewidth = 0.8) +
  geom_ribbon(data=dat_g, aes(ymin=V2, ymax=V3), fill="lightcyan3", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[7]]), linewidth = 0.8, color = "darkgreen", linetype = "dashed") + 
  geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  labs(title = "Plot of the time-varying direct effect - γ(t)",
       x = "Time Sequence",
       y = "γ(t) Curve") +
  scale_x_continuous(breaks = t.est)

ggplot(data = dat_ab, aes(t.est, V1)) +
  geom_line(color = "darkgreen", linewidth = 0.75) +
  geom_ribbon(data=dat_ab, aes(ymin=V2, ymax=V3), fill="azure4", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[12]]), linewidth = 0.8, color = "firebrick", linetype = "dashed") + 
  # geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  labs(title = "Plot of the time-varying mediation effect - POC",
       x = "Time Sequence",
       y = "IE POC Curve") +
  scale_x_continuous(breaks = t.est)

plot_POF11 <- ggplot(data = dat_POF, aes(t.est, OIE11)) +
  geom_line(color = "red", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_OIE11, ymax=UL_OIE11), fill="lightpink", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[13]]), linewidth = 0.8, color = "darkred", linetype = "dashed") + 
  
  geom_line(aes(t.est, IES11), color = "blue", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IES11, ymax=UL_IES11), fill="lightblue", alpha=0.5) + 
  geom_line(aes(t.est, true_coeff[[17]]), linewidth = 0.8, color = "darkblue", linetype = "dashed") + 
  
  geom_line(aes(t.est, IENS11), color = "darkgreen", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IENS11, ymax=UL_IENS11), fill="yellow", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[21]]), linewidth = 0.8, color = "gold4", linetype = "dashed") +
  
  geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  
  labs(title = "Plot of the time-varying mediation effect - POF",
       x = "Time Sequence",
       y = "Indirect Effect Curves") +
  scale_x_continuous(breaks = t.est)

plot_POF01 <- ggplot(data = dat_POF, aes(t.est, OIE01)) +
  geom_line(color = "red", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_OIE01, ymax=UL_OIE01), fill="lightpink", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[13]]), linewidth = 0.8, color = "darkred", linetype = "dashed") + 
  
  geom_line(aes(t.est, IES01), color = "blue", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IES01, ymax=UL_IES01), fill="lightblue", alpha=0.5) + 
  geom_line(aes(t.est, true_coeff[[17]]), linewidth = 0.8, color = "darkblue", linetype = "dashed") + 
  
  geom_line(aes(t.est, IENS01), color = "darkgreen", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IENS01, ymax=UL_IENS01), fill="yellow", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[21]]), linewidth = 0.8, color = "gold4", linetype = "dashed") +
  
  geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  
  labs(title = "Plot of the time-varying mediation effect - POF",
       x = "Time Sequence",
       y = "Indirect Effect Curves") +
  scale_x_continuous(breaks = t.est)

plot_POF10 <- ggplot(data = dat_POF, aes(t.est, OIE10)) +
  geom_line(color = "red", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_OIE10, ymax=UL_OIE10), fill="lightpink", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[13]]), linewidth = 0.8, color = "darkred", linetype = "dashed") + 
  
  geom_line(aes(t.est, IES10), color = "blue", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IES10, ymax=UL_IES10), fill="lightblue", alpha=0.5) + 
  geom_line(aes(t.est, true_coeff[[17]]), linewidth = 0.8, color = "darkblue", linetype = "dashed") + 
  
  geom_line(aes(t.est, IENS10), color = "darkgreen", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IENS10, ymax=UL_IENS10), fill="yellow", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[21]]), linewidth = 0.8, color = "gold4", linetype = "dashed") +
  
  geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  
  labs(title = "Plot of the time-varying mediation effect - POF",
       x = "Time Sequence",
       y = "Indirect Effect Curves") +
  scale_x_continuous(breaks = t.est)


plot_POF00 <- ggplot(data = dat_POF, aes(t.est, OIE00)) +
  geom_line(color = "red", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_OIE00, ymax=UL_OIE00), fill="lightpink", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[13]]), linewidth = 0.8, color = "darkred", linetype = "dashed") + 
  
  geom_line(aes(t.est, IES00), color = "blue", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IES00, ymax=UL_IES00), fill="lightblue", alpha=0.5) + 
  geom_line(aes(t.est, true_coeff[[17]]), linewidth = 0.8, color = "darkblue", linetype = "dashed") + 
  
  geom_line(aes(t.est, IENS00), color = "darkgreen", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IENS00, ymax=UL_IENS00), fill="yellow", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[21]]), linewidth = 0.8, color = "gold4", linetype = "dashed") +
  
  geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  
  labs(title = "Plot of the time-varying mediation effect - POF",
       x = "Time Sequence",
       y = "Indirect Effect Curves") +
  scale_x_continuous(breaks = t.est)

# for (tp in 1:20) {
#   ### Direct Effect ###
#   datg <- effect_interest[[1]][,tp]
#   
#   tp1_datg <- matrix(nrow = 2000, ncol = 3)
#   
#   for (i in 1:length(datg)) {
#     tp1_datg[i,] <- c(datg[i],CI[[7]][[i]][tp,])
#   }
#   
#   g_data <- as.data.frame(tp1_datg)
#   g_data$SampleID <- rownames(g_data)
#   
#   g_data1 <- g_data[1:500,]
#   g_data2 <- g_data[501:1000,]
#   g_data3 <- g_data[1001:1500,]
#   g_data4 <- g_data[1501:2000,]
#   
#   ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
#     geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_DE1_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
#     geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_DE2_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
#     geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_DE3_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
#     geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_DE4_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   
#   ### POC-IE ###
#   datg <- effect_interest[[2]][,tp]
#   
#   tp1_datg <- matrix(nrow = 2000, ncol = 3)
#   
#   for (i in 1:length(datg)) {
#     tp1_datg[i,] <- c(datg[i],CI[[12]][[i]][tp,])
#   }
#   
#   g_data <- as.data.frame(tp1_datg)
#   g_data$SampleID <- rownames(g_data)
#   
#   g_data1 <- g_data[1:500,]
#   g_data2 <- g_data[501:1000,]
#   g_data3 <- g_data[1001:1500,]
#   g_data4 <- g_data[1501:2000,]
#   
#   ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
#     geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_POC1_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
#     geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_POC2_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
#     geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_POC3_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
#     geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_POC4_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   
#   ### POF-IE11 ###
#   datg <- effect_interest[[3]][,tp]
#   
#   tp1_datg <- matrix(nrow = 2000, ncol = 3)
#   
#   for (i in 1:length(datg)) {
#     tp1_datg[i,] <- c(datg[i],CI[[13]][[i]][tp,])
#   }
#   
#   g_data <- as.data.frame(tp1_datg)
#   g_data$SampleID <- rownames(g_data)
#   
#   g_data1 <- g_data[1:500,]
#   g_data2 <- g_data[501:1000,]
#   g_data3 <- g_data[1001:1500,]
#   g_data4 <- g_data[1501:2000,]
#   
#   ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
#     geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE11_1_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
#     geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE11_2_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
#     geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE11_3_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
#     geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE11_4_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ### POF-IE01 ###
#   datg <- effect_interest[[4]][,tp]
#   
#   tp1_datg <- matrix(nrow = 2000, ncol = 3)
#   
#   for (i in 1:length(datg)) {
#     tp1_datg[i,] <- c(datg[i],CI[[14]][[i]][tp,])
#   }
#   
#   g_data <- as.data.frame(tp1_datg)
#   g_data$SampleID <- rownames(g_data)
#   
#   g_data1 <- g_data[1:500,]
#   g_data2 <- g_data[501:1000,]
#   g_data3 <- g_data[1001:1500,]
#   g_data4 <- g_data[1501:2000,]
#   
#   ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
#     geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE01_1_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
#     geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE01_2_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
#     geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE01_3_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
#     geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE01_4_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ### POF-10 ###
#   datg <- effect_interest[[5]][,tp]
#   
#   tp1_datg <- matrix(nrow = 2000, ncol = 3)
#   
#   for (i in 1:length(datg)) {
#     tp1_datg[i,] <- c(datg[i],CI[[15]][[i]][tp,])
#   }
#   
#   g_data <- as.data.frame(tp1_datg)
#   g_data$SampleID <- rownames(g_data)
#   
#   g_data1 <- g_data[1:500,]
#   g_data2 <- g_data[501:1000,]
#   g_data3 <- g_data[1001:1500,]
#   g_data4 <- g_data[1501:2000,]
#   
#   ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
#     geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE10_1_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
#     geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE10_2_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
#     geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE10_3_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
#     geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE10_4_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   
#   ### POF-00 ###
#   datg <- effect_interest[[6]][,tp]
#   
#   tp1_datg <- matrix(nrow = 2000, ncol = 3)
#   
#   for (i in 1:length(datg)) {
#     tp1_datg[i,] <- c(datg[i],CI[[16]][[i]][tp,])
#   }
#   
#   g_data <- as.data.frame(tp1_datg)
#   g_data$SampleID <- rownames(g_data)
#   
#   g_data1 <- g_data[1:500,]
#   g_data2 <- g_data[501:1000,]
#   g_data3 <- g_data[1001:1500,]
#   g_data4 <- g_data[1501:2000,]
#   
#   ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
#     geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE00_1_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
#     geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE00_2_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
#     geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE00_3_",tp,".png"), width = 20, height = 12, units = "in")
#   
#   
#   ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
#     geom_point(color='blue') +
#     geom_errorbarh(height=.2) +
#     labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
#     geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
#     theme_minimal() + coord_flip() +
#     theme(axis.text.x = element_blank(),  #remove y axis labels
#           axis.ticks.x = element_blank())
#   ggsave(paste("CI_cov_IE00_4_",tp,".png"), width = 20, height = 12, units = "in")
#   
# }

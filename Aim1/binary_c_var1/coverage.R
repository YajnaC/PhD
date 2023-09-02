setwd("~/Aim1/binary_c_var1")

load("~/Aim1/aim1_dat.rda")

CI_full <- vector("list", 5)
CI_full[[1]] <- vector("list",2000)
CI_full[[2]] <- vector("list",2000)
CI_full[[3]] <- vector("list",2000)
CI_full[[4]] <- vector("list",2000)
CI_full[[5]] <- vector("list",2000)
# CI_full[[6]] <- vector("list",2000)
# CI_full[[7]] <- vector("list",2000)
# CI_full[[8]] <- vector("list",2000)
# CI_full[[9]] <- vector("list",2000)
# CI_full[[10]] <- vector("list",2000)
# CI_full[[11]] <- vector("list",2000)
# CI_full[[12]] <- vector("list",2000)
# CI_full[[13]] <- vector("list",2000)
# CI_full[[14]] <- vector("list",2000)
# CI_full[[15]] <- vector("list",2000)
# CI_full[[16]] <- vector("list",2000)

for (k in 1:2000) {
  CI_part <- get(load(paste("~/Aim1/binary_c_var1/CI_binaryc_500_",k,".rda", sep = "")))
  for(z in 1:5){
    CI_full[[z]][[k]] <- CI_part[[z]][[k]] 
  }
}

save(CI_full, file = paste("CI_binary_500_common.rda", sep = ""))

CI <- get(load("~/Aim1/binary_c_var1/CI_binary_500_common.rda"))

### True coefficients of association ###

j <- seq(from = 1, to = 20, length.out = 20)

alpha0 <- -0.4
alpha1 <- exp((-0.25 -(2/j)))

lm <- -1.2
lo <- 0.36

x1m <- 0.05
x1o <- -0.02

m <- 1.25
o <- 0.25

beta0 <- -0.1 - 0.74
gamma_A <- -0.1*exp(1/j)
beta1_c <- -(1/exp(1/(j^3)))

true_coeff <- vector("list",16)
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
    ## At baseline L(t+1)=1; C=2.25
    orNIE_l1 <- ((1 + exp(beta1_c[i] + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25)))
    
    ## At baseline L(t+1)=0; C=2.25
    orNIE_l0 <- ((1 + exp(beta1_c[i] + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25)))
    
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
    orNIE_m1l1 <- ((1 + exp(beta1_c[i] + alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(alpha0 + lm + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25 + m))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=1; C=2.25
    orNIE_m0l1 <- ((1 + exp(beta1_c[i] + alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(alpha0 + lm + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + lm + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + lm + x1m*2.25)))
    
    ## At t+1 M(t)=1; L(t+1)=0; C=2.25
    orNIE_m1l0 <- ((1 + exp(beta1_c[i] + alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(alpha0 + x1m*2.25 + m)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25 + m))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25 + m)))
    
    ## At t+1 M(t)=0; L(t+1)=0; C=2.25
    orNIE_m0l0 <- ((1 + exp(beta1_c[i] + alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(alpha0 + x1m*2.25)))/((1 + exp(alpha0 + alpha1[i] + x1m*2.25))*(1 + exp(beta1_c[i] + alpha0 + x1m*2.25)))
    
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

CI_cov <- vector("list", 5)

for(z in 1:length(CI)){
  for (k in 1:length(CI[[z]])) {
    mat_interest <- CI[[z]][[k]]
    
    coverage <- ifelse(mat_interest[,1] < true_coeff[[z+11]] & 
                         mat_interest[,2] > true_coeff[[z+11]], 1, 0)
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

save(CI_cov, file = paste("CIcoverage_binc_500.rda", sep = ""))

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

var_names <- c("α1(t)×β1(t)", "POF IE11", "POF IE01", "POF IE10", "POF IE00")

colnames(cov_perc) <- var_names
cov_perc_select <- cov_perc

write.csv(cov_perc_select, paste("coverage_binr.csv"), row.names=FALSE)

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
id_med = 384
### Sample ID 511 rendered a product of coefficients estimate closest to the median ###

t.seq <- j

effect_interest <- vector("list",6)
effect_interest[[1]] <- est_smooth[[7]]
effect_interest[[2]] <- est_smooth[[12]]

effect_interest[[3]] <- est_smooth[[13]]
effect_interest[[4]] <- est_smooth[[14]]
effect_interest[[5]] <- est_smooth[[15]]
effect_interest[[6]] <- est_smooth[[16]]

t.est <- seq(from = 1, to = 20, length.out = 20)

# dat_g <- as.data.frame(cbind(effect_interest[[1]][id_med,], CI[[7]][[id_med]], t.est))
dat_ab <- as.data.frame(cbind(effect_interest[[2]][id_med,], CI[[1]][[id_med]], t.est))
dat_IE11 <- as.data.frame(cbind(effect_interest[[3]][id_med,], CI[[2]][[id_med]]))
dat_IE01 <- as.data.frame(cbind(effect_interest[[4]][id_med,], CI[[3]][[id_med]]))
dat_IE10 <- as.data.frame(cbind(effect_interest[[5]][id_med,], CI[[4]][[id_med]]))
dat_IE00 <- as.data.frame(cbind(effect_interest[[6]][id_med,], CI[[5]][[id_med]]))
dat_POF <- as.data.frame(cbind(dat_IE11, dat_IE01, dat_IE10, dat_IE00, t.est))

var_names <- c("IE11", "LL_IE11", "UL_IE11",
               "IE01", "LL_IE01", "UL_IE01",
               "IE10", "LL_IE10", "UL_IE10",
               "IE00", "LL_IE00", "UL_IE00", "t.est")

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

plot_POF1 <- ggplot(data = dat_POF, aes(t.est, IE11)) +
  geom_line(color = "red", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IE11, ymax=UL_IE11), fill="lightpink", alpha=0.5) +
  geom_line(aes(t.est, true_coeff[[13]]), linewidth = 0.8, color = "darkred", linetype = "dashed") + 
  
  geom_line(aes(t.est, IE01), color = "blue", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IE01, ymax=UL_IE01), fill="lightblue", alpha=0.5) + 
  geom_line(aes(t.est, true_coeff[[14]]), linewidth = 0.8, color = "darkblue", linetype = "dashed") + 
  
  geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  
  labs(title = "Plot of the time-varying mediation effect - POF",
       x = "Time Sequence",
       y = "Indirect Effect Curves") +
  scale_x_continuous(breaks = t.est)

plot_POF2 <- ggplot(data = dat_POF, aes(t.est, IE10)) +
  geom_line(color = "darkred", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IE10, ymax=UL_IE10), fill="yellow", alpha=0.5) + 
  geom_line(aes(t.est, true_coeff[[15]]), linewidth = 0.8, color = "red", linetype = "dashed") + 
  
  
  geom_line(aes(t.est, IE00), color = "darkgreen", linewidth = 0.75) +
  geom_ribbon(data=dat_POF, aes(ymin=LL_IE00, ymax=UL_IE00), fill="lightgreen", alpha=0.5) + 
  geom_line(aes(t.est, true_coeff[[16]]), linewidth = 0.8, color = "navy", linetype = "dashed") +
  
  geom_line(aes(t.est, 1), linewidth = 0.7, color = "black") +
  
  labs(title = "Plot of the time-varying mediation effect - POF",
       x = "Time Sequence",
       y = "Indirect Effect Curves") +
  scale_x_continuous(breaks = t.est)

for (tp in 1:20) {
  # ### Direct Effect ###
  # datg <- effect_interest[[1]][,tp]
  # 
  # tp1_datg <- matrix(nrow = 2000, ncol = 3)
  # 
  # for (i in 1:length(datg)) {
  #   tp1_datg[i,] <- c(datg[i],CI[[7]][[i]][tp,])
  # }
  # 
  # g_data <- as.data.frame(tp1_datg)
  # g_data$SampleID <- rownames(g_data)
  # 
  # g_data1 <- g_data[1:500,]
  # g_data2 <- g_data[501:1000,]
  # g_data3 <- g_data[1001:1500,]
  # g_data4 <- g_data[1501:2000,]
  # 
  # ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
  #   geom_point(color='blue') +
  #   geom_errorbarh(height=.2) +
  #   labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
  #   geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
  #   theme_minimal() + coord_flip() +
  #   theme(axis.text.x = element_blank(),  #remove y axis labels
  #         axis.ticks.x = element_blank())
  # ggsave(paste("CI_cov_DE1_",tp,".png"), width = 20, height = 12, units = "in")
  # 
  # ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
  #   geom_point(color='blue') +
  #   geom_errorbarh(height=.2) +
  #   labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
  #   geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
  #   theme_minimal() + coord_flip() +
  #   theme(axis.text.x = element_blank(),  #remove y axis labels
  #         axis.ticks.x = element_blank())
  # ggsave(paste("CI_cov_DE2_",tp,".png"), width = 20, height = 12, units = "in")
  # 
  # 
  # ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
  #   geom_point(color='blue') +
  #   geom_errorbarh(height=.2) +
  #   labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
  #   geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
  #   theme_minimal() + coord_flip() +
  #   theme(axis.text.x = element_blank(),  #remove y axis labels
  #         axis.ticks.x = element_blank())
  # ggsave(paste("CI_cov_DE3_",tp,".png"), width = 20, height = 12, units = "in")
  # 
  # 
  # ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
  #   geom_point(color='blue') +
  #   geom_errorbarh(height=.2) +
  #   labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
  #   geom_vline(xintercept=true_coeff[[7]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
  #   theme_minimal() + coord_flip() +
  #   theme(axis.text.x = element_blank(),  #remove y axis labels
  #         axis.ticks.x = element_blank())
  # ggsave(paste("CI_cov_DE4_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  
  ### POC-IE ###
  datg <- effect_interest[[2]][,tp]
  
  tp1_datg <- matrix(nrow = 2000, ncol = 3)
  
  for (i in 1:length(datg)) {
    tp1_datg[i,] <- c(datg[i],CI[[1]][[i]][tp,])
  }
  
  g_data <- as.data.frame(tp1_datg)
  g_data$SampleID <- rownames(g_data)
  
  g_data1 <- g_data[1:500,]
  g_data2 <- g_data[501:1000,]
  g_data3 <- g_data[1001:1500,]
  g_data4 <- g_data[1501:2000,]
  
  ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
    geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_POC1_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
    geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_POC2_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
    geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_POC3_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
    geom_vline(xintercept=true_coeff[[12]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_POC4_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  
  ### POF-IE11 ###
  datg <- effect_interest[[3]][,tp]
  
  tp1_datg <- matrix(nrow = 2000, ncol = 3)
  
  for (i in 1:length(datg)) {
    tp1_datg[i,] <- c(datg[i],CI[[2]][[i]][tp,])
  }
  
  g_data <- as.data.frame(tp1_datg)
  g_data$SampleID <- rownames(g_data)
  
  g_data1 <- g_data[1:500,]
  g_data2 <- g_data[501:1000,]
  g_data3 <- g_data[1001:1500,]
  g_data4 <- g_data[1501:2000,]
  
  ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
    geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE11_1_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
    geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE11_2_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
    geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE11_3_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
    geom_vline(xintercept=true_coeff[[13]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE11_4_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ### POF-IE01 ###
  datg <- effect_interest[[4]][,tp]
  
  tp1_datg <- matrix(nrow = 2000, ncol = 3)
  
  for (i in 1:length(datg)) {
    tp1_datg[i,] <- c(datg[i],CI[[3]][[i]][tp,])
  }
  
  g_data <- as.data.frame(tp1_datg)
  g_data$SampleID <- rownames(g_data)
  
  g_data1 <- g_data[1:500,]
  g_data2 <- g_data[501:1000,]
  g_data3 <- g_data[1001:1500,]
  g_data4 <- g_data[1501:2000,]
  
  ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
    geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE01_1_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
    geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE01_2_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
    geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE01_3_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
    geom_vline(xintercept=true_coeff[[14]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE01_4_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ### POF-10 ###
  datg <- effect_interest[[5]][,tp]
  
  tp1_datg <- matrix(nrow = 2000, ncol = 3)
  
  for (i in 1:length(datg)) {
    tp1_datg[i,] <- c(datg[i],CI[[4]][[i]][tp,])
  }
  
  g_data <- as.data.frame(tp1_datg)
  g_data$SampleID <- rownames(g_data)
  
  g_data1 <- g_data[1:500,]
  g_data2 <- g_data[501:1000,]
  g_data3 <- g_data[1001:1500,]
  g_data4 <- g_data[1501:2000,]
  
  ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
    geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE10_1_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
    geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE10_2_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
    geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE10_3_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
    geom_vline(xintercept=true_coeff[[15]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE10_4_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  
  ### POF-00 ###
  datg <- effect_interest[[6]][,tp]
  
  tp1_datg <- matrix(nrow = 2000, ncol = 3)
  
  for (i in 1:length(datg)) {
    tp1_datg[i,] <- c(datg[i],CI[[5]][[i]][tp,])
  }
  
  g_data <- as.data.frame(tp1_datg)
  g_data$SampleID <- rownames(g_data)
  
  g_data1 <- g_data[1:500,]
  g_data2 <- g_data[501:1000,]
  g_data3 <- g_data[1001:1500,]
  g_data4 <- g_data[1501:2000,]
  
  ggplot(data=g_data1, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1-500') +
    geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE00_1_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data2, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 501-1000') +
    geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE00_2_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data3, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1001-1500') +
    geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE00_3_",tp,".png"), width = 20, height = 12, units = "in")
  
  
  ggplot(data=g_data4, aes(y=SampleID, x=V1, xmin=V2, xmax=V3)) +
    geom_point(color='blue') +
    geom_errorbarh(height=.2) +
    labs(title='Effect Size by Study Sample', x='Effect Size', y = 'Sample ID 1501-2000') +
    geom_vline(xintercept=true_coeff[[16]][tp], color='red', linetype='dashed', alpha=1, linewidth=1) +
    theme_minimal() + coord_flip() +
    theme(axis.text.x = element_blank(),  #remove y axis labels
          axis.ticks.x = element_blank())
  ggsave(paste("CI_cov_IE00_4_",tp,".png"), width = 20, height = 12, units = "in")
  
}

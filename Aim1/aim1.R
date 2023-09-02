library(dplyr)
library(MASS)
library(VGAM)

set.seed(111722)
N <- 20000

j <- seq(from = 1, to = 20, length.out = 20)
## Coefficient curves
alpha0 <- -0.4
alpha1 <- exp((-0.25 -(2/j)))

beta0 <- -0.1
beta0_cont <- 2.3

gamma_A <- -0.1*exp(1/j)

beta1_c <- -(1/exp(1/(j^3)))
beta1_r <- -(1/exp(1/j))-1.9

treat  <- rbinom(N,1,0.5)
x1 <- rnorm(N,2.25,0.5)

dat1 <- as.data.frame(cbind(treat, x1))
dat1 <- tibble::rownames_to_column(dat1, "id")

dat2 <- dat1[rep(seq_len(nrow(dat1)), each = 20), ]
dat2 <- dat2 %>% 
  dplyr::group_by(id) %>% # build grouping by combination of variables
  dplyr::mutate(time_pt = dplyr::row_number()) %>% # add row number which works per group due to prior grouping
  dplyr::ungroup()  # ungroup to prevent unexpected behavior down stream

coeff <- as.data.frame(cbind(alpha0,
                             cbind(alpha1,
                                   cbind(beta0,
                                         cbind(beta0_cont,
                                               cbind(gamma_A,
                                                     cbind(beta1_c, beta1_r)))))))
coeff <- coeff %>%
  dplyr::mutate(time_pt = dplyr::row_number())

dat3 <- merge(dat2, coeff, by = "time_pt")
dat4 <- dat3 %>%
  arrange(as.numeric(id))
dat4$id <- as.numeric(dat4$id)

for(i in 1:nrow(dat4)){
  ## Create the time-varying co-variate data
  if(dat4$time_pt[i] == 1){
    dat4$prob_L[i] <- exp(-0.2*dat4$x1[i])/(1 + exp(-0.2*dat4$x1[i]))
  }else if(dat4$time_pt[i] > 1){
    dat4$prob_L[i] <- exp(-0.2*dat4$x1[i] + 0.25*dat4$L[i-1])/(1 + exp(-0.2*dat4$x1[i] + 0.25*dat4$L[i-1]))
  }
  dat4$L[i] <- rbinom(1,1,dat4$prob_L[i])

  ## Create mediator data
  if(dat4$time_pt[i] == 1){
    num <- exp(dat4$alpha0[i] + dat4$alpha1[i]*dat4$treat[i] - 1.2*dat4$L[i] + 0.05*dat4$x1[i])
  }else if(dat4$time_pt[i] > 1){
    num <- exp(dat4$alpha0[i] + dat4$alpha1[i]*dat4$treat[i] + 1.25*dat4$M[i-1] - 1.2*dat4$L[i] + 0.05*dat4$x1[i])
  }
  dat4$prob_M[i] <- num/(1 + num)
  dat4$M[i] <- rbinom(1,1,dat4$prob_M[i])
  print(paste("M_",i))
}

for(i in 1:nrow(dat4)){
  ## Create normally distributed continuous outcome data
  if(dat4$time_pt[i] == 1){
    dat4$mu.y[i] <- dat4$beta0_cont[i] + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_c[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i] + rnorm(1, 0, 1)
  }else if(dat4$time_pt[i] > 1){
    dat4$mu.y[i] <- dat4$beta0_cont[i] + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_c[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i] + 0.5*dat4$Ycrav[i-1] + rnorm(1, 0, 1)
  }
  dat4$Ycrav[i] <- rnorm(1, dat4$mu.y[i], 1)

  ## Create binary outcome data
  if(dat4$time_pt[i] == 1){
    dat4$num_c[i] <- (dat4$beta0[i]-0.74) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_c[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i]
    dat4$num_r[i] <- (dat4$beta0[i]-2) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_r[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i]
  }else if(dat4$time_pt[i] > 1){
    dat4$num_c[i] <- (dat4$beta0[i]-0.74) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_c[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i] + 0.25*dat4$Yc[i-1]
    dat4$num_r[i] <- (dat4$beta0[i]-2) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_r[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i] + 0.9*dat4$Yr[i-1]
  }
  dat4$prob_Yc[i] <- exp(dat4$num_c[i])
  dat4$prob_Yr[i] <- exp(dat4$num_r[i])/(1 + exp(dat4$num_r[i]))
  dat4$Yc[i] <- rbinom(1,1,dat4$prob_Yc[i])
  dat4$Yr[i] <- rbinom(1,1,dat4$prob_Yr[i])

  ## Create count outcome data
  if(dat4$time_pt[i] == 1){
    dat4$num_pr1[i] <- (dat4$beta0[i]-1.17) + 0.2*dat4$treat[i] + 0.9*dat4$M[i] - 0.36*dat4$L[i] + 0.02*dat4$x1[i]
    dat4$mu.count.zi[i] <- exp((dat4$beta0[i] + 0.5) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_c[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i])
    dat4$mu.count.nzi[i] <- exp((dat4$beta0[i] + 0.5) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_r[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i])
  }else if(dat4$time_pt[i] > 1){
    dat4$num_pr1[i] <- (dat4$beta0[i]-1.17) + 0.2*dat4$treat[i] + 0.9*dat4$M[i] - 0.36*dat4$L[i] + 0.02*dat4$x1[i] - 0.05*dat4$Cigc[i-1]
    dat4$mu.count.zi[i] <- exp((dat4$beta0[i] + 0.5) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_c[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i] + 0.05*dat4$Cigc[i-1])
    dat4$mu.count.nzi[i] <- exp((dat4$beta0[i] + 0.5) + dat4$gamma_A[i]*dat4$treat[i] + dat4$beta1_r[i]*dat4$M[i] + 0.36*dat4$L[i] - 0.02*dat4$x1[i] + 0.05*dat4$Cigr[i-1])
  }
  dat4$prob0_zi[i] <- exp(dat4$num_pr1[i])

  dat4$Cigc[i] <- rzipois(1, dat4$mu.count.zi[i], pstr0 = dat4$prob0_zi[i]) ## zero-inflated count ##
  dat4$Cigr[i] <- rzipois(1, dat4$mu.count.nzi[i])  ## not zero-inflated count ##
  print(paste("Y_",i))
}

aim1_dat <- dat4
save(aim1_dat, file = "aim1_dat.rda")

########hierarchical.sex########
##Léo Le Gall
##created in 15/05/2025
##last modified in 15/05/2025
################################
growth <- function(linf, k, t0, age) {
  linf * (1 - exp(-k * (age - t0)))
}

## Tau function
tau <- function(theta) {
  linfF <- theta[1]
  linfM <- theta[2]
  kF <- theta[3]
  kM <- theta[4]
  t0F <- -abs(theta[5])
  t0M <- -abs(theta[6])
  sdF <- exp(theta[7])
  sdM <- exp(theta[8])
  
  muF <- log(growth(linfF, kF, t0F, ages))
  muM <- log(growth(linfM, kM, t0M, ages))
  
  #density function in log
  densF0 <- sapply(1:m, function(z) {
    log(mixp_s[pos_y]) + log(mixp[pos_y,z]) + dlnorm(dat$l, meanlog = muF[z] - sdF^2 / 2, sdlog = sdF, log = TRUE)
  })
  densM0 <- sapply(1:m, function(z) {
    log(1 - mixp_s[pos_y]) + log(mixp[pos_y,z]) + dlnorm(dat$l, meanlog = muM[z] - sdM^2 / 2, sdlog = sdM, log = TRUE)
  })
  
  #!key point! increase the numbers for computation (subtract the max to the log-dens)
  row_max <- apply(cbind(densF0,densM0), 1, max)
  
  densF <- densF0 - row_max
  densM <- densM0 - row_max
  
  tauF1 <- exp(densF) / rowSums(cbind(exp(densF),exp(densM)))
  
  tauM1 <- exp(densM) / rowSums(cbind(exp(densF),exp(densM)))
  
  ## classified ind
  tauF1[-idxUU,] <- 0
  tauM1[-idxUU,] <- 0
  
  vecaF <- dat$a[idxF] + 1
  vecaM <- dat$a[idxM] + 1
  
  tauF1[cbind(idxF,vecaF)] <- 1
  tauM1[cbind(idxM,vecaM)] <- 1
  
  ##Unsexed ind (datU)
  vecaU <- dat$a[idxU] + 1
  
  tauU1 <- exp(densM) / (exp(densF) + exp(densM))
  
  tauF1[cbind(idxU,vecaU)] <- tauU1[cbind(idxU,vecaU)]
  tauM1[cbind(idxU,vecaU)] <- 1 - tauU1[cbind(idxU,vecaU)]
  
  # ret function
  ret<-list(Female = tauF1, Male = tauM1)
  
  return(ret)
}

cll <- function(theta) {
  linfF <- theta[1]
  linfM <- theta[2]
  kF <- theta[3]
  kM <- theta[4]
  t0F <- -abs(theta[5])
  t0M <- -abs(theta[6])
  sdF <- exp(theta[7])
  sdM <- exp(theta[8])
  
  muF <- log(growth(linfF, kF, t0F, ages))
  muM <- log(growth(linfM, kM, t0M, ages))
  
  densF<- sapply(1:m, function(z) {
    log(mixp_s[pos_y]) + log(mixp[pos_y,z]) + dlnorm(dat$l, meanlog = muF[z] - sdF^2 / 2, sdlog = sdF, log = TRUE)
  })
  densM <- sapply(1:m, function(z) {
    log(1 - mixp_s[pos_y]) + log(mixp[pos_y,z]) + dlnorm(dat$l, meanlog = muM[z] - sdM^2 / 2, sdlog = sdM, log = TRUE)
  })
  
  # complete log likelihood equation
  cllF <- sapply(1:m, function(z) {dat$n * tauF.cll[, z] * densF[, z]})
  cllM <- sapply(1:m, function(z) {dat$n * tauM.cll[, z] * densM[, z]})
  
  return(-sum(cllF) - sum(cllM))
}

oll <- function(theta) {
  linfF <- theta[1]
  linfM <- theta[2]
  kF <- theta[3]
  kM <- theta[4]
  t0F <- -abs(theta[5])
  t0M <- -abs(theta[6])
  sdF <- exp(theta[7])
  sdM <- exp(theta[8])
  
  muF <- log(growth(linfF, kF, t0F, ages))
  muM <- log(growth(linfM, kM, t0M, ages))
  
  densF <- sapply(1:m, function(z) {
    log(mixp_s[pos_y]) + log(mixp[pos_y,z]) + dlnorm(dat$l, meanlog = muF[z] - sdF^2 / 2, sdlog = sdF, log = TRUE)
  })
  densM <- sapply(1:m, function(z) {
    log(1 - mixp_s[pos_y]) + log(mixp[pos_y,z]) + dlnorm(dat$l, meanlog = muM[z] - sdM^2 / 2, sdlog = sdM, log = TRUE)
  })
  
  #classified individuals
  ollF <- sapply(1:m, function(z){dat$n * tauF.oll[,z] * densF[,z]})
  ollM <- sapply(1:m, function(z){dat$n * tauM.oll[,z] * densM[,z]})
  
  ##unsexed
  vecaU <- dat$a[idxU] + 1
  tauU<-matrix(data=0,nrow=nrow(dat),ncol=m)
  tauU[cbind(idxU,vecaU)] <- 1
  
  ollU <- sapply(1:m, function(z){dat$n * tauU[,z] * log(exp(densF[,z]) + exp(densM[,z]))})
  
  ##unsexed and unaged
  ollUU <- dat$n * log(rowSums(exp(densF) + exp(densM)))
  
  ##calculate the observed log-likelihood
  olltot <- sum(ollF[idxF,]) + sum(ollM[idxM,]) + sum(ollU[idxU,]) + sum(ollUU[idxUU])
  return(olltot)
}

##### EM #####
nloop <- 500
##set starting values
para <- c(rep(400,2), rep(0.4,2), rep(1,2), rep(log(0.1),2))
mixp <- matrix(rep(rep(1/m,m)),nrow= length(years), ncol=m)
mixp_s <- rep(rep(1/2,length(years)))
ollvec<-rep(NA,nloop)
maxiter<-1e-10

for (i in 1:nloop) {
  cat("Iteraction:", i, "\n")
  
  tau.cll <- tau(para)
  tauF.cll <- tau.cll$Female
  tauM.cll <- tau.cll$Male
  
  opt <- optim(
    para, cll, 
    method = "L-BFGS-B",
    lower = c(rep(100,2), rep(0.01,2), rep(0.1,2), rep(log(0.01),2)),
    upper = c(rep(2000,2), rep(1,2), rep(5,2), rep(log(1),2))
  )
  
  para <- opt$par
  
  tau.oll <- tau(para)
  tauF.oll <- tau.oll$Female
  tauM.oll <- tau.oll$Male
  
  ntauF <- sapply(1:m, function(z) {dat$n * tauF.oll[, z]})
  ntauM <- sapply(1:m, function(z) {dat$n * tauM.oll[, z]})
  
  ollvec[i] <- oll(para)
  
  ## stop the loop if stab
  if(i > 1){if(abs(ollvec[i]-ollvec[i-1])<ollvec[i]*maxiter){break}}
  
  ### update mix. prop.
  for (a in seq_along(years)){
    mixp[a,] <- colSums(rbind(ntauF[pos_y==a,],ntauM[pos_y==a,])) / sum(ntauF[pos_y==a,],ntauM[pos_y==a,])
    mixp_s[a] <- sum(ntauF[pos_y==a,])/sum(ntauF[pos_y==a,],ntauM[pos_y==a,])
  }
  ## plot
  plot(1:length(ollvec), ollvec, xlab = "Interaction number", ylab = "Observed log-likelihood")
}
dev.off()
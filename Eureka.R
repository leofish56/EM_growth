########Eureka2_log#######
##Léo Le Gall
##created in 09/05/2024
##last modified in 12/05/2025
####################
#### Functions####
growth <- function(linf, k, t0, age) {
  linf * (1 - exp(-k * (age - t0)))
}

## Tau function
tau <- function(theta) {
  linf <- theta[1]
  k <- theta[2]
  t0 <- -abs(theta[3])
  sd <- exp(theta[4])
  
  mu <- growth(linf, k, t0, ages)
  
  dens <- sapply(1:m, function(z) {
    log(mixp[z]) + dnorm(dat$l, mean = mu[z], sd = sd, log = TRUE)
  })
  
  dens <- (dens - apply(dens, 1, max))
  
  tau1 <- exp(dens) / rowSums(exp(dens))
  
  if (any(tau1 == 0)) {
    warning("Near-zero probabilities detected, indicating numerical instability.")
  }
  
  tau1[-idxUU, ] <- 0
  
  vecaF <- dat$a[idxF] + 1
  vecaM <- dat$a[idxM] + 1
  vecaU <- dat$a[idxU] + 1
  
  tau1[cbind(idxF, vecaF)] <- 1
  tau1[cbind(idxM, vecaM)] <- 1
  tau1[cbind(idxU, vecaU)] <- 1
  
  return(tau1)
}

cll <- function(theta) {
  linf <- theta[1]
  k <- theta[2]
  t0 <- -abs(theta[3])
  sd <- exp(theta[4])
  
  mu <- growth(linf, k, t0, ages)
  
  dens <- sapply(1:m, function(z) {
    log(mixp[z]) + dnorm(dat$l, mean = mu[z], sd = sd, log = TRUE)
  })
  
  cll <- sapply(1:m, function(z) {dat$n * tau.cll[, z] * dens[, z]})
  
  return(-sum(cll))
}

oll <- function(theta) {
  linf <- theta[1]
  k <- theta[2]
  t0 <- -abs(theta[3])
  sd <- exp(theta[4])
  
  mu <- growth(linf, k, t0, ages)
  
  dens <- sapply(1:m, function(z) {
    log(mixp[z]) + dnorm(dat$l, mean = mu[z], sd = sd, log = TRUE)
  })
  
  ollCA <- sapply(1:m, function(z) {dat$n * tau.oll[, z] * dens[, z]})
  
  ollUU <- dat$n * log(apply(exp(dens), 1, sum))
  
  olltot <- sum(c(ollCA[-idxUU, ], ollUU[idxUU]))
  return(olltot)
}

update_mixp <- function(ntau) {
  mixp <- colSums(ntau) / sum(ntau)
  if (any(mixp == 0)) {
    mixp <- mixp + 1e-100
  }
  mixp <- mixp / sum(mixp)
  return(mixp)
}

##### EM #####
para <- c(1000, 0.1,2, log(0.5))

##set starting mixing proportion
if(miss_CA){
  mixp <- rep(1/m,m)
}else{
  mixp <- rep(1e-100, m)
  prop <- prop.table(table(dat$a))
  vec_a <- unique(sort(dat$a))
  mixp[vec_a + 1]<- prop
}
mixp <- rep(1/m,m)
nloop <- 200
ollvec<-rep(NA,nloop)
maxiter<-1e-3
for (i in 1:nloop) {
  cat("Iteraction:", i, "\n")
  
  tau.cll <- tau(para)
  
  opt <- optim(
    para, cll, 
    method = "L-BFGS-B",
    lower = c(200, 0.1, 0.1, log(0.1)),
    upper = c(1000, 1.5, 10, log(20))
  )
  
  para <- opt$par
  
  tau.oll <- tau(para)
  ntau <- sapply(1:m, function(z) {dat$n * tau.oll[, z]})
  
  ollvec[i] <- oll(para)
  if(i > 1){if(abs(ollvec[i]-ollvec[i-1])<maxiter){break}}
  
  mixp <- update_mixp(ntau)
  
  growth_values <- growth(para[1], para[2], -abs(para[3]), 0:(m-1))
  point_sizes <- mixp * 10
  
  #plot(0:(m-1), growth_values, type="b", cex=point_sizes, pch=1, col="black",xlab="Age", ylab="Length (mm)")
  #plot(1:nloop, ollvec, xlab = "Interaction number", ylab = "Observed log-likelihood")
}


####final plot####
linf <- para[1]
k <- para[2]
t0 <- -abs(para[3])
sd <- exp(para[4])

# Aggregate length frequencies
lf <- aggregate(n ~ l, data = dat, FUN = sum)

# Calculate densities
mu <- growth(linf, k, t0, ages)

dens <- sapply(1:m, function(z) {
  dnorm(lf$l, mean = mu[z], sd = sd)
})

tau_p0 <- t(dens) / colSums(dens)
tau_p <- t(tau_p0)

total_n <- sum(lf$n)

tot_tau <- t(tau_p) * total_n * mixp
# Create barplot
bp <- barplot(height = lf$n, names.arg = lf$l, col = "grey", border = "white",
              xlab = "Length", ylab = "Number")
abline(h=0, lwd=2,col="white")
legend("topright", legend=years, bty="n")
# Overlay model components
lines(x = bp, y = colSums(tot_tau), lwd = 2, lty = 2)
for (i in 1:m) {
  pred <- tau_p[, i] * total_n * mixp[i]
  lines(x = bp, y = pred, lwd = 2)
  points(y=max(pred), x=bp[pred==max(pred)], pch= 20, cex =5)
  text(ages[i], y=max(pred), x=bp[pred==max(pred)], col="white")
}

########Eureka_log_coilin#######
##Léo Le Gall
##created in 09/05/2024
##last modified in 13/05/2025
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
  
  mu <- log(growth(linf, k, t0, ages))
  
  #density function in log
  dens <- sapply(1:m, function(z) {
    log(mixp[z]) + dlnorm(dat$l, meanlog = mu[z] - sd^2 / 2, sdlog = sd, log = TRUE)
  })
  
  #!key point! increase the numbers for computation (substract the max to the log-dens)
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
  
  mu <- log(growth(linf, k, t0, ages))
  
  dens <- sapply(1:m, function(z) {
    log(mixp[z]) + dlnorm(dat$l, meanlog = mu[z] - sd^2 / 2, sdlog = sd, log = TRUE)
  })
  
  cll <- sapply(1:m, function(z) {dat$n * tau.cll[, z] * dens[, z]})
  
  return(-sum(cll))
}

oll <- function(theta) {
  linf <- theta[1]
  k <- theta[2]
  t0 <- -abs(theta[3])
  sd <- exp(theta[4])
  
  mu <- log(growth(linf, k, t0, ages))
  
  dens <- sapply(1:m, function(z) {
    log(mixp[z]) + dlnorm(dat$l, meanlog = mu[z] - sd^2 / 2, sdlog = sd, log = TRUE)
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

##set starting mixing proportion
para <- c(500, 0.4,2, log(0.1))
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
    lower = c(50, 0.05, 0.05, log(0.05)),
    upper = c(1500, 0.8, 10, log(1))
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
  plot(1:nloop, ollvec, xlab = "Interaction number", ylab = "Observed log-likelihood")
}


####final plot####
linf <- para[1]
k <- para[2]
t0 <- -abs(para[3])
sd <- exp(para[4])

mu <- log(growth(linf, k, t0, ages))

# Aggregate length frequencies
lf <- aggregate(n ~ l, data = dat, FUN = sum)

# Calculate densities
seq <- seq(0,max(lf$l), by=10)
dens <- sapply(1:m, function(z) {
  dlnorm(seq, meanlog = mu[z] - sd^2 / 2, sdlog = sd)
})

tau_p0 <- t(dens) / colSums(dens)
tau_p <- t(tau_p0)

total_n <- sum(lf$n)

tot_tau <- t(tau_p) * total_n * mixp

# Create barplot
plot(lf$l, lf$n, type= "h", lwd=8, col="grey", xlab = "Lenght (mm)", ylab = "Number", ylim=c(0,max(lf$n)),yaxs="i")
points(lf$l, lf$n, pch=15, cex =1.5, col="white")
legend("topright", legend=years, bty="n")
box()


# Overlay model components
lines(x = seq, y = colSums(tot_tau), lwd = 2, lty = 2)
palette_10 <- c(
  "#88CCEE", "#44AA99", "#117733", "#DDCC77",
  "#CC6677", "#AA4499", "#882255", "#661100", "#6699CC", "#332288"
)
for (i in 1:m) {
  pred <- tau_p[, i] * total_n * mixp[i]
  polygon(x=seq, y= pred, col = adjustcolor(palette_10[i], alpha.f = 0.2))
  lines(x = seq, y = pred, lwd = 2, col=palette_10[i])
  points(y=max(pred), x=seq[pred==max(pred)], pch= 20, cex =5, col=palette_10[i])
  text(ages[i], y=max(pred), x=seq[pred==max(pred)], col="white")
}

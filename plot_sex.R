####final plot####
#jpeg("plotsex.jpeg", width = 800, height = 600, quality = 100)
linfF <- para[1]
linfM <- para[2]
kF <- para[3]
kM <- para[4]
t0F <- -abs(para[5])
t0M <- -abs(para[6])
sdF <- exp(para[7])
sdM <- sdF

muF <- log(growth(linfF, kF, t0F, ages))
muM <- log(growth(linfM, kM, t0M, ages))

# Aggregate length frequencies
lf <- aggregate(n ~ l, data = dat, FUN = sum)

# Calculate densities
seq <- seq(0,max(lf$l), by=10)
densF <- sapply(1:m, function(z) {
  dlnorm(seq, meanlog = muF[z] - sdF^2 / 2, sdlog = sdF)
})
densM <- sapply(1:m, function(z) {
  dlnorm(seq, meanlog = muM[z] - sdM^2 / 2, sdlog = sdM)
})

tauF_p0 <- t(densF) / colSums(densF)
tauF_p <- t(tauF_p0)
tauM_p0 <- t(densM) / colSums(densM)
tauM_p <- t(tauM_p0)

total_n <- sum(lf$n)

tot_tauF <- sapply(1:m, function(z){tauF_p[,z] * mixp[z] * mixp_s * total_n})
tot_tauM <- sapply(1:m, function(z){tauM_p[,z] * mixp[z] * (1-mixp_s) * total_n})
tot_tau <- tot_tauF + tot_tauM

# Create barplot
plot(lf$l, lf$n, type= "h", lwd=8, col="grey", xlab = "Lenght (mm)", ylab = "Number", ylim=c(0,max(lf$n)),yaxs="i")
points(lf$l, lf$n, pch=15, cex =1.5, col="white")
legend("topright", legend=years, bty="n")


# Overlay model components
lines(x = seq, y = rowSums(tot_tau), lwd = 2, lty = 2)
palette_10 <- c("pink", "turquoise")
for (i in 1:m) {
  predF <- tauF_p[, i] * total_n * mixp[i] * mixp_s
  predM <- tauM_p[, i] * total_n * mixp[i] * (1-mixp_s)
  
  polygon(x=seq, y= predF, col = adjustcolor(palette_10[1], alpha.f = 0.2))
  lines(x = seq, y = predF, lwd = 2, col=palette_10[1])
  points(y=max(predF), x=seq[predF==max(predF)], pch= 20, cex =5, col=palette_10[1])
  text(ages[i], y=max(predF), x=seq[predF==max(predF)], col="white")
  
  polygon(x=seq, y= predM, col = adjustcolor(palette_10[2], alpha.f = 0.2))
  lines(x = seq, y = predM, lwd = 2, col=palette_10[2])
  points(y=max(predM), x=seq[predM==max(predM)], pch= 20, cex =5, col=palette_10[2])
  text(ages[i], y=max(predM), x=seq[predM==max(predM)], col="white")
}
box()
#dev.off()

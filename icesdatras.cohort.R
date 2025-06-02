########icesdatras.year#######
##Léo Le Gall
##created in 14/10/2024
##last modified in 12/05/2025
################################
#install.packages("icesDatras")
setwd("C:/Users/legall/Documents/Coilin/code")
library(icesDatras)
library(nnet)
# select the years (2003 to 2023)
years<-2003
spHL<-data.frame()
spCA<-data.frame()
for (i in seq_along(years)){
  print(years[i])
  HLdata <- getHLdata(survey = "IE-IGFS", year = years[i], quarter = 4)
  CAdata <- getCAdata(survey = "IE-IGFS", year = years[i], quarter = 4)
  ##select one specie (ex: 126439=WHB; 126437=HAD; 126417=HER; 126438=MER; 127143=Plaice; 126822=H MAC; 126436 = Cod)
  spHL<-rbind(spHL,subset(HLdata, SpecCode==126437))
  spCA<-rbind(spCA,subset(CAdata, SpecCode==126437))
}
HL <- spHL
CA <- spCA
miss_CA <- all(is.na(CA$Age)) || length(CA$Age) == 0
if (miss_CA) {
  ages <- 0:5
  warning("Warning: no age information in CA")
} else {
  ages <- seq(0, max(CA$Age, na.rm = TRUE))
}
m <- length(ages)
## set functions
growth <- function(linf, k, t0, age) {
  linf * (1 - exp(-k * (age - t0)))
}

## raised number at length 
HL$n <- HL$HLNoAtLngt * HL$SubFactor

## aggregate across hauls
nl <- aggregate(n ~ LngtClass + Year, FUN = sum, data = HL)

names(nl)[names(nl) == "Year"] <- "y"
names(nl)[names(nl) == "LngtClass"] <- "l"
ynl <- nl$y

## aged data
CA$n <- CA$CANoAtLngt

if(miss_CA){
  na <- aggregate(n ~ LngtClass + Sex + Year, FUN = sum, data = CA)
  na$a <- NA
}else{
  na <- aggregate(n ~ LngtClass + Age + Sex + Year, FUN = sum, data = CA)
}
names(na)[names(na) == "LngtClass"] <- "l"
names(na)[names(na) == "Age"] <- "a"
names(na)[names(na) == "Sex"] <- "s"
names(na)[names(na) == "Year"] <- "y"

## remove the aged lengths from the unaged
nla <- merge(nl, aggregate(n ~ l + y, FUN = sum, data = na), by = c("y","l"), all.x = TRUE)
nla[is.na(nla)] <- 0
nla$n <- NA
nla$n <- with(nla, ifelse(n.x > n.y, n.x - n.y, n.y))
nl <- nla[, c("l","y","n")]
nl$a <- NA
nl$s <- NA

na$source <- "age"
nl$source <- "length"

dat <- rbind(na, nl)

dat$source <- factor(dat$source, levels = c("length", "age"))
idxF <- which(dat$s == "F")
idxM <- which(dat$s == "M")
idxU <- which(dat$s == "U")
idxUU <- which(dat$source == "length")
idxCA <- which(dat$source == "age")

pos_y <- match(dat$y, years)
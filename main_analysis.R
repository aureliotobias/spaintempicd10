###########################################################################################################
###########################################################################################################
### Tobías A, Íñiguez C, Royé D, Hashizume M, Madaniyazi L.
### Cause-specific mortality burden attributable to ambient temperature and seasonal variations in Spain. 
### To be submitted.
###########################################################################################################
### REPRODUCIBLE EXAMPLE AS IN THE MANUSCRIPT-
### Last update: 2024.08.07
### main_analysis.R
###########################################################################################################
###########################################################################################################

# Remove objects.
rm(list = ls())

# Load libraries.
library(foreign) 
library(mgcv) 
library(tsModel)
library(dlnm)
library(splines)
library(mondate)
library(ggplot2)
library(devtools)

# Load functions.
source("cyclic.R")
source("attrs.R")
source("attrdl.R")

###########################################################################################################
###  Load Valencia dataset
###########################################################################################################

# Load and inspect the dataset.
data <- read.csv('valencia.csv')
names(data)
head(data)

# Formatting date.
date_string <- data$date
date_format <- "%d%b%Y"
data$date <- as.Date(date_string, format=date_format)
head(data)

###########################################################################################################
### Data management                                                    
###########################################################################################################

# Generating time variables.
data$year <- as.factor(as.character(data$year))
years <- unique(data$year)
data$dow <- as.factor(as.character(data$dow))

# Generating stratum with year and dow to control for long term trend and effect of day of week.
data$stratum <- as.factor(data$year:data$dow)

# Generating seasonality indicator (day of year).
for(t in years){
  tempind <- data$year == t
  temp <- as.numeric( strftime( data$date[tempind] , format = "%j") )
  if( length(temp) == 365 ){ temp[60:365] <- 61:366 }
  else if( length(temp) == 364 ){ temp[60:364] <- 61:365 }
  
  data$seasonal[tempind] <- temp
}

# Derive the cyclic spline for seasonality.
spline.season <- cyclic(data$seasonal, df=4)

# Crossbasis for temperature.
cb.temp <- crossbasis(data$tmean, lag=21,
                      argvar=list(fun="ns", knots=quantile(data$tmean, c(.10,.75,.90), na.rm=T)) ,
                      arglag=list(fun="ns", knots=logknots(21,3)))

###########################################################################################################
### Fit the time-series regression model                                                  
###########################################################################################################

model <- glm(all ~ cb.temp + spline.season + factor(dow) + factor(stratum), 
             data=data, family=quasipoisson(), na.action="na.exclude")

###########################################################################################################
### Temperature-response
###########################################################################################################

# Get prediction centered at the MMT.
pred  <- crosspred(cb.temp, model, by=1) 
mmt <- pred$predvar[which.min(pred$allRRfit)] 
predcen <- crosspred(cb.temp, model, cen=mmt, by=.1) 

# Figure.
par(mex=0.8,mfrow=c(1,2))
  
###### Plot temperature-response.
dt <- cbind(predcen$predvar, predcen$allRRfit, predcen$allRRlow, predcen$allRRhigh)

# Heat.
dt.heat <- subset(dt, dt[,1] > mmt)
plot(dt.heat[,1], dt.heat[,2], type="l", lwd = 2,
     main="Temperature-response", 
     ylim=c(0.9,2), xlim=c(0,32), 
     ylab="Relative Risk", xlab="Temperature (ºC)",
     col="red")

polygon(c(dt.heat[,1], rev(dt.heat[,1])), c(dt.heat[,4], rev(dt.heat[,3])), 
        col=rgb(1, 0, 0, 0.2), border=NA)
#Cold.
dt.cold <- subset(dt, dt[,1] <= mmt)
lines(dt.cold[,1], dt.cold[,2], col="blue", lwd = 2)
polygon(c(dt.cold[,1], rev(dt.cold[,1])), c(dt.cold[,4], rev(dt.cold[,3])), 
        col=rgb(0, 0, 1, 0.2), border=NA)

# Reference lines.
abline(h=1, col="black")
abline(v=c(mmt), lty=(c(2)), col="grey")
box(lty = 1)

# Legend.
legend("topright", legend = c("Cold", "Heat"), 
       col=c("blue", "red"), lty=1, lwd=2, bty="n", cex=.8)

###### Attributable Deaths. 
min <- min(data$tmean, na.rm=T)
max <- max(data$tmean, na.rm=T)  

an.temp <- attrdl(data$tmean, cb.temp, data$all, model, cen=mmt, dir="forw", type="an", range=c(min,max))
simaf <- attrdl(data$tmean, cb.temp, data$all, model, cen=mmt, dir="forw", type="an", range=c(min,max), tot=T, sim=T)
an.low.temp <- quantile(simaf,c(2.5)/100)
an.upp.temp <- quantile(simaf,c(97.5)/100)
cbind(an.temp, an.low.temp, an.upp.temp)

###########################################################################################################
### Seasonality.
###########################################################################################################

# Obtaining coefficients from cyclic spline.
coef <- coef(model)[ 2:5 ]
vcov <- vcov(model)[ 2:5, 2:5 ]
  
# Obtaining coefficients for each day of year (logrr, se).
logrr <- matrix(NA, 1, 366 )
logrr_se <- matrix(NA, 1, 366 )

Trange <- 1:366
bvar <- cyclic(data$seasonal, df=4)
  
# logrr.
l <- length(Trange)
for(j in 1:l){ logrr[1,j] <- t(bvar[j,]) %*% as.numeric(coef) }
  
# standard error for logrr.
for(j in 1:l){ logrr_se[1,j] <- sqrt(as.numeric( t(bvar[j,]) %*% vcov%*% (bvar[j,]) ) ) }
  
# Centering logrr at trough.
maxcen <- apply(logrr, 1, which.min)
bvar_cen <- bvar[maxcen,]
  
# Centered logrr.
l <- length(Trange)
for(j in 1:l){ logrr[1,j] <- as.numeric(t(bvar[j,]-bvar_cen) %*% as.numeric(coef)) }
  
# sd for centered logrr.
for(j in 1:l){ logrr_se[1,j] <- sqrt( as.numeric(t(bvar[j,]-bvar_cen) %*% vcov%*% (bvar[j,]-bvar_cen)) ) }

###### Plot seasonality.

# Create day of year for x axis.
data$date <- as.Date(data$date, "%Y-/%m-/%d")
DY <- as.Date(data$date[data$year == 2004])
dateY <- format(DY, format="%b %d")
M <- mondate("1-1-2004")
FDM <- as.Date( rev(M - 1:12) )
FDM <- c( FDM, as.Date( DY[length(DY)] ) )
dateM <- format(FDM, format="%b %d")
monIND <- which( as.character(dateY) %in% as.character(dateM) )
  
# Make the plot.
plot(Trange, exp(logrr), type='l', xaxt='n',
     main="Seasonality",
     xlab="Day of year", ylab="Relative Risk", 
     ylim=c(0.9,2), 
     lwd=2, 
      col="black")
axis(1, at=Trange[monIND], labels=FALSE)
text(x=Trange[monIND], y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), 
     labels=dateM, srt=45, adj=1, xpd=TRUE, cex=0.8)
polygon(c(Trange, rev(Trange)), 
          c(exp(logrr-1.96*logrr_se), rev(exp(logrr+1.96*logrr_se)) ),
          col=rgb(0.128,0.128,0.128,alpha=0.1), border=NA)

# Reference line.
abline(h=1, col="black")

# Legend.
legend("topright", legend=c("Seasonal variation"), 
       col=c("black"), lty=1, lwd=2, bty="n", cex=.8)

# Close figure.
layout(1)

####### Attributable Deaths.
data$death <- data$all
an.season <- attrs(data$seasonal, spline.season, data, model, type="an", tot=T) 
simaf2 <- attrs(data$seasonal, spline.season, data, model, type="an", tot=T, sim=T)
an.low.season <- quantile(simaf2,c(2.5)/100)
an.upp.season <- quantile(simaf2,c(97.5)/100)
cbind(an.season, an.low.season, an.upp.season)

###########################################################################################################
###########################################################################################################
###                                       End of script file 
###########################################################################################################
###########################################################################################################


#########################################
###  Script for Dominick Retail Data
#########################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\oneDrive\\Documents\\GitHub\\sigex")
load_all(".")

######################
### Part I: load data

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting")
dom <- read.table("bullwhipCase.dat")

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting\\Figures")


#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

# tti is Bathroom tissues (product 65)
# ptw is Paper towels (product 41)
# No information on start date, somewhere between 1989 and 1994
start.date <- c(1,1)
period <- 1

## create ts object and plot
dataALL.ts <- sigex.load(dom,start.date,period,c("tti","ptw"),TRUE)

#############################
## select span and transforms 

N <- dim(dataALL.ts)[2]
T <- dim(dataALL.ts)[1]
## replace meager values by NA  
for(k in 1:N)
{
  times.zero <- seq(1,T)[dataALL.ts[,k] <= 0]
  times.zero <- times.zero[-seq(1,length(times.zero))[is.na(times.zero)]]
  dataALL.ts[times.zero,k] <- NA
}

## full span with log transform
transform <- "log"
aggregate <- FALSE
subseries <- c(1,2)
begin.date <- start.date
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)



###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"varma",c(3,0),NULL,"process",1)
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)


##################################
### PART IV: Model Fitting

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting: commented out, this takes a while 
#fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## input parameter from previous fit (MLE on entire span)
#  divergence:  -317.8978
psi.mle <- c(0.00774525185660389, -0.98822758684062, -1.82365193777049, 
             0.157588913641002, -0.0493570539692489, 0.0298417041626058, 
             0.513614768414021, 0.101277279570936, 0.0194128259398051, 0.0805986605317239, 
             0.268243218791683, 0.0754100600644741, 0.0107627589935598, 0.0145065018011489, 
             0.0911281986829464, 4.12465264379121, 3.75849659831345)
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## get fixed effects
reg.trend <- NULL
for(k in 1:N) {
  reg.trend <- cbind(reg.trend,sigex.fixed(data.ts,mdl,k,par.mle,"Trend"))
}


############################################
## checks for ragged lik calculation
#  (this block can be skipped for final results)
#  divergence should be -317.8978
data.vec <- matrix(t(data.ts - reg.trend),ncol=1)
na.flag <- seq(1,length(data.vec))[is.na(data.vec)==FALSE]
x.acf <- sigex.resid(psi.mle,mdl,data.ts)[[2]]
x.acf2 <- aperm(x.acf,c(1,3,2))
x.acf2[,,1] <- x.acf[,1,]/2
d <- 0
gamma.bigmat <- matrix(sigex.blocktoep(x.acf2[,,1:(T-d),drop=FALSE]),N*(T-d),N*(T-d))
gamma.bigmat <- gamma.bigmat + t(gamma.bigmat)
likvar.mat <- gamma.bigmat[na.flag,na.flag,drop=FALSE]
likvar.chol <- t(chol(likvar.mat))
x.eps <- solve(likvar.chol,data.vec[na.flag,drop=FALSE])
alt.lik <- 2*log(det(likvar.chol)) + t(x.eps) %*% x.eps
print(alt.lik)
print(sigex.lik(psi.mle,mdl,data.ts))

############################################

##  model checking 
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(Re(resid.mle)),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
sigex.portmanteau(resid.mle,80,length(psi.mle))
sigex.gausscheck(resid.mle)
acf(resid.mle,lag.max=80)

# bundle for default span
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)


##########################################
### Part V: Missing Value Imputation  
 
## load up the fitted model  
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## get imputations
tol <- 1e-10
times.na <- NULL
for(k in 1:N) { times.na <- union(times.na,seq(1,T)[is.na(data.ts)[,k]]) }
times.na <- sort(times.na)
data.casts <- sigex.midcast(psi,mdl,data.ts,0)
x.casts <- data.casts[[1]] + t(reg.trend[times.na,])
se.casts <- matrix(sqrt(tol+diag(data.casts[[2]])),nrow=N)
data.casts <- list()
data.casts[[1]] <- data.ts
data.casts[[1]][times.na,] <- t(x.casts)
data.casts[[2]] <- data.ts
data.casts[[2]][times.na,] <- t(x.casts) + 2*t(se.casts)
data.casts[[3]] <- data.ts
data.casts[[3]][times.na,] <- t(x.casts) - 2*t(se.casts)

# plot of series with imputations as dashed lines
#pdf(file="DomCasts.pdf")
par(mfrow=c(N,1))
for(k in 1:N) 
{
  plot(data.casts[[1]][,k],ylab=colnames(data.ts)[k],col=grey(.7))
  lines(data.casts[[2]][,k],col=grey(.7),lty=2)
  lines(data.casts[[3]][,k],col=grey(.7),lty=2)
  lines(data.ts[,k])
}
dev.off()


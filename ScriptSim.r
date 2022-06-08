########################################
##  Script for Casting Paper Simulations
########################################


## wipe
rm(list=ls())

library(devtools)
library(tictoc)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\sigex")
load_all(".")

###############################
### Part I: Simulation settings

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting\\Figures")

#set.seed(1234)

# define process
phi.matrix <- rbind(c(1,.5),c(-.2,.3))
innovar.matrix <- diag(2)

# user settings
#na_prop <- 0
na_prop <- .1
#na_prop <- .2
#na_prop <- .3
#na_prop <- .4
#na_prop <- .5
T <- 50
#T <- 100
#T <- 200
#T <- 400
#T <- 800
tol <- 1e-10
Monte <- 10000

############################
### Part II: Run Simulations


# compute model
data.ts <- ts(matrix(0,nrow=T,ncol=2))
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,2),"varma",c(1,0),NULL,"process",1)
mdl <- sigex.meaninit(mdl,data.ts,0)
par.true <- sigex.default(mdl,data.ts,NULL)
psi.true <- sigex.par2psi(par.true,mdl)
psi.true[4:7] <- var2.par2pre(array(phi.matrix,c(2,2,1)))
#sigex.psi2par(psi.true,mdl,data.ts)
gamma <- VARMAauto(array(phi.matrix,c(2,2,1)),NULL,innovar.matrix,10)
gamma.0 <- gamma[,,1]

times <- NULL
mses_na <- NULL
# run simulation
for(i in 1:Monte)
{

  na_inds <- sort(sample(seq(1,2*T),size=ceiling(na_prop*2*T),replace=FALSE))
  x.init <- t(chol(gamma.0)) %*% rnorm(2)
  x.next <- x.init
  x.sim <- NULL
  for(t in 1:T)
  {
    x.next <- phi.matrix %*% x.next + rnorm(2)
    x.sim <- cbind(x.sim,x.next)
  }
  x.sim_na <- matrix(t(x.sim),ncol=1)
  x.sim <- ts(t(x.sim),start=c(1,1))
  if(length(na_inds)>0) x.sim_na[na_inds] <- NA
  x.sim_na <- ts(matrix(x.sim_na,ncol=2),start=c(1,1))
  times.na <- NULL
  for(k in 1:2) { times.na <- union(times.na,seq(1,T)[is.na(x.sim_na)[,k]]) }
  times.na <- sort(times.na)
#  plot(x.sim,col=1)
#  lines(x.sim_na,col=2)
  
  tic()  
  likval <- sigex.lik(psi.true,mdl,x.sim,debug=FALSE)
  time <- toc(quiet=TRUE)
  times <- c(times,time$toc)
  
  data.casts <- sigex.midcast(psi.true,mdl,x.sim_na,0)
  x.casts <- data.casts[[1]]
  se.casts <- matrix(sqrt(tol+diag(data.casts[[2]])),nrow=2)
  data.casts <- list()
  data.casts[[1]] <- x.sim_na
  if(length(na_inds)>0) { data.casts[[1]][times.na,] <- t(x.casts) }
#  data.casts[[2]] <- x.sim_na
#  data.casts[[2]][times.na,] <- t(x.casts) + 2*t(se.casts)
#  data.casts[[3]] <- x.sim_na
#  data.casts[[3]][times.na,] <- t(x.casts) - 2*t(se.casts)
  
  mse_na <- sum((data.casts[[1]] - x.sim)^2)
  if(length(na_inds)>0) { mse_na <- mse_na/length(na_inds) }
  mses_na <- c(mses_na,mse_na)
  
  if(i %% 100 == 0) print(i)
}

mean(times)*10^{-6}
mean(mses_na)


#######################################
###  Script for all Daily Retail Data 
#	 used for "Casting" project
#######################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\oneDrive\\Documents\\GitHub\\sigex")
load_all(".")

######################
### Part I: load data

load("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SigExNew\\retail.RData")

#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date <- c(10,1,2012)
end.date <- day2date(dim(retail)[1]-1,start.date)
period <- 365

# calendar calculations
start.day <- date2day(start.date[1],start.date[2],start.date[3])
end.day <- date2day(end.date[1],end.date[2],end.date[3])
begin <- c(start.date[3],start.day) 
end <- c(end.date[3],end.day)

## create ts object and plot
dataALL.ts <- sigex.load(retail,begin,period,colnames(retail),TRUE)

#############################
## select span and transforms 

## series 4482 with no transform
transform <- "none"
aggregate <- FALSE
subseries <- 5
range <- NULL
dataONE.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)


#######################
## spectral exploratory

## levels
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
  sigex.specar(dataONE.ts,FALSE,i,7)
}
dev.off()

## growth rates
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
  sigex.specar(dataONE.ts,TRUE,i,7)
}
dev.off()


###########################
## embed as a weekly series

first.day <- 1
data.ts <- sigex.daily2weekly(dataONE.ts,first.day,start.date)
plot(data.ts)


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

############################## 
## Generate holiday regressors

easter.dates <- read.table("data\\easter500.txt")
easter.reg1 <- gethol(easter.dates,0,0,start.date,end.date)
easter.reg2 <- gethol(easter.dates,8,-1,start.date,end.date)

cyber.dates <- read.table("data\\cyber400.txt")
cyber.reg <- gethol(cyber.dates,0,0,start.date,end.date)

black.dates <- read.table("data\\black400.txt")
black.reg <- gethol(black.dates,0,0,start.date,end.date)

super.dates <- read.table("data\\super400.txt")
super.reg <- gethol(super.dates,0,0,start.date,end.date)

labor.dates <- read.table("data\\labor400.txt")
labor.reg <- gethol(labor.dates,0,0,start.date,end.date)

cny.dates <- read.table("data\\cny200.txt")
cny.reg <- gethol(cny.dates,0,0,start.date,end.date)

###########################
## Embed holiday regressors

easter.reg1 <- sigex.daily2weekly(easter.reg1,first.day,start.date)
easter.reg2 <- sigex.daily2weekly(easter.reg2,first.day,start.date)
cyber.reg <- sigex.daily2weekly(cyber.reg,first.day,start.date)
black.reg <- sigex.daily2weekly(black.reg,first.day,start.date)
super.reg <- sigex.daily2weekly(super.reg,first.day,start.date)
labor.reg <- sigex.daily2weekly(labor.reg,first.day,start.date)
cny.reg <- sigex.daily2weekly(cny.reg,first.day,start.date)

# replace ragged NA with zero
easter.reg1[is.na(easter.reg1)] <- 0
easter.reg2[is.na(easter.reg2)] <- 0
cyber.reg[is.na(cyber.reg)] <- 0
black.reg[is.na(black.reg)] <- 0
super.reg[is.na(super.reg)] <- 0
labor.reg[is.na(labor.reg)] <- 0
cny.reg[is.na(cny.reg)] <- 0


##############
## Basic Model

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"svarma",c(1,0,1,0,52),NULL,"process",1)
mdl <- sigex.meaninit(mdl,data.ts,0)	

for(i in 1:N) {
  mdl <- sigex.reg(mdl,i,ts(as.matrix(easter.reg1[,i]),
                            start=start(easter.reg1),
                            frequency=frequency(easter.reg1),
                            names="Easter-day"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(easter.reg2[,i]),
                            start=start(easter.reg2),
                            frequency=frequency(easter.reg2),
                            names="Easter-pre"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(cyber.reg[,i]),
                            start=start(cyber.reg),
                            frequency=frequency(cyber.reg),
                            names="Cyber"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(black.reg[,i]),
                            start=start(black.reg),
                            frequency=frequency(black.reg),
                            names="Black"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(super.reg[,i]),
                            start=start(super.reg),
                            frequency=frequency(super.reg),
                            names="Super"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(labor.reg[,i]),
                            start=start(labor.reg),
                            frequency=frequency(labor.reg),
                            names="Labor"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(cny.reg[,i]),
                            start=start(cny.reg),
                            frequency=frequency(cny.reg),
                            names="CNY"))
}


##################################
### PART IV: Model Fitting

#constraint <- cbind(rep(0,132), diag(133)[-1,])
constraint <- NULL
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(2,2,2,2,2,2,2),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(3,3,3,3,3,3,3),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(4,4,4,4,4,4,4),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(5,5,5,5,5,5,5),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(6,6,6,6,6,6,6),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(7,7,7,7,7,7,7),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(8,8,8,8,8,8,8),NULL))

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting: commented out, this took 3 weeks!
#fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
#psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint) 
#hess <- fit.mle[[1]]$hessian
#par.mle <- fit.mle[[2]]


## input parameter from previous fit (MLE on entire span)
#  divergence:    -4296.088 lik
psi.mle <- c(0.363244177168849, -0.0348124553411078, 0.491812935546558, 
             0.298664579485942, 0.958839706417251, 0.649058647979278, 0.608776644746655, 
             0.478831183040553, 0.129174845416729, 0.357441000025807, -0.175296271876546, 
             0.276005377856337, 0.269077062659677, 0.305839799443437, 0.465325619348851, 
             0.309989084472692, 1.00266650620707, 0.945284838668533, -0.351837870657996, 
             0.41272938717394, 0.591105620527945, -5.28994657588153, -4.91910496349852, 
             -4.27036055397041, -4.82248752232916, -4.65943436732755, -3.87588417776183, 
             -3.92941606603775, 
             0.243727358314956, 0.0451694016769132, 0.0895626865293538, 
               0.162027670053483, 0.097924889690983, 0.362718824863377, 0.154661073532565, 
               0.148850357971705, 0.377687048488826, 0.576611082980269, 0.0502380240838354, 
               0.0556123667223894, 0.0387204203264784, 0.465568884949394, 0.195482242662665, 
               0.171665647586408, 0.34261570013354, 0.475569362113592, 0.240527708472822, 
               0.23292421775797, 0.143311259387335, 0.216106756080234, 0.0826450255875768, 
               -0.418938435334888, -0.132447487084855, 0.12502877549537, -0.616449163936332, 
               -0.407399804578616, 0.239130444756953, 0.317841467345823, 0.11784445553256, 
               0.0443702198243786, 0.414400471006841, -0.000821672769191997, 
               -0.00591846616264505, 0.31081818878246, 0.368737311049801, 0.0251387704788187, 
               0.15218901786631, -0.00317463573051935, 0.231341416230008, -0.13897051162688, 
               0.560683226480955, 0.482827480083078, 0.592262387170034, 0.620213187459803, 
               0.483963576216744, 0.67673170378826, 1.15995146261087, 0.508079038102276, 
               0.230167938429291, 0.330870363588793, 0.0103709750110208, -0.177376178690009, 
               -0.0305357378118018, 0.140757998889117, 0.162124168372868, 0.731310717012488, 
               0.424920156388882, 0.218116682409072, 0.164674822291419, -0.241357660845099, 
               -0.377736226451398, -0.00991464966155541, -0.241091345031796, 
               0.167748645279208, 0.991028050342076, 0.0965846360554733, -0.215247547586599, 
               0.058431055194623, 0.0154142068904882, -0.0530377951073669, -0.208251834324292, 
               -0.235115963110736, 1.33842548154317, -0.252227481164646, 0.143774065072053, 
               -0.0435454492555431, -0.0161572013865816, 0.194520674009415, 
               -0.0994169729831187, -0.0961572958752303, 0.863454725859323, 
               -0.0540881437240361, 0.067578728713169, 0.381984426294581, 0.306580399438896, 
               0.444867384143703, -0.0389549582087286, 0.58728269678549, 0.417476853573103, 
               -0.108337739770599, -0.283307949032328, -0.222222674555547, -0.267920620663949, 
               -0.0630734590904744, -0.0729714585645782, 
               0.322687034086722,
             0.950787256576206, -0.852207175545125, 0.0584615417056083, 0.0532149301262836, 
             2.70711727317259, -0.118431624115301, 0.171386079598921, 0.0380331425134629, 
             0.938879481938871, -0.852207175545123, 0.0584615417056076, 0.0532149301262835, 
             2.70711727317259, -0.118431624115301, 0.17138607959892, 0.0380331425134629, 
             0.919716685878865, -0.852207175545125, 0.0584615417056085, 0.0532149301262835, 
             2.70711727317259, -0.118431624115301, 0.171386079598921, 0.038033142513463, 
             0.954766752731184, -0.852207175545125, 0.0584615417056093, 0.053214930126284, 
             2.70711727317259, -0.118431624115301, 0.171386079598921, 0.0380331425134627, 
             1.01434200276068, -0.852207175545126, 0.0584615417056081, 0.0532149301262839, 
             2.70711727317259, -0.118431624115301, 0.171386079598921, 0.0380331425134631, 
             1.31309017375681, -0.852207175545125, 0.0584615417056079, 0.0532149301262842, 
             2.70711727317259, -0.1184316241153, 0.171386079598921, 0.0380331425134629, 
             1.73161669258616, -0.852207175545125, 0.0584615417056081, 0.0532149301262825, 
             2.70711727317259, -0.118431624115301, 0.17138607959892, 0.0380331425134629)
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

##  model checking 
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(Re(resid.mle)),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=2*53,plot=FALSE)$acf

#pdf(file="retResidAcf.pdf",height=10,width=10)
par(mfrow=c(N,N),mar=c(3,2,2,0)+0.1,cex.lab=.8,cex.axis=.5,bty="n")
for(j in 1:N)
{
  for(k in 1:N)
  {
    plot.ts(resid.acf[,j,k],ylab="",xlab="Lag",ylim=c(-1,1),cex=.5)
    abline(h=1.96/sqrt(T),lty=3)
    abline(h=-1.96/sqrt(T),lty=3)
  }
}
#dev.off()

# bundle 
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)



##########################################
### Part V: Signal Extraction  

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting\\Figures")

## load up the fitted model for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)


## embed daily low-pass filter as a weekly filter
#mu <- 2*pi/365
#len <- 50000
#mu <- pi/7
#len <- 5000
#lp.hifilter <- c(mu/pi,sin(seq(1,len)*mu)/(pi*seq(1,len)))
#lp.hifilter <- c(rev(lp.hifilter),lp.hifilter[-1])

## embed daily SA filter as a weekly filter
sa.hifilter <- c(1,rep(2,365),1)/(2*365)
len <- 183
hi.freq <- 7
low.freq <- 1
shift.hi <- len
out <- sigex.hi2low(sa.hifilter,hi.freq,low.freq,shift.hi)
sa.lowfilter <- out[[1]]  
shift.low <- out[[2]]

sa.low <- sigex.adhocextract(psi,mdl,data.ts,sa.lowfilter,shift.low,0,TRUE)
sa.low.daily <- list()
sa.low.daily[[1]] <- sigex.weekly2daily(ts(sa.low[[1]],start=start(data.ts),frequency=frequency(data.ts)),first.day)
sa.low.daily[[2]] <- sigex.weekly2daily(ts(sa.low[[2]],start=start(data.ts),frequency=frequency(data.ts)),first.day)
sa.low.daily[[3]] <- sigex.weekly2daily(ts(sa.low[[3]],start=start(data.ts),frequency=frequency(data.ts)),first.day)

## embed daily TD filter as a weekly filter
#td.hifilter <- c(1,rep(2,7),1)/(2*7)
td.hifilter <- rep(1,7)/7
len <- 3
hi.freq <- 7
low.freq <- 1
shift.hi <- len
out <- sigex.hi2low(td.hifilter,hi.freq,low.freq,shift.hi)
td.lowfilter <- out[[1]]  
shift.low <- out[[2]]

td.low <- sigex.adhocextract(psi,mdl,data.ts,td.lowfilter,shift.low,0,TRUE)
td.low.daily <- list()
td.low.daily[[1]] <- sigex.weekly2daily(ts(td.low[[1]],start=start(data.ts),frequency=frequency(data.ts)),first.day)
td.low.daily[[2]] <- sigex.weekly2daily(ts(td.low[[2]],start=start(data.ts),frequency=frequency(data.ts)),first.day)
td.low.daily[[3]] <- sigex.weekly2daily(ts(td.low[[3]],start=start(data.ts),frequency=frequency(data.ts)),first.day)


## get fixed effects
reg.trend <- NULL
reg.trend <- cbind(reg.trend,mean(matrix(param[[4]],ncol =N)[1,])*rep(1,length(sa.low.daily[[1]])))

## plotting
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

#pdf(file="retSignals.pdf",height=8,width=10)
plot(dataONE.ts,xlab="Year")
sigex.graph(sa.low.daily,reg.trend,start(sa.low.daily[[1]]),
            period,1,0,trendcol,fade)
sigex.graph(td.low.daily,reg.trend,start(td.low.daily[[1]]),
            period,1,0,sacol,fade)
#dev.off()

## spectral diagnostics: seasonal adjustment
sigex.specar(sa.low.daily[[1]],FALSE,1,period)
#dev.off()

## spectral diagnostics: non-weekly effect
sigex.specar(td.low.daily[[1]],FALSE,1,7)
#dev.off()

 
 
   
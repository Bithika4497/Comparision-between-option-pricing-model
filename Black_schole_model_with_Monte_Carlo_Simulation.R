
#Comparision between the actual Black-Scholes option price model 
#with Monte Carlo simulated option price model
 
#initialize parameters

s=100   #initial stock price     
k=60    #strike price
r=0.05  #risk-free interest rate
T=1     #Maturity period
sig=0.4 #Risk

#create a stock price path by Geometric Brownian Motion

set.seed(1)
library(sde)
P<-GBM(x=100,r=0,N=3000) 

#Pricing option by Blackschole

BlackScholes <- function(S, K, r, T, sig, type){
  if(type=="C"){
    d1 <- (log(S/K) + (r + sig^2/2)*T) / (sig*sqrt(T))
    d2 <- d1 - sig*sqrt(T)
    value <- S*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
    return(value)}
  if(type=="P"){
    d1 <- (log(S/K) + (r + sig^2/2)*T) / (sig*sqrt(T))
    d2 <- d1 - sig*sqrt(T)
    value <- (K*exp(-r*T)*pnorm(-d2) - S*pnorm(-d1))
    return(value)}
}

#Analytical solution of Black-Scholes for 3000 different stock price
for(i in 1:3000){
  B[k]<-BlackScholes(S[k],100,0.05,1,0.4,"C")
}
plot(density(B))

#option price by Monte-Carlo simulation

set.seed(2)
n <- 3000
r <- sigma^2/2
K <- 60
T_minus_t <- 1

## the value of the put option by Black-Scholes

d2 <- (log(S/K) + (r - sigma^2/2) * (T_minus_t))/(sigma * sqrt(T_minus_t))
d1 <- (log(S/K) + (r + sigma^2/2) * (T_minus_t))/(sigma * sqrt(T_minus_t))
VP <- K * exp(-r * T_minus_t) * pnorm(-d2) - S * pnorm(-d1)
x <- rlnorm(n) #Note use of default meanlog=0, sdlog=1
y1 <- sapply(x, function(z) max(0, K - exp(z)))
y2 <- sapply(-x, function(z) max(0, K - exp(z)))
y <- (y1 + y2)/2
density(y)
plot(density(y))
mc100 <- t.test(exp(-r * T_minus_t) * y[1:100]) # first 100 simulations
mc1000 <- t.test(exp(-r * T_minus_t) * y[1:1000]) # first 1000 simulations
mcall <- t.test(exp(-r * T_minus_t) * y) # all simulation results
type <- c("Blacksholes", "Monte Carlo for 100", "Monte Carlo for 1000", "Monte Carlo for all")
putestimate <- c(VP, mc100$estimate, mc1000$estimate, mcall$estimate)
putconfintleft <- c(NA, mc100$conf.int[1], mc1000$conf.int[1], mcall$conf.int[1])
putconfintright <- c(NA, mc100$conf.int[2], mc1000$conf.int[2], mcall$conf.int[2])
d <- data.frame(type, putestimate, putconfintleft, putconfintright)
print(d)

x <- rnorm(n)
y <- sapply(x, function(z) max(0, K - z))
anti100 <- t.test(exp(-r * T_minus_t) * y[1:100]) # first 100 simulations
anti1000 <- t.test(exp(-r * T_minus_t) * y[1:1000]) # first 1000 simulations
antiall <- t.test(exp(-r * T_minus_t) * y) # all simulation results
type <- c("Blacksholes", "Antithetic 100", "Antithetic 1000", "Antithetic
all")
putestimate <- c(VP, anti100$estimate, anti1000$estimate, antiall$estimate)
putconfintleft <- c(NA, anti100$conf.int[1], anti1000$conf.int[1],
                    antiall$conf.int[1])
putconfintright <- c(NA, anti100$conf.int[2], anti1000$conf.int[2],
                     antiall$conf.int[2])
d <- data.frame(type, putestimate, putconfintleft, putconfintright)
print(d)

xN <- rnorm(n)
yN <- sapply(xN, function(z) max(0, K - exp(z)))
plot(density(yN))
mcN100 <- t.test(exp(-r * T_minus_t) * yN[1:100]) # first 100 simulations
mcN1000 <- t.test(exp(-r * T_minus_t) * yN[1:1000]) # first 1000 simulations
mcNall <- t.test(exp(-r * T_minus_t) * yN) # all simulation results
type <- c("Blacksholes", "Monte Carlo 100", "Monte Carlo 1000", "Monte
Carlo all")
putestimate <- c(VP, mcN100$estimate, mcN1000$estimate, mcNall$estimate)
putconfintleft <- c(NA, mcN100$conf.int[1], mcN1000$conf.int[1],
                    mcNall$conf.int[1])
putconfintright <- c(NA, mcN100$conf.int[2], mcN1000$conf.int[2],
                     mcNall$conf.int[2])
d <- data.frame(type, putestimate, putconfintleft, putconfintright)
print(d)

#Compare the simulated probability density function of Black-Schole SDE
#with an analytical probability density function

h1<-hist(B,breaks=12,col="orange",xlab="sim_id",ylab="closed_form_price",main="Histogram with probabilty density function for closed price")
xfit<-seq(min(B),max(B),length=40)
yfit<-dnorm(xfit,mean=mean(B),sd=sd(B))
yfit<-yfit*diff(h1$mids[1:2])*length(B)
lines(xfit,yfit,col="black",lwd=2)

h2<-hist(y,breaks=12,col="pink",xlab="sim_id",ylab="Monte-Carlo_simulation_price",main="Histogram with probability density for simulation price")
xfit<-seq(min(y),max(y),length=40)
yfit<-dnorm(xfit,mean=mean(y),sd=sd(y))
yfit<-yfit*diff(h2$mids[1:2])*length(y)
lines(xfit,yfit,col="brown",lwd=2)

#Output-
                    type  putestimate putconfintleft putconfintright
1            Blacksholes  0.28126114             NA              NA
2    Monte Carlo for 100  0.29584862       50.03873        53.34373
3   Monte Carlo for 1000  0.29165719       51.68897        52.71445
4    Monte Carlo for all  0.31183482       52.20828        52.77618

              type   putestimate putconfintleft putconfintright
1       Blacksholes  0.28126114             NA              NA
2    Antithetic 100  0.29584862       57.35991        57.73106
3   Antithetic 1000  0.29165719       57.37837        57.50294
4   Antithetic\nall  0.31183482       57.36371        57.43358

              type   putestimate putconfintleft putconfintright
1        Blacksholes  0.28126114             NA              NA
2    Monte Carlo 100  0.29584862       55.28169        56.02255
3   Monte Carlo 1000  0.29165719       55.69892        55.91098
4   Monte\nCarlo all  0.31183482       55.66868        55.81800




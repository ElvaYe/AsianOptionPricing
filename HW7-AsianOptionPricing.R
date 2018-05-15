AsianOptionPrice <- function(n) {
# parameters:
T <- 1
r <- 0.01
sigma <- 0.2
K <- 100
m <- 12
delta_t <- (T/m)
S0 <- 110
################### stock path function by MC ###################################
stock <- function(S0, r, sigma, delta_t)
{
  S <- matrix(0, nr=n, nc=m+1)
  S[,1] <- S0

  for (i in 1:n)
  {
    for (j in 1:m)
    {
      S[i,j+1] <- S[i,j]*exp((r-1/2*sigma^2)*delta_t+sigma*sqrt(delta_t)*rnorm(1))
    }
  }
  return(S)
}
################### Arithmaric Algorithm #########################################
S <- stock(S0,r,sigma,delta_t)
Y <- matrix(0,n,1)
for (i in 1:n) { Y[i] <- exp(-r*T)*max(mean(S[i,])-K,0) }

################### Geometric Algorithm ##########################################
t <- numeric(m)
for (i in 1:m){t[i]<-i*delta_t}
T_bar <- mean(t)
q <- numeric(m)
for (i in 1:m) {q[i]<-(2*i-1)*t[m+1-i]}
sigma_bar <- sigma/m*sqrt(sum(q)/T_bar)
delta <- 1/2*sigma^2-1/2*sigma_bar^2
d <- (log(S0/K)+(r-delta+1/2*sigma_bar^2)*T_bar)/(sigma_bar*sqrt(T_bar))
# expected value of X by BGM
X_mean <- exp(-delta*T_bar-r*(T-T_bar))*S0*pnorm(d)-exp(-r*T)*K*pnorm(d-sigma_bar*sqrt(T_bar))

# geometric average 
Sg <- matrix(0,n,1)
X <- matrix(0,n,1)
for (i in 1:n)
{
  Sg[i]<-(prod(S[i,])/S0)^(1/m)
  X[i]<-exp(-r*T)*max(Sg[i]-K,0)
}

# Optimal value of b
b <- cov(X,Y)/var(X)
# Estimator Y(b*)
Y_b <- Y-b[1]*(X-X_mean)
################### Result #######################################################
# By first algorithm 
price1 <- mean(Y)
error1 <- sd(Y)/sqrt(n)
# By second algorithm
price2 <- mean(Y_b)
error2 <- sd(Y_b)/sqrt(n)
# Correlation coefficient
corr <- cor(Y,X)[1]
c("Number of replications"=n,"Price by 1st Algo"=price1
  ,"Price by 2nd algorithm"=price2
  ,"Error est. by 1st Algo"=error1,"Error est. by 2nd Algo"=error2
  ,"Correlation"=corr)
}

################### Output #######################################################
# call AsianOptionPrice function with different number of replications
output <- do.call(rbind,lapply(10^(2:6), AsianOptionPrice))
output
write.csv(output,"/Users/yeegoo/Desktop/6233 Stocastic Calculus/HW/HW7_output.csv",row.names = F)


# > output
# Number of replications Price by 1st Algo Price by 2nd algorithm Error est. by 1st Algo Error est. by 2nd Algo Correlation
# [1,]                  1e+02          11.51104               11.70511             1.08093349           0.0361514700   0.9994406
# [2,]                  1e+03          11.27481               11.62910             0.34024901           0.0123245111   0.9993438
# [3,]                  1e+04          11.56960               11.62456             0.10963160           0.0038212786   0.9993924
# [4,]                  1e+05          11.64962               11.62298             0.03441584           0.0011978680   0.9993941
# [5,]                  1e+06          11.61550               11.62316             0.01089438           0.0003793298   0.9993936

# Conclusion:
# The price calculated by both algorithm are similar but the second algorithm (variance reduction) has much smaller error estimator than the first one. 
#  And when the number of replication  increases, both error estimators decrease.
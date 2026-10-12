#
#
#   MATH 2650: Computational Statistics
#   Lab 06
#   16 October 2026
# 
#   Second Lab Group 5:
#
#

install.packages("ISwR")
library(ISwR)

#### (1.) Jackknife Bernoulli Empirical Bayes ####
x <- malaria$mal
n <- nrow(malaria)

##### (1b.) Jackknife Estimate #####
a_i    <- 1/2
b_i    <- 1/2

tx_hat <- (sum(x) + a_i)/(n + b_i) # Empirical Bayes estimator

# We obtain an analytical Jackknife estimator of 0.27002, so we estimate that the proportion 
# of patients in the population who present with symptoms consistent with malaria is roughly 27%.

# Analytical Jackknifed Estimate
Jtx_a <- n*tx_hat - ((n-1) / (n*(n-1+b_i)))*((n-1)*sum(x) + n*a_i)
print(paste("Analytical Jackknifed estimate:", round(Jtx_a, 6)))

##### (1c.) Jackknifed Estimator, Variance, Bias #####
thetas <- vector('numeric', length = n)

for(j in 1:n){
  thetas[j] <- (sum(x[-j]) + a_i) / (n - 1 + b_i) # Use n - 1 instead of n
}

Jtx <- n*tx_hat - (n-1)*mean(thetas) # Numerical Jackknifed t(x), J[t(x)]
print(paste("Numerical Jackknifed estimate of t(x):", round(Jtx, 6)))

t_stars <- n*tx_hat - (n-1)*thetas # Pseudovalues
mean(t_stars)

# Use Var[J{t(x)}] as an estimator of variance of t(x)
Var_Jtx <- (1/(n*(n-1)))*sum((t_stars - Jtx)^2)
print(paste("Jackknife estimate of variance of t(x):", round(Var_Jtx, 6)))

# Jackknife estimate of bias in t(x)
tbar_x <- mean(thetas)
Btx <- (n-1)*(tbar_x - tx_hat)
print(paste("Jackknife estimate of bias of t(x):", round(Btx, 6)))

c(Jtx - 2*sqrt(Var_Jtx), Jtx + 2*sqrt(Var_Jtx)) # 75% interval

#### (2.) Exponential Jackknife Estimator ####
y <- malaria$ab
ty_hat <- 1/mean(y)

yj <- vector('numeric', length = n)

for(j in 1:n){
  yj[j] <- 1/mean(y[-j]) 
}

Jty <- n*ty_hat - (n-1)*mean(yj) # Numerical Jackknifed t(y), J[t(y)]
print(paste("Numerical Jackknifed estimate of t(y):", round(Jty, 6)))
print(paste("Original estimator t(y):", round(ty_hat, 6)))

ty_stars <- n*ty_hat - (n-1)*yj # Pseudovalues

# Use Var[J{t(y)}] as an estimator of variance of t(y)
Var_Jty <- sum((ty_stars - Jty)^2)/(n*(n-1))
print(paste("Jackknife estimate of variance of t(y):", Var_Jty))

c(Jty - 4*sqrt(Var_Jty), Jty + 4*sqrt(Var_Jty)) # 95% interval

# Jackknife estimate of bias in t(y)
tbar_y <- mean(yj)
Bty <- (n-1)*(tbar_y - ty_hat)
print(paste("Jackknife estimator of bias of t(y):", Bty))
print(paste("Estimated bias of t(y):", ty_hat - Jty))


#### (3.) Bias of Gaussian Variance Estimator ####
z <- cystfibr$fev1
m <- length(z)

sigmasq_hat <- mean((z-mean(z))^2)
print(paste("Estimated variance of FEV1:", sigmasq_hat))

sigmasqs <- vector('numeric', length = m)

for(i in 1:m){
  sigmasqs[i] <- mean((z[-i] - mean(z[-i]))^2)
}

Jtz <- m*sigmasq_hat - (m-1)*mean(sigmasqs)
print(paste("Jackknife estimator of variance of FEV1:", Jtz))

Btz <- sigmasq_hat - Jtz
print(Btz)

print(paste("Sample variance of FEV1:", var(z)))

#
#
#   MATH 2650: Computational Statistics
#   Lab 08
#   05 November 2026
# 
#   Second Lab Group 5:
#   Destiny Agubata, Vonn Russell, Angie Wang
#
#

# install.packages("ISwR")
# library(ISwR)

#### (1.) Parametric & Nonparametric Bootstrap Interval ####
B       <- 10000
malaria <- ISwR::malaria
n       <- nrow(malaria)
theta_0 <- sum(malaria$mal == 1)/n # Initial estimate 

# Nonparametric bootstrap
thetas_nonp <- thetas_p <- vector('numeric', length = B)

set.seed(1924)
for(i in 1:B){
  mal_sample     <- sample(malaria$mal, n, replace = TRUE)
  thetas_nonp[i] <- sum(mal_sample == 1)/n
}

hist(thetas_nonp, freq = FALSE)
abline(v = theta_0, col = 'darkblue', lwd = 2)

# 99% bootstrap interval
alpha <- 0.01
print("99% Nonparametric Bootstrap interval for mean")
quantile(thetas_nonp, probs = c(alpha/2, 1-alpha/2))

# Parametric bootstrap
theta_0 <- sum(malaria$mal == 1)/n # Success probability from original distribution

set.seed(1924)
for(i in 1:B){
  mal_sample  <- rbinom(n = 1, size = n, prob = theta_0)
  thetas_p[i] <- mal_sample/n
}

hist(thetas_p, freq = FALSE)
abline(v = theta_0, col = 'darkblue', lwd = 2)

# 99% bootstrap interval
alpha <- 0.01
print("99% Parametric Bootstrap interval for theta")
quantile(thetas_p, probs = c(alpha/2, 1-alpha/2))

#### (2.) Nonparametric Intervals for Mean & Median ####
B        <- 10000
antibody <- malaria$ab

##### (2a.) 90% Interval for Median #####
med_0    <- median(antibody) # "Population" median
med_vec  <- vector('numeric', length = B)

set.seed(2005)
for(i in 1:B){
  ab_med     <- median(sample(antibody, n, replace = TRUE))
  med_vec[i] <- ab_med
}

hist(med_vec, freq = FALSE)
abline(v = med_0, col = 'darkblue', lwd = 2)

# 90% bootstrap interval
alpha <- 0.1
print("90% Bootstrap interval for median")
quantile(med_vec, probs = c(alpha/2, 1-alpha/2))

##### (2b.) 90% Interval for Mean #####
mean_0   <- mean(antibody) # "Population" mean
mean_vec <- vector('numeric', length = B)

set.seed(2005)
for(i in 1:B){
  ab_mean     <- mean(sample(antibody, n, replace = TRUE))
  mean_vec[i] <- ab_mean
}

hist(mean_vec, freq = FALSE)
abline(v = mean_0, col = 'darkblue', lwd = 2)

# 90% bootstrap interval
alpha <- 0.1
print("90% Bootstrap interval for mean")
quantile(mean_vec, probs = c(alpha/2, 1-alpha/2))

#### (3.) Mean Absolute Deviation ####
mad_0   <- mad(antibody) # MAD for antibody levels
mad_0
mad_vec <- vector('numeric', length = B)

set.seed(1998)
for(i in 1:B){
  ab_mad     <- mad(sample(antibody, n, replace = TRUE))
  mad_vec[i] <- ab_mad
}

hist(mad_vec, freq = FALSE)
abline(v = mad_0, col = 'darkblue', lwd = 2)

# 95% bootstrap interval
alpha <- 0.05
print("95% Bootstrap interval of MAD")
quantile(mad_vec, probs = c(alpha/2, 1-alpha/2))


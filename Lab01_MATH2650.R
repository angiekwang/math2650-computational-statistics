###############################################################################
#
#
#   MATH 2650: Computational Statistics
#   Lab 01
#   10 September, 2026
# 
#   First Lab Group 1
#   Angie Wang
#
#
###############################################################################

# Question 1: Use uniroot() function to calculate MLE for Bernoulli and Exponential
#             distributions. Compare results to the traditional analytical approach.

#### (1a.) uniroot Bernoulli ####
mag <- 1*(attenu$mag >6) # Obtain all magnitudes greater than 6

##### Analytical MLE #####
analytic_theta <- mean(mag)
analytic_theta

##### Numerical MLE #####
calc_deriv_bern <- function(theta, y){
  # Derivative of log-likelihood function of the Bernoulli pdf
  y <- y[complete.cases(y)]
  n <- length(y)
  
  sum(y)/theta - (n - sum(y))*(1/(1-theta))
}

bern_root <- uniroot(calc_deriv_bern, y = mag, lower = 0, upper = max(mag))
numeric_theta <- bern_root$root
numeric_theta


#### (1b.) uniroot Exponential ####
accel <- attenu$accel

##### Analytical MLE #####
analytic_lambda <- length(accel)/sum(accel)
analytic_lambda

##### Numerical MLE #####
calc_deriv_exp <- function(lambda, x){
  # Derivative of log-likelihood function of the Exponential pdf
  
  x <- x[complete.cases(x)]
  n <- length(x)
 
  (n/lambda) - sum(x)
}

exp_root <- uniroot(calc_deriv_exp, x = accel, lower = 0, upper = 10)
numeric_lambda <- exp_root$root
numeric_lambda

###############################################################################

# Question 2: Use the Probability Integral Transform method to generate random 
#             numbers from Gumbel and Lomax distributions.

#### (2a.) PIT Gumbel ####

# Constants:
N <- 20000 # Number of random draws
BETA_1 <- 1 # Testing with beta = 1 and beta = 3
BETA_3 <- 3

z <- runif(N) # Obtain N draws from uniform distribution

calc_gumbel_inv <- function(z, beta){
  # Inverse gumbel cdf
  -beta*log(-log(z))
}

# Test with parameters beta = 1 and beta = 3
gumbel_dist_beta1 <- calc_gumbel_inv(z = z, beta = BETA_1)
gumbel_dist_beta3 <- calc_gumbel_inv(z = z, beta = BETA_3)

# Generate histograms
hist(gumbel_dist_beta1, freq = FALSE, main = paste(N, ' samples', sep = ''))
lines(density(gumbel_dist_beta1), lwd = 2, col = 'darkblue')

hist(gumbel_dist_beta3, freq = FALSE, main = paste(N, ' samples', sep = ''))
lines(density(gumbel_dist_beta3), lwd = 2, col = 'darkblue')

#### (2b.) PIT Lomax ####

# Constants:
M <- 20000 # Number of random draws
ALPHA <- 20 # Test with 20
LAMBDA_5 <- 5 # Test with lambda = 5 and lambda = 10
LAMBDA_10 <- 10 


w <- runif(M) # Obtain M draws from uniform distribution

calc_lomax_inv <- function(w, lambda, alpha){
  # Inverse lomax cdf
  lambda*((1-w)^(-1/alpha) - 1) 
}

# Test with parameters lambda = 5 and lambda = 10
lomax_dist_lambda5 <- calc_lomax_inv(w = w, lambda = LAMBDA_5, alpha = ALPHA) 
lomax_dist_lambda10 <- calc_lomax_inv(w = w, lambda = LAMBDA_10, alpha = ALPHA)

# Generate histograms
hist(lomax_dist_lambda5, freq = FALSE, main = paste(M, ' samples', sep = ''))
lines(density(lomax_dist_lambda5), lwd = 2, col = 'darkblue')

hist(lomax_dist_lambda10, freq = FALSE, main = paste(M, ' samples', sep = ''))
lines(density(lomax_dist_lambda10), lwd = 2, col = 'darkblue')

###############################################################################

# Question 3: Generate a random number using the Middle-Square / Von Neumann
#             method using a predetermined starting seed.

#### (3.) Middle-Square / Von Neumann ####

seed <- 812016 # Starting seed

seed_str <- as.character(seed)
seed_len <- nchar(seed_str)
seed_sq <- as.character(seed^2)
seed_sq_len <- nchar(seed_sq)

start <- as.integer(seed_sq_len/2 - seed_len/2 + 1) # Add 1 since R starts counting at 1
end <- as.integer(start - 1 + seed_len)

mid_square <- as.numeric(substring(seed_sq, start, end))
mid_square

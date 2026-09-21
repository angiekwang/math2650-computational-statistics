#
#
#   MATH 2650: Computational Statistics
#   Lab 03
#   24 September 2026
# 
#   First Lab Group 1
#
#
#

#### (1.) Finding Quantiles ####

first_quart <- 0.25
median      <- 0.5
third_quart <- 0.75

##### (1a.) Bisection Method for Geometric Distribution #####

P           <- 1/20 # Set parameter value
TOL_1a      <- 1    # Set tolerance
MAX_IT      <- 100  # Set maximum number of iterations
x_l         <- 1    # Initial lower bound
x_u         <- 50   # Initial upper bound
change_in_x <- x_u - x_l
counter     <- 0    # Counts number of steps
x_grid <- seq(1,50) # For plotting

g <- function(a, q){ # Geometric cdf minus quartile
  x <- floor(a)
  1 - (1 - P)^floor(a) - q
}

estimate_geom_quantile <- function(q){
  
  while(abs(change_in_x) > TOL_1a){
    
    y <- (x_l + x_u)/2
      
    if(sign(g(y,q)) == sign(g(x_l,q))){
      x_l <- y
    } else{
      x_u <- y
    }
    
    counter     <- counter + 1
    change_in_x <- x_u - x_l
    
    if(counter >= MAX_IT){
      warning('Max number of iterations reached.')
      break
    }
  }
  cat('Reached estimate of', floor(y), 'in', counter, 'steps.')
  plot(x_grid, g(x_grid,q))
}

estimate_geom_quantile(first_quart) 
estimate_geom_quantile(median)      
estimate_geom_quantile(third_quart) 


##### (1b.) Secant Method for Type II Generalized Logistic Distribution #####

GAMMA   <- 2
TOL_1b  <- 1e-05 # Set tolerance
MAX_IT  <- 100   # Set maximum number of iterations
x <- vector('numeric', length = MAX_IT)
x[1]     <- 0     # First initial point
x[2]     <- 1     # Second initial point
dx      <- x[2] - x[1]
i <- 2     # Counts number of steps

h <- function(a, q){ # cdf minus quartile
  1 - ( exp(-a) / (1 + exp(-a)) )^GAMMA - q
}

estimate_logistic_quantile <- function(q){
  
  while(abs(dx) > TOL_1b){
    
    num <- h(x[i], q) * (x[i] - x[i-1])
    den <- h(x[i], q) - h(x[i-1], q)
    
    x[i+1] <- x[i] - num/den
    
    dx <- x[i+1] - x[i] # Update checking condition
    i  <- i + 1         # Update counter

    if(i >= MAX_IT+2){
      warning('Max number of iterations reached.')
      break
    }
  }
  cat('Reached estimate of', x[i], 'in', i-2, 'steps.')
  plot(1:i, x[1:i], type = 'b', xlab = 'steps')
}


estimate_logistic_quantile(first_quart) 
estimate_logistic_quantile(median)      
estimate_logistic_quantile(third_quart) 

#### (2.) Fisher Scoring for MLE ####

##### (2a.) Poisson Distribution #####
z <- as.numeric(InsectSprays$count)
TOL_2a <- 1e-05
MAX_ITER <- 100
LAMBDA_0 <- 1 # Starting value
lambda <- vector('numeric', length = MAX_ITER)
lambda[1] <- LAMBDA_0
change_in_lambda <- 1 # Initialize to be greater than tolerance
a <- 1 # Initialize counter at 1

pois_derivloglik <- function(lambda){ # First derivative of log-likelihood
  (1/lambda)*sum(z) - length(z)
}

while(change_in_lambda > TOL_2a){
  
  score2a <- (1/lambda[a])*sum(z) - length(z)
  info2a <- (length(z))*(mean(z))/lambda[a]^2
  
  update <- score2a/info2a # Fisher update
  
  lambda[a+1] <- lambda[a] + update
  
  change_in_lambda <- abs(lambda[a+1] - lambda[a]) # Update check against tolerance
  a <- a + 1 # Update counter
  
  if(a >= MAX_ITER){
    warning('Max number of iterations reached.')
    break
  }
}

pois_uni <- uniroot(pois_derivloglik, lower = 0, upper = max(z))

lambda[a] # Fisher scoring estimate
mean(z) # Analytical MLE estimate
pois_uni$root # Uniroot estimate


##### (2b.) Exponential Distribution #####
w <- as.numeric(sunspot.year) 
TOL_2b <- 1e-05 # Set tolerance
MAX_ITER <- 100
THETA_0 <- 0.01 # Set initial starting value
theta <- vector('numeric', length = MAX_ITER)
theta[1] <- THETA_0 
change_in_theta <- 1 # Initialize to be greater than tolerance
b <- 1 # Initialize counter


while(change_in_theta > TOL_2b){
  
  score2b <- length(w)/theta[b] - sum(w)
  info2b <- length(w)/theta[b]^2
  
  update <- score2b/info2b # Fisher update
  
  theta[b+1] <- theta[b] + update
  
  change_in_theta <- abs(theta[b+1] - theta[b]) # Update check against tolerance
  b <- b + 1 # Update counter
  
  if(b >= MAX_ITER){
    warning('Max number of iterations reached.')
    break
  }
}

hist(w, freq = FALSE)
curve(dexp(x, theta[b]), add = TRUE, lwd = 2,)


##### (2c.) Rayleigh Distribution #####
s <- airquality$Wind
BETA_0 <- 0.01
TOL_2c <- 1e-05 # Set tolerance
MAX_ITER <- 100
change_in_beta <- 1 # Initialize to be greater than tolerance
beta <- vector('numeric', length = MAX_ITER)
beta[1] <- BETA_0
c <- 1 # Initialize counter


while(change_in_beta > TOL_2c){
  
  score2c <- (length(s)/beta[c]) - sum(s^2)
  info2c <- length(s)/beta[c]^2
  
  update <- score2c/info2c # Fisher update
  
  beta[c+1] <- beta[c] + update
  
  change_in_beta <- abs(beta[c+1] - beta[c]) # Update check against tolerance
  c <- c + 1 # Update counter
  
  if(c >= MAX_ITER){
    warning('Max number of iterations reached.')
    break
  }
}

hist(s, freq = FALSE)
curve(beta[c] * x * exp(-beta[c] * x^2), add = TRUE, lwd = 2)


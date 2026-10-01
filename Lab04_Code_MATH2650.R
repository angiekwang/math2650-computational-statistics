#
#
#   MATH 2650: Computational Statistics
#   Lab 04
#   01 October 2026
# 
#   First Lab Group 1:
#   Matt O'Connor, Vonn Russell, Angie Wang
#
#

#### (1.) Rejection Sampling ####
SEED1 <- 1979
B1    <- 10000
B2    <- 50000

##### (1a.) Kumaraswamy Distribution #####
THETA1 <- 2
GAMMA1 <- 5

dkumaraswamy <- function(x, theta, gamma){
    ifelse(x > 0 & x < 1, 
           (theta*gamma)*x^(theta-1)*(1-x^theta)^(gamma-1), 0)
}

g <- function(x){ # Derivative of reciprocal of importance ratio
  (1/x) - (8*x)/(1-x^2)} 

guni <- uniroot(g, lower = 0.01, upper = 0.99) 
M1a  <- dkumaraswamy(guni$root, THETA1, GAMMA1)/dunif(guni$root) 

###### (1a.i) Obtain B1 = 10000 samples ######
xb_i <- vector('numeric', length = B1)
U_i  <- vector('numeric', length = B1)
fX_i <- vector('numeric', length = B1)
k    <- 1 

total_draws <- 0 # For calculating efficiency

set.seed(SEED1)

while(k <= B1){
  
  U_i[k]  <- runif(1)
  xb_i[k] <- runif(1)
  
  if(U_i[k] < dkumaraswamy(xb_i[k], THETA1, GAMMA1)/(M1a*dunif(xb_i[k]))){
    
    fX_i[k]   <- xb_i[k]
    k         <- k + 1
  }
  total_draws <- total_draws + 1
}

hist(fX_i[1:B1], main = '10000 samples from Kumaraswamy', xlab = 'x', freq = FALSE)
curve((THETA1*GAMMA1)*x^(THETA1-1)*(1-x^THETA1)^(GAMMA1-1), add = TRUE, lwd = 2)

efficiency <- round(B1/total_draws, digits = 4)
efficiency_percent <- efficiency*100
print(paste("The efficiency is", efficiency, "or", efficiency_percent, "percent for",
            B1, "samples."))

###### (1a.ii) Obtain B2 = 50000 samples ######
xb_ii <- vector('numeric', length = B2)
U_ii  <- vector('numeric', length = B2)
fX_ii <- vector('numeric', length = B2)
j     <- 1 # Counter

total_draws <- 0 # Reset 

set.seed(SEED1)

while(j <= B2){
  
  U_ii[j]  <- runif(1)
  xb_ii[j] <- runif(1)
  
  if(U_ii[j] < dkumaraswamy(xb_ii[j], THETA1, GAMMA1)/(M1a*dunif(xb_ii[j]))){
    
    fX_ii[j]  <- xb_ii[j]
    j         <- j + 1
  }
  total_draws <- total_draws + 1
}

hist(fX_ii[1:B2], main = '50000 samples from Kumaraswamy', xlab = 'x', freq = FALSE)
curve((THETA1*GAMMA1)*x^(THETA1-1)*(1-x^THETA1)^(GAMMA1-1), add = TRUE, lwd = 2)

efficiency <- round(B2/total_draws, digits = 4)
efficiency_percent <- efficiency*100

print(paste("The efficiency is", efficiency, "or", efficiency_percent, "percent for",
            B2, "samples."))

##### (1b.) U-quadratic Distribution #####
ALPHA1 <- 0
BETA1  <- 3

duquadratic <- function(y, alpha, beta){
  
  ifelse(y >= alpha & y <= beta,
         12/(beta - alpha)^3 * (y - (beta + alpha)/2)^2, 
         0)
}

# Find M1b graphically
curve(12/(BETA1 - ALPHA1)^3 * (x - (BETA1 + ALPHA1)/2)^2, lwd = 2, col = "blue",
      xlim = c(ALPHA1, BETA1))
curve(3.01*dunif(x, min = ALPHA1, max = BETA1), lwd = 2, add = TRUE, col = "red")

# Zoom in
curve(12/(BETA1 - ALPHA1)^3 * (x - (BETA1 + ALPHA1)/2)^2, lwd = 2, col = "blue",
      xlim = c(-0.01, 0.01), ylim = c(0.99, 1.01))
curve(3.01*dunif(x, min = ALPHA1, max = BETA1), lwd = 2, add = TRUE, col = "red")

M1b <- 3.01

###### (1b.i) Obtain B1 = 10000 samples ######
yb_i <- vector('numeric', length = B1)
U_i  <- vector('numeric', length = B1)
fY_i <- vector('numeric', length = B1)
k    <- 1 # Counter

total_draws <- 0 # Reset 

set.seed(SEED1)

while(k <= B1){
  
  U_i[k]  <- runif(1)
  yb_i[k] <- runif(1, min = 0, max = 3)
  
  if(U_i[k] < duquadratic(yb_i[k], ALPHA1, BETA1)/(M1b*dunif(yb_i[k], min = 0, max = 3))){
    
    fY_i[k]   <- yb_i[k]
    k         <- k + 1
  }
  total_draws <- total_draws + 1
}

hist(fY_i[1:B1], main = '10000 samples from U-quadratic', xlab = 'y', freq = FALSE)
curve(12/(BETA1 - ALPHA1)^3 * (x - (BETA1 + ALPHA1)/2)^2, add = TRUE, lwd = 2)

efficiency <- round(B1/total_draws, digits = 4)
efficiency_percent <- efficiency*100
print(paste("The efficiency is", efficiency, "or", efficiency_percent, "percent for",
            B1, "samples."))


###### (1b.ii) Obtain B2 = 50000 samples ######
yb_ii <- vector('numeric', length = B2)
U_ii  <- vector('numeric', length = B2)
fY_ii <- vector('numeric', length = B2)
j     <- 1 # Counter

total_draws <- 0 # Reset 

set.seed(SEED1)

while(j <= B2){
  
  U_ii[j]  <- runif(1)
  yb_ii[j] <- runif(1, min = 0, max = 3)
  
  if(U_ii[j] < duquadratic(yb_ii[j], ALPHA1, BETA1)/(M1b*dunif(yb_ii[j], min = 0, max = 3))){
    
    fY_ii[j]  <- yb_ii[j]
    j         <- j + 1
  }
  total_draws <- total_draws + 1
}

hist(fY_ii[1:B2], main = '50000 samples from U-quadratic', xlab = 'y', freq = FALSE)
curve(12/(BETA1 - ALPHA1)^3 * (x - (BETA1 + ALPHA1)/2)^2, add = TRUE, lwd = 2)

efficiency <- round(B2/total_draws, digits = 4)
efficiency_percent <- efficiency*100
print(paste("The efficiency is", efficiency, "or", efficiency_percent, "percent for",
            B2, "samples."))



#### (2.) Importance Sampling/Resampling ####
SEED2 <- 1789
B1 <- 10000
B2 <- 50000 

##### (2a.) Kumaraswamy Distribution #####
THETA2 <- 2
GAMMA2 <- 2

###### (2a.i) Uniform Proposal Density ######
xb_2ai <- vector('numeric', length = B1)
w_2ai <- vector('numeric', length = B1)

set.seed(SEED2)

for(n in 1:B1){
  xb_2ai[n] <- runif(1)
  w_2ai[n] <- dkumaraswamy(xb_2ai[n], THETA2, GAMMA2)/dunif(xb_2ai[n])
}

# Calculate mean
first_mom <- sum(w_2ai*xb_2ai)/sum(w_2ai)
print(paste("Estimated Kumaraswamy mean using std uniform proposal density is", round(first_mom, 4)))

# Calculate variance
second_mom <- sum(w_2ai*xb_2ai^2)/sum(w_2ai)
variance <- second_mom - first_mom^2
print(paste("Estimated Kumaraswamy variance using std uniform proposal density is", round(variance, 4)))


###### (2a.ii) Kumaraswamy w/ Beta Proposal Density ######
xb_2aii <- vector('numeric', length = B1)
w_2aii <- vector('numeric', length = B1)

set.seed(SEED2)

for(n in 1:B1){
  xb_2aii[n] <- rbeta(1, 2, 2)
  w_2aii[n] <- dkumaraswamy(xb_2aii[n], THETA2, GAMMA2)/dbeta(xb_2aii[n], 2, 2)
}

# Calculate mean
first_mom <- sum(w_2aii*xb_2aii)/sum(w_2aii)
print(paste("Estimated Kumaraswamy mean using beta proposal density is", round(first_mom, 4)))

# Calculate variance
second_mom <- sum(w_2aii*xb_2aii^2)/sum(w_2aii)
variance <- second_mom - first_mom^2
print(paste("Estimated Kumaraswamy variance using beta proposal density is", round(variance, 4)))


##### (2b.) U-quadratic Distribution #####
ALPHA2 <- -2
BETA2  <- 2

###### (2b.i) U-quadratic w/ Uniform Proposal Density ######
xb_2bi <- vector('numeric', length = B2)
w_2bi <- vector('numeric', length = B2)

set.seed(SEED2)

for(n in 1:B2){
  xb_2bi[n] <- runif(1, -2, 2)
  w_2bi[n] <- duquadratic(xb_2bi[n], ALPHA2, BETA2)/dunif(xb_2bi[n], -2, 2)
}

# Calculate mean and variance
first_mom <- sum(w_2bi[1:B1]*xb_2bi[1:B1])/sum(w_2bi[1:B1])
second_mom <- sum(w_2bi[1:B1]*xb_2bi[1:B1]^2)/sum(w_2bi[1:B1])
variance <- second_mom - first_mom^2
print(paste("Estimated U-quadratic mean using uniform proposal density is", round(first_mom, 4)))
print(paste("Estimated U-quadratic variance using uniform proposal density is", round(variance, 4)))


###### (2b.ii) U-quadratic w/ Normal Proposal Density ######
xb_2bii <- vector('numeric', length = B2)
w_2bii <- vector('numeric', length = B2)

set.seed(SEED2)

for(n in 1:B2){
  xb_2bii[n] <- rnorm(1)
  w_2bii[n] <- duquadratic(xb_2bii[n], ALPHA2, BETA2)/dnorm(xb_2bii[n])
}

# Calculate mean, excluding importance ratios equal to zero 
keep    <- w_2bii != 0
w_2bii  <- w_2bii[keep]
xb_2bii <- xb_2bii[keep]

first_mom <- sum(w_2bii[1:B1]*xb_2bii[1:B1])/sum(w_2bii[1:B1])
print(paste("Estimated U-quadratic mean using std normal proposal density is", round(first_mom, 4)))

# Calculate variance
second_mom <- sum(w_2bii[1:B1]*xb_2bii[1:B1]^2)/sum(w_2bii[1:B1])
variance <- second_mom - first_mom^2
print(paste("Estimated U-quadratic variance using std normal proposal density is", round(variance, 4)))


##### (2c.) U-quadratic w/ Importance Resampling #####
B1 <- 10000 # Number of samples to keep

###### (2c.i) U-quadratic Resampling w/ Uniform Proposal Density ######
xb_2ci <- vector('numeric', length = B1)
w_2ci <- vector('numeric', length = B1)

set.seed(SEED2)

xb_2ci <- sample(xb_2bi, B1, replace = FALSE, prob = w_2bi)

# Calculate mean and variance
first_mom <- mean(xb_2ci)
second_mom <- mean(xb_2ci^2)
variance <- second_mom - first_mom^2
print(paste("Estimated resampling U-quadratic mean using uniform proposal density is", round(first_mom, 4)))
print(paste("Estimated resampling U-quadratic variance using uniform proposal density is", round(variance, 4)))


###### (2c.ii) U-quadratic Resampling w/ Normal Proposal Density ######
xb_2cii <- vector('numeric', length = B1)
w_2cii <- vector('numeric', length = B1)

set.seed(SEED2)

xb_2cii <- sample(xb_2bii, B1, replace = FALSE, prob = w_2bii)

# Calculate mean and variance
first_mom <- mean(xb_2cii)
second_mom <- mean(xb_2cii^2)
variance <- second_mom - first_mom^2
print(paste("Estimated resampling U-quadratic mean using standard normal proposal density is", round(first_mom, 4)))
print(paste("Estimated resampling U-quadratic variance using standard normal proposal density is", round(variance, 4)))





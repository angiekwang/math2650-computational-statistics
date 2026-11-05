#
#
#   MATH 2650: Computational Statistics
#   Lab 09
#   12 November 2026
# 
#   Third Lab Group 2
#
#

##### (1b.) Beta posterior sampling function #####

betapost <- function(M, x, alpha, beta){
  a <- alpha + sum(x)
  b <- beta + length(x) - sum(x)
  
  rbeta(M, a, b)
}


##### (1c.) Sample from data #####
# install.packages("ISwR")
library(ISwR)

mal <- ISwR::malaria$mal

set.seed(857)
posterior1 <- betapost(M = 20000, x = mal, alpha = 0.5, beta = 0.5)

# 99% credible interval
alpha <- 0.01
quantile(posterior1, c(alpha/2, 1-alpha/2))

mean(posterior1)
sd(posterior1)

##### (2a.) Gamma posterior sampling function #####

gammapost <- function(M, z, alpha, beta){
  a <- alpha + length(z)
  b <- beta + sum(z)
  
  rgamma(M, a, b)
  
}

##### (2b.) Sample from data #####
ab <- ISwR::malaria$ab

set.seed(1557)
posterior2 <- gammapost(M = 20000, z = ab, alpha = 0.01, beta = 0.01)

# 90% credible interval
alpha <- 0.1
quantile(posterior2, c(alpha/2, 1-alpha/2))

mean(posterior2)
sd(posterior2)


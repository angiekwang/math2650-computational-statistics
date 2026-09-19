###############################################################################
#
#
#   MATH 2650: Computational Statistics
#   Lab 02
#   17 September 2026
# 
#   First Lab Group 1
#   Angie Wang
#
#
###############################################################################

# Question 1: Generate random numbers from Rademacher distribution by 
#             transforming Bernoulli R.V.

#### (1a.) Rademacher Random Numbers ####
B <- 20000 # Number of trials
P <- 1/2   # Probability of success 

gen_rademacher_vars <- function(n){
  # Generate Rademacher R.V.s by transforming Bernoulli R.V.s 
  x <- rbinom(n = n, size = 1, prob = P) # Generate Bernoulli R.V.s
  r <- 2*x - 1                           # Transform to Rademacher R.V.
}

rademacher_draws <- gen_rademacher_vars(B)
mean(rademacher_draws) 

# Interpretation: The sample mean is close to zero because the Rademacher R.V.
#                 takes on values -1 and 1, each with 1/2 probability. Therefore
#                 the theoretical expected value E(x) = (-1)(1/2) + (1)(1/2) = 0.


#### (1b.) Normal & Lognormal Random Numbers ####
N     <- 20000 # Number of trials
MU    <- 5     # Mean
SIGMA <- 1/2   # Standard deviation

##### (1b.i) Lognormal Random Numbers #####
gen_lognorm_vars <- function(mu, sigma, n){
  # Generate log normal R.V.s from transforming normal R.V.s
    z <- rnorm(n, mu, sigma) # Obtain normal R.V.s
    x <- exp(z)              # Exponentiate normal R.V.s
}

lognorm_draws <- gen_lognorm_vars(MU, SIGMA, N)
dens_lognorm  <- density(lognorm_draws)

hist(lognorm_draws, prob = TRUE)
lines(dens_lognorm, col = "red", lwd = 2)


##### (1b.ii) Normal R.V. from Lognormal #####

norm_draws <- log(lognorm_draws) # Undo log transformation
dens_norm  <- density(norm_draws)

hist(norm_draws, prob = TRUE)
lines(dens_norm, col = "blue", lwd = 2)

# Interpretation: This second plot does appear normally distributed and centered
#                 at 5. This makes sense, as undoing the log transformed variables
#                 results in our original random variables, which come from a 
#                 normal distribution with mean of 5.  


# Question 2: Explore Markov Chains for t distributions with various degrees of freedom.

#### (2a.) Markov Chain ####
DF_1  <- 2
DF_2  <- 5
STEPS <- 101 # 100 steps + 1 since we start at 0

##### (2a.i) First Random Walk with 2 df #####
rand_walk_1    <- vector('numeric', length = STEPS)
rand_walk_1[1] <- 0 # Initialize starting point at 0

for(i in 2:STEPS){
  # Generate random walk from t-distribution
  rand_walk_1[i] <- rand_walk_1[i-1] + rt(1, DF_1)
}

par(mfrow = c(1, 2))
plot(1:STEPS, rand_walk_1, type = 'l', xlab = 'Steps', ylab = 'Random Walk', main = 'df = 2')

##### (2a.ii) Second Random Walk with 5 df #####
rand_walk_2    <- vector('numeric', length = STEPS)
rand_walk_2[1] <- 0 # Initialize

for(i in 2:STEPS){
  # Generate random walk from t-distribution
  rand_walk_2[i] <- rand_walk_2[i-1] + rt(1, DF_2)
}

plot(1:STEPS, rand_walk_2, type = 'l', xlab = 'Steps', ylab = 'Random Walk', main = 'df = 5')

# Interpretation: Our first random walk with 2 degrees of freedom, in general,
#                 has a larger deviation from the initial value of 0 than the 
#                 random walk with 5 degrees of freedom. This makes sense because 
#                 the random walk with 2 df has more variability (step size comes 
#                 from distribution with fatter tails) and is more likely to 
#                 produce large steps. The actual behavior of the random walk
#                 itself is difficult to generalize.


#### (2b.) Markov Chain with Stopping Rule ####
EPSILON <- 5 # Set tolerance (range of random walk)

##### (2b.i) First Random Walk with 2 df #####
rand_walk_1a_chain <- 0 # "Vector" of random walks
rand_walk_1a_check <- 0 # To check against tolerance
i                  <- 2 # Counter

while(abs(rand_walk_1a_check) < EPSILON){
  rand_walk_1a_chain[i] <- rand_walk_1a_chain[i-1] + rt(1, DF_1)
  rand_walk_1a_check    <- rand_walk_1a_chain[i]
  i                     <- i + 1
}

plot(rand_walk_1a_chain, type = 'l', xlab = paste('Steps = ', i-2, sep = ''), ylab = 'Random Walk',
     main = paste('df = ', DF_1, '; Tolerance = ', EPSILON, sep = ''))

##### (2b.ii) Second Random Walk with 5 df #####
rand_walk_2a_chain <- 0
rand_walk_2a_check <- 0
j                  <- 2 # Counter

while(abs(rand_walk_2a_check) < EPSILON){
  rand_walk_2a_chain[j] <- rand_walk_2a_chain[j-1] + rt(1, DF_2)
  rand_walk_2a_check    <- rand_walk_2a_chain[j]
  j                     <- j + 1
}

plot(rand_walk_2a_chain, type = 'l', xlab = paste('Steps = ', j-2, sep = ''), ylab = 'Random Walk',
     main = paste('df = ', DF_2, '; Tolerance = ', EPSILON, sep = ''))

# Interpretation: The random walk with df = 5 usually has more steps than the 
#                 df = 2 random walk before it exceeds the tolerance threshold. 
#                 Again, this is likely because the steps in the df = 2 walk 
#                 are larger and therefore more likely to jump past the threshold
#                 with a fewer number of steps.

# Question 3: Generate pseudorandom numbers using the Von Neumann/Middle Square Method.
#             Test with different seeds.

#### (3.) Von Neumann / Middle Square Method ####
SEED_1 <- 812016
SEED_2 <- 4000
M      <- 6 # Generate 6 random numbers

gen_midsq_num <- function(seed_num, n){
  
  seed_str <- as.character(seed_num)
  seed_len <- nchar(seed_str)
  
  if(seed_len %% 2 != 0){
    stop('Error: Seed must have an even number of digits. Please use another number.')
  }
  
  seed <- seed_num
  midsq_nums <- vector('numeric', length = n)
  
  for(i in 1:n){
    seed_sq <- format(seed^2, scientific = FALSE)
    
    while(nchar(seed_sq) < 2 * seed_len){
      seed_sq <- paste0("0", seed_sq)
    }
    
    seed_sq_len <- nchar(seed_sq)
      
    start <- seed_sq_len/2 - seed_len/2 + 1 # Add 1 since R starts counting at 1
    end <- start - 1 + seed_len
    
    seed <- as.numeric(substr(seed_sq, start, end))
    
    midsq_nums[i] <- seed
  }
    midsq_nums
}

gen_midsq_num(SEED_1, M)
gen_midsq_num(SEED_2, M)

# Interpretation: Using the seed 4000 only generates one random number: 0. 
#                 This would be a bad seed to start with (using this method)
#                 because squaring 4000 results in a number (16000000) with the 
#                 four middle digits all being zero, which would then get stuck 
#                 in a "loop" of only generating random numbers of zero.






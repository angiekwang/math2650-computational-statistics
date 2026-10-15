#
#   MATH 2650: Computational Statistics
#   Lab 07
#


#### (1.) 4-Fold Cross Validation ####
ratings <- USJudgeRatings
n       <- nrow(ratings)
K       <- 4  
n_k     <- floor(n/K)

PRESS_k_4   <- matrix(NA, nrow = K, ncol = 3, 
                     dimnames = list(K = 1:K, Model = c('PREP', 'ORAL', 'WRIT')))

set.seed(2018)
perm_ids    <- sample(1:n, n)

for(k in 1:K){
  if(k != K){
    kth_set <- 1:n_k + n_k*(k-1)
  } else{
    min_k   <- min(1:n_k + (k-1)*n_k)
    kth_set <- min_k:n
  }
  
  hold_ids   <- perm_ids[kth_set]
  ratings_mk <- ratings[-hold_ids,] # Data minus holdout set
  ratings_k  <- ratings[hold_ids,] # Holdout set
  
  lm_PREP_k <- lm(RTEN ~ PREP, data = ratings_mk)
  lm_ORAL_k <- lm(RTEN ~ ORAL, data = ratings_mk)
  lm_WRIT_k <- lm(RTEN ~ WRIT, data = ratings_mk)
  
  pred_PREP_k <- predict(lm_PREP_k, newdata = ratings_k)
  pred_ORAL_k <- predict(lm_ORAL_k, newdata = ratings_k)
  pred_WRIT_k <- predict(lm_WRIT_k, newdata = ratings_k)
  
  PRESS_k_4[k,1] <- sum((ratings_k$RTEN - pred_PREP_k)^2)
  PRESS_k_4[k,2] <- sum((ratings_k$RTEN - pred_ORAL_k)^2)
  PRESS_k_4[k,3] <- sum((ratings_k$RTEN - pred_WRIT_k)^2)
  
}

apply(PRESS_k_4, 2, sum)

# Using 4-Fold cross-validation, we conclude that the simple regression model using ORAL 
# (sound oral rulings) as the explanatory variable is the most predictive model, as the 
# PRESS statistic is the lowest (2.03) out of the three models.

#### (2.) Leave-One-Out CV ####

PRESS_j_LOO   <- matrix(NA, nrow = n, ncol = 3, 
                      dimnames = list(n = 1:n, Model = c('PREP', 'ORAL', 'WRIT')))

for(j in 1:n){
  ratings_mj  <- ratings[-j,] # Data minus holdout value
  ratings_j   <- ratings[j,] # Data with holdout value only

  lm_PREP_j   <- lm(RTEN ~ PREP, data = ratings_mj)
  lm_ORAL_j   <- lm(RTEN ~ ORAL, data = ratings_mj)
  lm_WRIT_j   <- lm(RTEN ~ WRIT, data = ratings_mj)
  
  pred_PREP_j <- predict(lm_PREP_j, newdata = ratings_j)
  pred_ORAL_j <- predict(lm_ORAL_j, newdata = ratings_j)
  pred_WRIT_j <- predict(lm_WRIT_j, newdata = ratings_j)
  
  PRESS_j_LOO[j,1] <- sum((ratings_j$RTEN - pred_PREP_j)^2)
  PRESS_j_LOO[j,2] <- sum((ratings_j$RTEN - pred_ORAL_j)^2)
  PRESS_j_LOO[j,3] <- sum((ratings_j$RTEN - pred_WRIT_j)^2)
}

PRESS_4fLOO <- rbind(apply(PRESS_j_LOO, 2, sum), apply(PRESS_k_4, 2, sum))
rownames(PRESS_4fLOO) <- c('LOO', '4-Fold')
PRESS_4fLOO

# Using leave-one-out cross-validation, we also conclude that the ORAL model is the most predictive, 
# with the lowest PRESS statistic (2.00) out of the three. While the absolute and relative magnitudes 
# of the PRESS statistics differ between the 4-Fold and leave-one-out cross-validation methods, we 
# come to the same conclusion about which model is the best predictor of retention score.


#### (3.) Monte Carlo CV Function ####

MC_cross <- function(P, B){ # P = number of withheld obs; B = number of samples
  
  PRESS_i_MC   <- matrix(NA, nrow = B, ncol = 3, 
                          dimnames = list(B = 1:B, Model = c('PREP', 'ORAL', 'WRIT')))
  
  
  for(i in 1:B){
    hold_ids   <- sample(1:n, P)
    ratings_mi <- ratings[-hold_ids,] # Data minus holdout set
    ratings_i  <- ratings[hold_ids,] # Holdout set
    
    lm_PREP_i <- lm(RTEN ~ PREP, data = ratings_mi)
    lm_ORAL_i <- lm(RTEN ~ ORAL, data = ratings_mi)
    lm_WRIT_i <- lm(RTEN ~ WRIT, data = ratings_mi)

    pred_PREP_i <- predict(lm_PREP_i, newdata = ratings_i)
    pred_ORAL_i <- predict(lm_ORAL_i, newdata = ratings_i)
    pred_WRIT_i <- predict(lm_WRIT_i, newdata = ratings_i)
    
    PRESS_i_MC[i,1] <- sum((ratings_i$RTEN - pred_PREP_i)^2)
    PRESS_i_MC[i,2] <- sum((ratings_i$RTEN - pred_ORAL_i)^2)
    PRESS_i_MC[i,3] <- sum((ratings_i$RTEN - pred_WRIT_i)^2)
  }
  PRESS_MC <- apply(PRESS_i_MC, 2, mean) # Mean press statistics
  PRESS_MC
}


#### (4.) Monte Carlo CV Implementation ####
set.seed(1974)
MC_cross(P = 10, B = 10000)

# Using Monte Carlo CV with 10,000 random draws, we conclude (consistent with the 
# previous two CV tests) that the ORAL model is the most predictive with the lowest mean PRESS statistic of 0.472. 

#### (5.) Monte Carlo CV Rerun ####
set.seed(1974)
MC_cross(P = 10, B = 20000)

# Using Monte Carlo CV with 20,000 random draws, we also conclude that the ORAL model is the 
# most predictive with the lowest mean PRESS statistic of 0.472. The mean PRESS statistics 
# are very similar between the Monte Carlo CV tests with 10,000 versus 20,000 draws, so 
# 10,000 samples is most likely enough samples to obtain a consistent estimate of the PRESS statistics.




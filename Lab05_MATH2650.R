#
#
#   MATH 2650: Computational Statistics
#   Lab 05
#   08 October 2026
# 
#
#

install.packages('arrangements')
library(arrangements)

#### (1.) Sleep Data ####
extra1 <- sleep$extra[sleep$group == 1]
extra2 <- sleep$extra[sleep$group == 2]

sleep_paired <- data.frame(ID = sleep$ID[sleep$group == 1], extra1 = extra1, extra2 = extra2)

sleep_paired$diff <- sleep_paired$extra2 - sleep_paired$extra1

##### (1a.) Mean Difference #####
tmean_xy <- mean(sleep_paired$diff) # Test statistic

n_pairs <- nrow(sleep_paired)
perm_signs <- unique(permutations(x = c(-1, 1), k = n_pairs, freq = c(n_pairs, n_pairs)))

N <- nrow(perm_signs)
tmean_stars <- vector('numeric', N)

for (i in 1:N) {
  signs <- perm_signs[i, ]
  diff_star <- sleep_paired$diff * signs
  tmean_stars[i] <- mean(diff_star)
}


# Initial box-and-whisker plot
plot(group, extrasleep, xlab = 'Group', ylab = 'Extra Sleep (hrs)')

# Histogram with permuted t-stats
hist(tmean_stars, freq = FALSE, xlab = 'Mean differences')
abline(v = tmean_xy, col = 'blue', lwd = 2)

# Compare t-test to our permutation test
p_mean <- mean(tmean_stars <= -abs(tmean_xy) | tmean_stars >= abs(tmean_xy))
sleep_ttest <- t.test(sleep_paired$extra1, sleep_paired$extra2, paired = TRUE) # Compare to t-test

print(paste('P-value from permutation:', round(p_mean, 5)))
print(paste('P-value from t-test:', round(sleep_ttest$p.value, 5)))


##### (1b.) Median Difference #####

tmedian_xy <- median(sleep_paired$diff) # Test statistic
tmedian_stars <- vector('numeric', N)

for (i in 1:N) {
  signs <- perm_signs[i, ]
  diff_star <- sleep_paired$diff * signs
  tmedian_stars[i] <- median(diff_star)
}

# Histogram 
hist(tmedian_stars, freq = FALSE, xlab = 'Median differences')
abline(v = tmedian_xy, col = 'blue', lwd = 2)

# Compare median permutation test to mean permutation test
p_median <- mean(tmedian_stars <= -abs(tmedian_xy) | tmedian_stars >= abs(tmedian_xy))

print(paste('P-value from median permutation:', round(p_median, 5)))
print(paste('P-value from mean permutation:', round(p_mean, 5)))

#### (2.) Graft vs Host Disease (GVHD) Data ####
install.packages('ISwR')

##### (2a.) Difference in proportions #####
gvhd <- ISwR::graft.vs.host$gvhd
dead <- ISwR::graft.vs.host$dead

gvhd_0 <- subset(ISwR::graft.vs.host, gvhd == 0)
gvhd_1 <- subset(ISwR::graft.vs.host, gvhd == 1)

# Calculate sample difference in proportions 
p_dead_0 <- mean(gvhd_0$dead == 1)
p_dead_1 <- mean(gvhd_1$dead == 1)
p_diff_xy <- p_dead_1 - p_dead_0

M <- 50000 # Large enough value of M since num of all possible combinations is 15905368710
p_diff_star <- vector('numeric', M)

set.seed(1991)

for(i in 1:M){
  dead_star <- sample(dead, length(dead), replace = FALSE)
  df_star <- data.frame(gvhd, dead_star)
  
  gvhd_0_star <- subset(df_star, gvhd == 0)
  gvhd_1_star <- subset(df_star, gvhd == 1)
  
  p_diff_i <- mean(gvhd_1_star$dead_star == 1) - mean(gvhd_0_star$dead_star == 1)
  p_diff_star[i] <- p_diff_i
  
}

hist(p_diff_star, freq = FALSE)
abline(v = p_diff_0, col = 'blue', lwd = 2)

# Calculate two-sided empirical p-value
mean(p_diff_star <= -abs(p_diff_xy) | p_diff_star >= abs(p_diff_xy))
# Note: Fail to reject null hypothesis (so no stat sig diff between means)
# P value 0.10292

##### (2b.) Difference in proportions #####
p_diff_star_2b <- vector('numeric', M)
term_1 <- (sum(gvhd_1$dead == 1) + 1)/(nrow(gvhd_1) + 2)
term_2 <- (sum(gvhd_0$dead == 1) + 1)/(nrow(gvhd_0) + 2)
p_diff_xy_2b <- term_1 - term_2

set.seed(1991)

for(i in 1:M){
  dead_star_2b <- sample(dead, length(dead), replace = FALSE)
  df_star_2b <- data.frame(gvhd, dead_star_2b)
  
  gvhd_0_star_2b <- subset(df_star_2b, gvhd == 0)
  gvhd_1_star_2b <- subset(df_star_2b, gvhd == 1)
  
  term_1_i <- (sum(gvhd_1_star_2b$dead_star_2b == 1) + 1)/(nrow(gvhd_1_star_2b) + 2)
  term_2_i <- (sum(gvhd_0_star_2b$dead_star_2b == 1) + 1)/(nrow(gvhd_0_star_2b) + 2)
  
  p_diff_i_2b <- term_1_i - term_2_i
  p_diff_star_2b[i] <- p_diff_i_2b
  
}

hist(p_diff_star_2b, freq = FALSE)
abline(v = p_diff_xy_2b, col = 'blue', lwd = 2)

mean(p_diff_star_2b <= -abs(p_diff_xy_2b) | p_diff_star_2b >= abs(p_diff_xy_2b))


#### (3.) Air Quality Data: Test Linear Relationship ####
clean_airquality <- subset(airquality, !is.na(Ozone)) # Clean missing Ozone values
temp <- clean_airquality$Temp
ozone <- clean_airquality$Ozone

rho_0     <- cor(temp, ozone, method = "pearson")
B         <- 100000 
rho_star  <- vector('numeric', length = B)
n         <- length(temp)

set.seed(1982)

for(i in 1:B){
  temp_star  <- sample(temp, n, replace = FALSE) # Permute on temperature
  cor_i         <- cor(temp_star, ozone, method = "pearson") 
  rho_star[i]   <- cor_i
}

hist(rho_star, freq = FALSE)
hist(rho_star, freq = FALSE, xlim = c(-0.4, 0.8))
abline(v = rho_0, col = 'blue', lwd = 2) # Note rho_0 is so far beyond our values generated that it doesn't show up on the histogram
rho_0

mean(rho_star <= -abs(rho_0) | rho_star >= abs(rho_0)) 



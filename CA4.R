# Load necessary libraries
library(ggplot2)
library(MVN)

# Load the dataset
data <- CA1


# Filter for Period 2
data_period2 <- subset(data, TimePeriod == 2)[, -5]  # Remove TimePeriod column

# Marginal Normality Analysis
par(mfrow = c(2, 2))  # Set up plot layout

# QQ-Plots for each variable
for (var in names(data_period2)) {
  qqnorm(data_period2[[var]], main = paste("QQ-Plot of", var), pch = 19)
  qqline(data_period2[[var]], col = "red")
}

# Shapiro-Wilk Test for each variable
shapiro_results <- sapply(data_period2, function(x) shapiro.test(x)$p.value)
print("Shapiro-Wilk Test P-Values:")
print(shapiro_results)

# Multivariate Normality Analysis
# Compute Mahalanobis distances
center <- colMeans(data_period2)
cov_matrix <- cov(data_period2)
inv_cov <- solve(cov_matrix)
mahal_distances <- apply(data_period2, 1, function(row) {
  diff <- as.matrix(row - center)
  t(diff) %*% inv_cov %*% diff
})

# Chi-Square Plot
sorted_dist <- sort(mahal_distances)
theoretical_quantiles <- qchisq(ppoints(length(sorted_dist)), df = ncol(data_period2))
plot(theoretical_quantiles, sorted_dist, main = "Chi-Square Plot for Multivariate Normality",
     xlab = "Theoretical Chi-Square Quantiles", ylab = "Observed Squared Mahalanobis Distances",
     pch = 19)
abline(a = 0, b = 1, col = "red", lty = 2)

# Mardia's Test for Multivariate Normality
mardia_test <- mvn(data_period2, mvnTest = "mardia")
print("Mardia's Test Results:")
print(mardia_test$multivariateNormality)


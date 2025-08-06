# Factor analysis
# Models observed variables(X) as linear combinations of a smaller number of unobserved latent variables(Factors) plus unique error
# Goal : Explaining correlations/covariance between the observed variables via latent factor
# latent factors explain correlation between observed variables

# X = LF+e
# X <- observed standardized variables ( px1) p <- number of observed variables
# L <- loadings (pxm)
# F <- latent factors (mx1) m <- number of latent factors
# e <- unique error (px1) uncorrelated with f

# the model implies
# Cov(X) <- LL' + U
# U <- DIAGONAL MATRIX OF UNIQUE VARIANCES

# coding for this section essentially means showing loadings

# Coding; CA

# Step 1: generate correlation matrix
air <- read.csv("~/Downloads/Air Pollution Data.csv")
cor_matrix <- cor(air)
round(cor_matrix,3)

# Step 2: Calculate principal component solution to factor model with m=1 and m=2 (loading vectors for 2 factor model)
# Number of latent factors = 2
# use eigen - similar to pca
eigen_mat <- eigen(cor_matrix)
eigenvalues <- eigen_mat$values
eigenvectors <- eigen_mat$vectors

# factor model for m = 1, m =2, means show first 2 eigenvectors loadings (this is m=2 factor model)
eigen_m1_m2 <- eigenvectors[,1:2]
# loadings
loadings <- eigen_m1_m2 %*% diag(sqrt(eigenvalues[1:2]))
colnames(loadings) <- c("Factor 1", "Factor 2")
rownames(loadings) <- colnames(air)
print(loadings)

# Obtaining max likelihood factor model (m =2) using factanal
max_likelihood <- factanal(air, factors = 2, rotation = "none")
max_likelihood_factor_model <- max_likelihood$loadings
print(max_likelihood_factor_model, cutoff = 0)

# 2D: Comparing factorization obtained by principal component and max likelihood method
# Compare pca_loadings vs mle_loadings
comparison <- cbind(loadings,max_likelihood_factor_model)
colnames(comparison) <- c("PC Factor 1", "PC Factor 2", "ML Factor 1", "ML Factor 2")
round(comparison, 3)

# 2E: Varimax rotation
varimax_method <- varimax(loadings)
varimax_rotation <- varimax_method$loadings

# 2F factor scores from = 2 maximum likelihood approach using weighted squares
# weighted least squares = F = (L'(PSI)^-1L)^-1(L'(PSI)^-1Z)
# PSI - DIAGONAL MATRIX OF SPECIFIC VARIANCES
# Z <- scale(data)

# Factor scores
psi <- diag(max_likelihood$uniquenesses)
l <- max_likelihood_factor_model
z <- scale(air)
Z <- t(z)
f_hat <- solve(t(l)%*%solve(psi)%*%l)%*%t(l)%*%solve(psi)%*%Z
t(f_hat)
round(t(f_hat), 3)



#Slides example

# Data
corr <- matrix(c(1,0.02,0.96,0.42,0.01,0.02,1,0.13,0.71,
               0.85,0.96, 0.13,1,0.50,0.11,
             0.42,0.71,0.50,1,0.79,0.01,0.85,0.11,0.79,1), 5,5,byrow=T)
# Calculate principal component solution for m=1, m =2, m=3, factor model.

eigen_matrix <- eigen(corr)
eigenvectors <- eigen_matrix$vectors
eigenvalues <- eigen_matrix$values

# m = 1 factor model 
loading_1 <- eigenvectors[,1:1] * (sqrt(eigenvalues[1:1]))
round(loading_1,3)

# m = 2 factor model
loading_2 <- eigenvectors[,1:2] %*% diag(sqrt(eigenvalues[1:2]))
round(loading_2,3)

# m = 3 factor model
loading_3 <- eigenvectors[,1:3] %*% diag(sqrt(eigenvalues[1:3]))
colnames(loading_3) <- c("Factor 1", "Factor 2", "Factor 3")
round(loading_3,3)

# Calculation max likelihood factor model for m = 2 factors
# max_likelihood factor model - use factanal function
max_likelihood_1 <- factanal(covmat = corr, correlation, factors = 2, rotation = "none")
max_likelihood_factor <- max_likelihood_1$loadings
print(max_likelihood_factor, cutoff = 0)


# varimax rotation
varimax <- varimax(loading_2)
varimax_rotation <- varimax$loadings
final <- round(varimax_rotation,3)
print(final,cutoff = 0)


# Factor scores using weighted least squares approach
psi <- diag(max_likelihood_1$uniquenesses)
Z <- t(scale(corr))
L <- max_likelihood_factor
f_hat <- solve(t(L)%*%solve(psi)%*%L)%*%t(L)%*%solve(psi)%*%Z
t(f_hat)


# calculating reproduced correlation matrix from pc solution
# Cov(X) = LL' + psi

L <- loading_2
L_tran <- t(loading_2)
h <- diag(L%*%L_tran)
psi <- diag(1-h)
cor_x <- L%*%L_tran + psi
cor_x


corre <- matrix(c(1,0.632,0.511,0.115,0.155, 0.632, 1,0.574,
0.322,0.213,0.511,0.574,1,
0.183, 0.146,0.115,0.322, 0.183,1,0.683,0.155,0.213,0.146,
0.683, 11),5,5, byrow=TRUE)

eigen_corre <- eigen(corre)
eigenvector <- eigen_corre$vectors
eigenvalue <- eigen_corre$values

# Factor model m=1, m=2 using Principal comp <- do seperately:
loading_m1 <- eigenvector[,1:1] * sqrt(eigenvalue[1])

loading_m2 <- eigenvector[,1:2] %*% diag(sqrt(eigenvalue[1:2]))
colnames(loading_m2) <- c("Factor 1", "Factor 2")
round(loading_m2,3)

# recomputing correlation matrix
L <- loading_m2
L_transpose <- t(loading_m2)
h <- diag(L%*%L_transpose)
psi <- diag(1-h)
cor_x <- L%*%L_transpose + psi
cor_x

# using max_likelihood factor model m =2 
max_like <- factanal(covmat = corre, factors = 2, rotation = "none")
max_like_factor <- max_like$loadings

print(round(max_like_factor,3), cutoff = 0)

# varimax
varimax <- varimax(loading_m2)
varimax_loading <- varimax$loadings
print(round(varimax_loading,3), cutoff = 0)


# weighted least square methods using pc m =2 <- factor score
Z <- t(scale(corre))
f_hat <- solve(L_transpose%*%solve(psi)%*%L) %*% L_transpose%*%solve(psi)%*%Z
t(f_hat)


# for PC  factor model: psi = DIAG(1- h), where h= diag(LL')
# For max_likelihood factor model <- psi = diag(ml$uniqueness)
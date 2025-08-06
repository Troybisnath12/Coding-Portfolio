# CA 11:
# Correspondence Analysis:
# Technique to analyze relationships in categorical data - specifically 2 way contingency tables (e.g. education level vs movie preference)
# Similar to PCA for categorical data, but instead of analyzing variance in numeric variables, it analyzes relationships between categories
# Start with matrix of counts <- 2 way contingency table
# Correspondence analysis from first principals <- steps

# Step 1: data

ratings <- structure(
  c(
    50, 30, 10, 1, 60, 80, 40, 2,
    40, 60, 20, 1, 10, 30, 50, 4),
  dim = c(4L, 4L),
  dimnames = list(
    c("High School", "Bachelor's", "Master's", "Doctorate"),
    c("Action", "Drama", "Comedy", "Documentary"))
)

# step 2: P <- sum(xij)/N
N <- ratings
P <- N/sum(N)

# Step 3: calculate row and column profiles
row_profile <- apply(P,1,sum) # row values/row total # row is margin 1
column_profile <- apply(P,2,sum) # column values/column total # column is margin 2
# prop.table() <- express as proportions over 1

#chisquare distance between row profiles:
# chi square distance between row profile: ri and rj
# d^2(i,j) = sum ((rik-rjk)^2/ck)
# ck - column proportion
# diagonals of matrix always 0
chi_distance_row <- dist(row_profile, method = "euclidean")
chi_matrix <- as.matrix(chi_distance_row)
heatmap(chi_matrix)

# deviations - (observed - expected)^2/expected
# expected value = Ri.Cj/n
# observed value = ratings


# principal row and principal column coordinates <- use first principles of ca or ca function, including data
# first principles steps
# step 1: calculate P
# P = data/sum(data)
N <- ratings
P <- N/sum(N)


# step 2: row masses and column masses
row_mass <- apply(P, MARGIN = 1, sum)
column_mass <- apply(P, MARGIN = 2, sum)

# step 3: calculate Drm, Dcm <- diagonal of row mass and column mass
Dr <- diag(row_mass)
Dc <- diag(col_mass)

# step 4: calculate S
# S = Dr^-1/2(P-rc')Dc^-1,2
S <- sqrt(solve(Dr))%*%(as.matrix(P)-row_mass%*%t(column_mass))%*%sqrt(solve(Dc))

# step 5: apply svd to s
svd_S <- svd(S)
U <- svd_S$u #LEFT SINGULAR VECTOR
D <- svd_S$d # SINGULAR VALUES
V <- svd_S$v # RIGHT SINGULAR VECTOR

# Step 6: calculate pricipal row and column coordinates
# row coordinates <- F = Dr^-1 %*% U %*% diag(D)
# COLUMN COORDINATES <- G = Dc^-1/2 %*% V %*%* diag(D)
 
# ca(ratings) <- plot(ca(ratings)) <- biplot
# always plot for first 2 dimensions
# transforming data for multiple correspondence analysis
# essentially, apply(data, c(1,2), sum) # gets row and column masses at once for MCA <- creates 2 way contingency table
# use mjca(new data, lambda = "adjusted)

library(ppcor)

data <- data.frame(
    x = c(1, 2, 3, 4, 5),
    y = c(2, 4, 6, 8, 10),
    z = c(3, 6, 9, 12, 15)
)

weights <- c(0.5, 1, 0.5, 1, 0.5) # Example weights

# Calculate the unweighted partial correlation
partial_corr <- pcor.test(data$x, data$y, data$z)$estimate

# Apply the weights to the partial correlation coefficient
weighted_partial_corr <- partial_corr * sqrt(weights[1] * weights[3])

weighted_partial_corr

cov.wt(data, cor = TRUE)

# create some binary data
set.seed(123)
x <- rbinom(100, 1, 0.5)
y <- rbinom(100, 1, 0.5)
z <- rbinom(100, 1, 0.5)
d <- data.table(x, y, z)

# create a weight vector
w <- runif(100)

cor(d)
cov.wt(d, cor = TRUE)$cor
cov.wt(d, wt = w, cor = TRUE)$cor |> as.table() |> as.data.frame()
test <- cov.wt(d, wt = w, cor = TRUE)$cor
data.table(
    phe1 = rownames(test)[row(test)],
    phe2 = colnames(test)[col(test)],
    cor  = c(test)
)

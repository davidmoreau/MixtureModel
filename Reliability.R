# Generate data from mixed Gaussian
library(mixtools)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)

n <- 100 #number of data points
#lam1 <- seq(0, 1, .1)
#lam2 <- 1 - lam1
#lam <- cbind(lam1, lam2)
iterations <- 1000 #repetition of the EM algorithm
mean <- c(0, 1) #mean of distributions (~ equal to meta)
#mu <- c(par$mu[1], par$mu[2])
#sigma <- c(par$sigma[1], par$sigma[2])
sigma <- c(.1, .2, .3, .4) #to cycle through in simulation


# Create function to evaluate reliability
reliability <- function(lam, sigma, mean = c(0, 1), n = 100, rep = iterations){
  sd <- sqrt(sigma)
  distrib <- rnormmix(n, lam, mean, sd)
  results.lambda <- replicate(rep, normalmixEM(distrib)$lambda)
  results.loglik <- replicate(rep, normalmixEM(distrib)$loglik)
  results.final <- apply(results.lambda, 2, sort)
  results.mean <- apply(results.final, 1, mean)
  results.dev <- c(abs(lam[1] - results.mean[1]), abs(lam[2] - results.mean[2]))
  allresults <- cbind(results.mean, results.dev, results.loglik)
  return(allresults[1,])
}

# Run simulation
data <- data.frame(
  lam1.sig25 <- reliability(lam = c(.1, .9), sigma[1]),
  lam1.sig5 <- reliability(lam = c(.1, .9), sigma[2]),
  lam1.sig75 <- reliability(lam = c(.1, .9), sigma[3]),
  lam1.sig1 <- reliability(lam = c(.1, .9), sigma[4]),
  
  lam2.sig25 <- reliability(lam = c(.2, .8), sigma[1]),
  lam2.sig5 <- reliability(lam = c(.2, .8), sigma[2]),
  lam2.sig75 <- reliability(lam = c(.2, .8), sigma[3]),
  lam2.sig1 <- reliability(lam = c(.2, .8), sigma[4]),
  
  lam3.sig25 <- reliability(lam = c(.3, .7), sigma[1]),
  lam3.sig5 <- reliability(lam = c(.3, .7), sigma[2]),
  lam3.sig75 <- reliability(lam = c(.3, .7), sigma[3]),
  lam3.sig1 <- reliability(lam = c(.3, .7), sigma[4]),
  
  lam4.sig25 <- reliability(lam = c(.4, .6), sigma[1]),
  lam4.sig5 <- reliability(lam = c(.4, .6), sigma[2]),
  lam4.sig75 <- reliability(lam = c(.4, .6), sigma[3]),
  lam4.sig1 <- reliability(lam = c(.4, .6), sigma[4]),
  
  lam5.sig25 <- reliability(lam = c(.5, .5), sigma[1]),
  lam5.sig5 <- reliability(lam = c(.5, .5), sigma[2]),
  lam5.sig75 <- reliability(lam = c(.5, .5), sigma[3]),
  lam5.sig1 <- reliability(lam = c(.5, .5), sigma[4])
)
#write.csv(data, file = "table1_results_100x1000.csv")
data <- read.csv("table1_results_100x1000_diffmean.csv")
data <- data[,2:21]
# Tidy data
names <- c("lam1.sig1",
           "lam1.sig2",
           "lam1.sig3",
           "lam1.sig4",
           
           "lam2.sig1",
           "lam2.sig2",
           "lam2.sig3",
           "lam2.sig4",
           
           "lam3.sig1",
           "lam3.sig2",
           "lam3.sig3",
           "lam3.sig4",
           
           "lam4.sig1",
           "lam4.sig2",
           "lam4.sig3",
           "lam4.sig4",
           
           "lam5.sig1",
           "lam5.sig2",
           "lam5.sig3",
           "lam5.sig4"
)
names(data) <- names
datal <- gather(data)
stat <- rep(c("mean", "diff", "log"), 20)
datal <- cbind(datal, stat)

datal$key <- as.character(datal$key)
#Then turn it back into an ordered factor
datal$key <- factor(datal$key, levels=unique(datal$key))


a <- datal %>%
  filter(stat == "log") %>%
  ggplot(aes(key, value, color = stat)) + geom_point()

b <- datal %>%
  filter(stat == "diff") %>%
  ggplot(aes(key, value, color = stat)) + geom_point()

c <- datal %>%
  filter(stat == "mean") %>%
  ggplot(aes(key, value, color = stat)) + geom_point()


grid.arrange(a, b, c, ncol = 1)


# Equivalent of function as a for loop (nested)
## Outer loop
for (j in 1:length(n)) {
  N <- n[j]
  result <- rep(NA, iterations)
  ## Inner loop
  for (i in 1:iterations) {
    mean <- 0
    sigma <- .25
    sd <- sqrt(sigma)
    distrib <- rnormmix(n, lam, mean, sd)
    results <- replicate(rep, normalmixEM(distrib)$lambda)
    results.final <- apply(results, 2, sort)
    results.mean <- apply(results.final, 1, mean)
    results.dev <- c(abs(lam[1] - results.mean[1]), abs(lam[2] - results.mean[2]))
    result[i] <- results.dev
  }
}
# Not used
distrib19 <- rnormmix(n = 100, lambda = c(.1, .9), 2, 1)
lam19 <- replicate(10, normalmixEM(distrib19)$lambda)
apply(lam19, 2, sort)

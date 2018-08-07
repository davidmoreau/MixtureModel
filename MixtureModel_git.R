#######################################################################################
# Mixture Model
# R script to accompany
# When Averaging Goes Wrong: The Case for Mixture Model Estimation in Psychological Science
# Journal of Experimental Psychology: General
# David Moreau, Auckland, 2016
#######################################################################################


# About
# This script presents the EM model and its application 

### SETUP #############################################################################

# Colors
palette <- c("#58D9BF", "#46ABE3", "#4485C5", "#3276BE", "#2A5766")

# Load packages (+ install if required)
packages <- c("psych", "ggplot2", "gridExtra", "compute.es", "metafor", "pwr", "scales", "mixtools", "dplyr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(packages, suppressPackageStartupMessages(require), warn.conflicts=F, quietly=T, character.only=T)

# Set default ggplot parameters
default <-  theme_bw() + 
  theme(plot.title=element_text(family="Arial", face="bold", size=16)) +
  theme(axis.line.x = element_line(colour = "black", size=.5, linetype = 1),
        axis.line.y = element_line(colour = "black", size=.5, linetype = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# Load data
pps <- read.csv("PPS_MelbyLervag.csv", header=T)
es.d <- pps$Std.diff.in.means
es <- pps$Hedges.s.g


### SIMULATIONS ##########################################################################

set.seed(113)

# Set parameter values
n <- 100
#lambda <- rep(14, 2)/2
lam <- 10
mu <- c(0, 2)
sigma <- rep(1, 2)

# Generate data from single Gaussian
unimod <- rnorm(n, 
                mean=(mu[1] + mu[2])/2, 
                sd=(sigma[1] + sigma[2]/2))

# Generate data from two Gaussian
#multimod <- c(rnorm(n/2, mean=mu[1]-1, sd=sigma[1]), rnorm(n/2, mean=mu[2]+1, sd=sigma[2]))

# Generate data from mixed Gaussian
multimod <- rnormmix(n, 
                     lam, #Vector of mixture probabilities, with length equal to m, the desired number of components (subpopulations). 
                     #This is assumed to sum to 1; if not, it is normalized
                     mu, 
                     sigma)

# Save in a data frame
data <- data.frame(unimod, multimod)

# Mixture model
mixmdl <- normalmixEM(multimod, k=2) ; summary(mixmdl) ; plot(mixmdl)
#mixmdl <- data.frame(mixmdl$x) #change this line to mixmdl <- data.frame(mixmdl$x) for graph

# Data frame for multimod
post <- data.frame(iteration = 1:length(mixmdl$posterior[,1]), 
                   value1 = mixmdl$posterior[,1], 
                   value2 = mixmdl$posterior[,2])
post <- post %>%
  mutate(cumsum1 = cumsum(value1), cumsum2 = cumsum(value2)) %>%
  mutate(converg1 = cumsum1 / iteration, converg2 = cumsum2 / iteration)


### PLOT BOTH DENSITY DISTRIBUTIONS (k = 2)

# Function to plot mixture comp
plotcomp <- function(x, mu, sigma, lamda) {
  lamda * dnorm(x, mu, sigma)
}

# FIGURE 2
data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = .33, colour = "gray", 
                 fill = "gray") +
  #geom_density(aes(x, ..density..), color="#FF8F00", size=1) +
  #geom_density(color="#FF8F00", size=1) +
  stat_function(geom = "line", fun = dnorm, args = list(mean = 1), colour = "black", lwd = .6, linetype = "dashed") +
  stat_function(geom = "line", fun = plotcomp,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "#FF8F00", lwd = 1.5) +
  stat_function(geom = "line", fun = plotcomp,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "#4F94CD", lwd = 1.5) +
  xlab("Effect Sizes") +
  ylab("Density") +
  xlim(-4, 5) +
  default

### FIGURE 4

posterior.plot <- ggplot(data=post, aes(x=iteration, y=converg1)) +
  geom_line(color="black", size=1.5) +
  geom_line(data=post, aes(x=iteration, y=converg2), color="#A2CD5A", size=1.5) +
  geom_hline(yintercept=.5, color="gray50", linetype="solid", size=1) +
  geom_hline(yintercept=mixmdl$lambda[1], color="gray50", linetype="dashed", size=1) +
  geom_hline(yintercept=mixmdl$lambda[2], color="gray50", linetype="dashed", size=1) +
  xlab("Iteration") + ylab("Estimated ω") +
  ggtitle("") +
  default
posterior.plot


### ADDITIONAL PLOTS ##################################


# Combine plots
#par(mfrow=c(2,1)) 

# Base plot 
plot(hist(unimod, breaks=30), 
     col = "gray",
     border = "gray",
     freq = FALSE,
     xlab = "Simulated d value", 
     main = "Single Gaussian")
lines(density(unimod),lty=2)
#sapply(1:2, plot.normal.components, mixture=es.2k)

### FIGURE 5
# Base plot 
plot(hist(multimod, breaks=30), 
     col = "gray",
     border = "gray",
     freq = FALSE,
     xlab = "Simulated d value", 
     main="", 
     ylim=c(0,1.09))
lines(density(multimod),
      lty = 2, 
      col = "#FF8F00",
      lwd = 2)
sapply(1:2, 
       plot.normal.components, 
       mixture = es.2k,
       col = "#4F94CD",
       lwd=2)


p1 <- ggplot(data, aes(multimod)) + 
  geom_histogram(aes(y=..density..), 
                 color = "gray", 
                 fill = "gray", 
                 binwidth = .3) +
  geom_density(color = "#FF8F00", size = 1) +
  #ylim(0.01, 0.25) +
  geom_vline(aes(xintercept = mean(multimod)),
             color="black", linetype = "dashed", size=1) +
  geom_vline(aes(xintercept = mean(multimod)/2),
             color="#4F94CD", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = mean(multimod)+(mean(multimod)/2)),
             color="#4F94CD", linetype="dashed", size=1) +
  xlab("Simulated effect size (d)") +
  ylab("Density") +
  ggtitle("Gaussian Mixture") +
  default

p2 <- ggplot(data, aes(unimod)) + 
  geom_histogram(aes(y=..density..), color="gray", fill="gray", binwidth=.3) +
  geom_density(color="#A2CD5A", size=1) +
  #ylim(0.01, 0.25) +
  geom_vline(aes(xintercept=mean(unimod)),
             color="black", linetype="dashed", size=1) +
  xlab("Simulated effect size (d)") +
  ylab("Density") +
  ggtitle("Single Gaussian") +
  default


grid.arrange(p2, p1, ncol=2)

# Base plot 
plot(hist(
  multimod, 
  breaks=30), 
  col="grey",
  border="grey",
  freq=FALSE,
     xlab="Simulated d value", 
  main="")
lines(density(multimod),lty=2)


# Base plot 
plot(hist(
  unimod, 
  breaks=30),
  col="grey", 
  border="grey", 
  freq=FALSE,
  xlab="Simulated d value", 
  main="")
lines(density(unimod), lty=2)



### ANALYSES #############################################################################

# Histogram of Effect Sizes
plot(hist(es,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="Effect size (Cohen's d)",main="Effect Sizes")
lines(density(es),lty=2)

# 2-component Gaussian mixture

es.2k <- normalmixEM(es, k=2, maxit=100, epsilon=0.01)
plot.normal.components <- function(mixture, component.number,...) 
  {
  curve(mixture$lambda[component.number] *
          dnorm(x,
                mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), 
        add=TRUE, ...)
}

### Plot
plot(hist(es,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="Effect size (Cohen's d)", main="Effect Sizes")
lines(density(es),lty=2, col="#FF8F00", lwd=2)
sapply(1:2,plot.normal.components,mixture=es.2k, col="#4F94CD", lwd=2)

pnormmix <- function(x, mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,
                            mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k, pnorm.from.mix, x=x)
  return(rowSums(pnorms))
}

distinct.es <- sort(unique(es))
tcdfs <- pnormmix(distinct.es, mixture=es.2k)
ecdfs <- ecdf(es)(distinct.es)
plot(tcdfs,ecdfs,xlab="Theoretical CDF", 
     ylab="Empirical CDF", 
     xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)

# Probability density function
dnormalmix <- function(x, mixture, log=FALSE) {
  lambda <- mixture$lambda
  k <- length(lambda)
  like.component <- function(x,component) {
    lambda[component]*dnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  likes <- sapply(1:k, like.component, x=x)
  d <- rowSums(likes)
  if (log) {
    d <- log(d)
  }
  return(d)
}

# Log likelihood function
loglike.normalmix <- function(x, mixture) {
  loglike <- dnormalmix(x, mixture, log=TRUE)
  return(sum(loglike))
}

n <- length(es)
data.points <- 1:n
data.points <- sample(data.points) # Permute randomly
train <- data.points[1:floor(n/2)] # First random half is training
test <- data.points[-(1:floor(n/2))] # 2nd random half is testing
candidate.component.numbers <- 2:5
loglikes <- vector(length=1+length(candidate.component.numbers))
# k=1 needs special handling
mu<-mean(es[train]) # MLE of mean
sigma <- sd(es[train])*sqrt((n-1)/n) # MLE of standard deviation
loglikes[1] <- sum(dnorm(es[test],mu,sigma,log=TRUE))
for (k in candidate.component.numbers) {
  mixture <- normalmixEM(es[train],k=k,maxit=400,epsilon=1e-2)
  loglikes[k] <- loglike.normalmix(es[test],mixture=mixture)
}

# Evaluation with v-fold cross-validation
iterations <- 10^6 #adjust number of iterations
for(i in 1:iterations) {
  n <- length(es)
  data.points <- 1:n
  data.points <- sample(data.points) 
  train <- data.points[1:floor(n/2)] 
  test <- data.points[-(1:floor(n/2))] 
  candidate.component.numbers <- 2:5
  loglikes <- vector(length=1+length(candidate.component.numbers))
  mu<-mean(es[train]) 
  sigma <- sd(es[train])*sqrt((n-1)/n)
  loglikes[1] <- sum(dnorm(es[test],mu,sigma,log=TRUE))
  for (k in candidate.component.numbers) {
    mixture <- normalmixEM(es[train],k=k,maxit=400,epsilon=1e-2)
    loglikes[k] <- loglike.normalmix(es[test],mixture=mixture)
  }
  v[i] <- cbind(mixture, loglikes)
}


### Figure 4
plot(x=1:5, y=loglikes, xlab="Number of mixture components",
     ylab="Log-likelihood on testing data")
lines(loglikes)

### FIGURE 6 ###
loglikesdata <- data.frame(x=1:5, loglikes)
p3 <- ggplot(loglikesdata, aes(x, loglikes)) + 
  geom_line(color="gray45", size=.75) +
  geom_point() +
  geom_smooth(color="#4F94CD") +
  xlab("Number of mixture components") +
  ylab("Log-likelihood (testing data)") +
  ggtitle("A") +
  default

# Plot
es.k2 <- normalmixEM(es, k=2, maxit=400, epsilon=1e-2)
plot(hist(es, breaks=101), 
     col="grey", border="gray", 
     freq=FALSE,
     xlab="Effect size (C d)",
     main="Effect Sizes", 
     xlim=c(-2.5, 2.5), 
     ylim=c(0, 1.1))
lines(density(es),lty=2, 
      col="#FF8F00")
sapply(1:2, 
       plot.normal.components, 
       mixture=es.k2)

# 3 components
es.k3 <- normalmixEM(es, 
                     k=3, 
                     maxit=400, 
                     epsilon=1e-2)
plot(hist(es, 
          breaks=101), 
     col="gray", 
     border="gray", 
     freq=FALSE,
     xlab="Effect size (C d)",
     main="3 components solution", 
     xlim=c(-2.5, 2.5), ylim=c(0, 1.1))
lines(density(es),
      lty=2, 
      col="#FF8F00")
sapply(1:2, plot.normal.components, mixture=es.k3)

distinct.es <- sort(unique(es))
tcdfs <- pnormmix(distinct.es, mixture=es.2k)
ecdfs <- ecdf(es)(distinct.es)
plot(tcdfs,ecdfs,xlab="Theoretical CDF", 
     ylab="Empirical CDF", 
     xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)


# Comp estim function
comp.estim <- function (y, 
                        x = NULL, 
                        N = NULL, 
                        max.comp = 2, 
                        B = 100, 
                        sig = 0.05, 
                        arbmean = TRUE, 
                        arbvar = TRUE, 
                        mix.type = c("logisregmix", 
                                     "multmix", 
                                     "mvnormalmix", 
                                     "normalmix", 
                                     "poisregmix",
                                     "regmix", 
                                     "regmix.mixed", 
                                     "repnormmix"), 
                        hist = TRUE, 
                        ...) 
{
  mix.type <- match.arg(mix.type)
  k = max.comp
  p = 0
  sigtest = 1
  Q.star = list()
  i = 0
  if (mix.type == "regmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          H0.fit = lm(y ~ x)
          beta = coef(H0.fit)
          Q0[i] = as.numeric(logLik(H0.fit))
          H1.fit = try(regmixEM(y = y, x = x, k = (i + 
                                                     1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(H1.fit) == "try-error") {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) 
              w = 1
            else w = 2
            beta = coef(H0.fit)
            xbeta = cbind(1, x) %*% beta
            xy.sd = sqrt(sum(H0.fit$res^2)/(length(y) - 
                                              2))
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rnorm(length(y), mean = xbeta, sd = xy.sd)
          xy.simout = lm(y.sim ~ x)
          em.out = try(regmixEM(y = y.sim, x = x, k = (i + 
                                                         1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(em.out) == "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - as.numeric(logLik(xy.simout)))
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(regmixEM(y = y, x = x, k = i, 
                                arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          H1.fit = try(regmixEM(y = y, x = x, k = (i + 
                                                     1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(H0.fit) == "try-error" || class(H1.fit) == 
              "try-error") {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE) {
              scale = H0.fit$scale
              beta = matrix(rep(H0.fit$beta, i), ncol = i)
            }
            else {
              scale = 1
            }
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) 
              w = 1
            else w = 2
          }
          beta.new = H0.fit$beta
          xbeta.new = cbind(1, x) %*% beta.new
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = sapply(1:length(y), function(i) rnorm(1, 
                                                          mean = xbeta.new, sd = ((scale * H0.fit$sigma)[wt[, 
                                                                                                            i] == 1])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = sapply(1:length(y), function(i) rnorm(1, 
                                                            mean = xbeta.new[i, (wt[, i] == 1)], 
                                                            sd = H0.fit$sigma))
            }
            else {
              y.sim = sapply(1:length(y), function(i) rnorm(1, 
                                                            mean = xbeta.new[i, (wt[, i] == 1)], 
                                                            sd = H0.fit$sigma[wt[, i] == 1]))
            }
          }
          em.out.0 = try(regmixEM(y = y.sim, x = x, k = i, 
                                  arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          em.out.1 = try(regmixEM(y = y.sim, x = x, k = (i + 
                                                           1), arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          if (class(em.out.0) == "try-error" || class(em.out.1) == 
              "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik - em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "repnormmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          dens = dnorm(y, mean = mean(y), sd = sd(y))
          Q0[i] = sum(log(dens[dens > 0]))
          H1.fit = try(repnormmixEM(x = y, k = (i + 1), 
                                    arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(H1.fit) == "try-error") {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rmvnormmix(nrow(y), mu = rep(mean(y), 
                                               ncol(y)), sigma = rep(sd(y), ncol(y)))
          dens.sim = dnorm(y.sim, mean = mean(y), sd = sd(y))
          x.simout = sum(log(dens.sim[dens.sim > 0]))
          em.out = try(repnormmixEM(x = y.sim, k = (i + 
                                                      1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(em.out) == "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - x.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(repnormmixEM(x = y, k = i, arbmean = arbmean, 
                                    arbvar = arbvar, ...), silent = TRUE)
          H1.fit = try(repnormmixEM(x = y, k = (i + 1), 
                                    arbmean = arbmean, arbvar = arbvar, ...), 
                       silen = TRUE)
          if (class(H0.fit) == "try-error" || class(H1.fit) == 
              "try-error") {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE) 
              scale = H0.fit$scale
            else scale = 1
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = sapply(1:ncol(y), function(i) rnorm(nrow(y), 
                                                        mean = H0.fit$mu, sd = ((scale * H0.fit$sigma)[wt[, 
                                                                                                          i] == 1])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = sapply(1:ncol(y), function(i) rnorm(nrow(y), 
                                                          mean = H0.fit$mu[wt[, i] == 1], sd = H0.fit$sigma))
            }
            else {
              y.sim = sapply(1:ncol(y), function(i) rnorm(nrow(y), 
                                                          mean = H0.fit$mu[wt[, i] == 1], sd = H0.fit$sigma[wt[, 
                                                                                                               i] == 1]))
            }
          }
          em.out.0 = try(repnormmixEM(x = y.sim, k = i, 
                                      arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          em.out.1 = try(repnormmixEM(x = y.sim, k = (i + 
                                                        1), arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          if (class(em.out.0) == "try-error" || class(em.out.1) == 
              "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik - em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "regmix.mixed") {
    if (is.list(y)) {
      if (length(y) != length(x)) 
        stop("Number of elements in lists for x and y must match!")
    }
    tt = sapply(1:length(x), function(i) x[[i]][, 1])
    beta = t(sapply(1:length(y), function(i) lm(y[[i]] ~ 
                                                  x[[i]])$coef))
    y = beta
    mix.type = "mvnormalmix"
  }
  if (mix.type == "mvnormalmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          y.mean = apply(y, 2, mean)
          y.cov = cov(y)
          dens = dmvnorm(y, mu = y.mean, sigma = y.cov)
          Q0[i] = sum(log(dens[dens > 0]))
          H1.fit = try(mvnormalmixEM(x = y, k = (i + 
                                                   1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(H1.fit) == "try-error") {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) 
              w = 1
            else w = 2
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rmvnorm(nrow(y), mu = apply(y, 2, mean), 
                          sigma = y.cov)
          dens.sim = dmvnorm(y.sim, mu = y.mean, sigma = y.cov)
          y.simout = sum(log(dens.sim[dens.sim > 0]))
          em.out = try(mvnormalmixEM(x = y.sim, k = (i + 
                                                       1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(em.out) == "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - y.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(mvnormalmixEM(x = y, k = i, ...), 
                       silent = TRUE)
          H1.fit = try(mvnormalmixEM(x = y, k = (i + 
                                                   1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(H0.fit) == "try-error" || class(H1.fit) == 
              "try-error") {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE) {
              H0.fit$mu = lapply(1:i, function(l) H0.fit$mu)
            }
            if (arbvar == FALSE) {
              H0.fit$sigma = lapply(1:i, function(l) H0.fit$sigma)
            }
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j <- 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(nrow(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = t(sapply(1:nrow(y), function(i) rmvnorm(1, 
                                                            mu = H0.fit$mu, sigma = H0.fit$sigma[wt[, 
                                                                                                    i] == 1][[1]])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = t(sapply(1:nrow(y), function(i) rmvnorm(1, 
                                                              mu = H0.fit$mu[wt[, i] == 1][[1]], sigma = H0.fit$sigma)))
            }
            else {
              y.sim = t(sapply(1:nrow(y), function(i) rmvnorm(1, 
                                                              mu = H0.fit$mu[wt[, i] == 1][[1]], sigma = H0.fit$sigma[wt[, 
                                                                                                                         i] == 1][[1]])))
            }
          }
          em.out.0 = try(mvnormalmixEM(x = y.sim, k = i, 
                                       arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          em.out.1 = try(mvnormalmixEM(x = y.sim, k = (i + 
                                                         1), arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          if (class(em.out.0) == "try-error" || class(em.out.1) == 
              "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik - em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "normalmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          dens = dnorm(y, mean = mean(y), sd = sd(y))
          Q0[i] = sum(log(dens[dens > 0]))
          H1.fit = try(normalmixEM(x = y, k = (i + 1), 
                                   arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(H1.fit) == "try-error") {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rnorm(length(y), mean = mean(y), sd = sd(y))
          dens.sim = dnorm(y.sim, mean = mean(y), sd = sd(y))
          x.simout = sum(log(dens.sim[dens.sim > 0]))
          em.out = try(normalmixEM(x = y.sim, k = (i + 
                                                     1), arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(em.out) == "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - x.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(normalmixEM(x = y, k = i, arbmean = arbmean, 
                                   arbvar = arbvar, ...), silent = TRUE)
          H1.fit = try(normalmixEM(x = y, k = (i + 1), 
                                   arbmean = arbmean, arbvar = arbvar, ...), 
                       silent = TRUE)
          if (class(H0.fit) == "try-error" || class(H1.fit) == 
              "try-error") {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE) 
              scale = H0.fit$scale
            else scale = 1
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = sapply(1:length(y), function(i) rnorm(1, 
                                                          mean = H0.fit$mu, sd = ((scale * H0.fit$sigma)[wt[, 
                                                                                                            i] == 1])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = sapply(1:length(y), function(i) rnorm(1, 
                                                            mean = H0.fit$mu[(wt[, i] == 1)], 
                                                            sd = H0.fit$sigma))
            }
            else {
              y.sim = sapply(1:length(y), function(i) rnorm(1, 
                                                            mean = H0.fit$mu[(wt[, i] == 1)], 
                                                            sd = H0.fit$sigma[wt[, i] == 1]))
            }
          }
          em.out.0 = try(normalmixEM(x = y.sim, k = i, 
                                     arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          em.out.1 = try(normalmixEM(x = y.sim, k = (i + 
                                                       1), arbmean = arbmean, arbvar = arbvar, ...), 
                         silent = TRUE)
          if (class(em.out.0) == "try-error" || class(em.out.1) == 
              "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik - em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "multmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          m = apply(y, 1, sum)
          n.i = apply(y, 2, sum)
          theta = n.i/sum(n.i)
          Q0[i] = sum(log(exp(apply(y, 1, ldmult, theta = theta))))
          H1.fit = try(multmixEM(y = y, k = (i + 1), 
                                 ...), silent = TRUE)
          if (class(H1.fit) == "try-error") {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = matrix(0, ncol = ncol(y), nrow = nrow(y))
          for (l in 1:length(m)) {
            y.sim[l, ] <- rmultinom(1, size = m[l], prob = theta)
          }
          theta.sim = apply(y.sim, 2, sum)/sum(apply(y.sim, 
                                                     2, sum))
          y.simout = sum(log(exp(apply(y.sim, 1, ldmult, 
                                       theta = theta))))
          em.out = try(multmixEM(y = y.sim, k = (i + 
                                                   1), ...), silent = TRUE)
          if (class(em.out) == "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - y.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(multmixEM(y = y, k = i, ...), 
                       silent = TRUE)
          H1.fit = try(multmixEM(y = y, k = (i + 1), 
                                 ...), silent = TRUE)
          if (class(H0.fit) == "try-error" || class(H1.fit) == 
              "try-error") {
            w = 1
          }
          else {
            theta = H0.fit$theta
            Q0[i] = H0.fit$loglik
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(nrow(y), size = 1, prob = H0.fit$lambda)
          y.sim = t(sapply(1:nrow(y), function(i) rmultinom(1, 
                                                            size = n.i[i], prob = H0.fit$theta[(wt[, 
                                                                                                   i] == 1), ])))
          new.y.sim = 0
          em.out.0 = try(multmixEM(y = new.y.sim, k = i, 
                                   ...), silent = TRUE)
          em.out.1 = try(multmixEM(y = new.y.sim, k = (i + 
                                                         1), ...), silent = TRUE)
          if (class(em.out.0) == "try-error" || class(em.out.1) == 
              "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik - em.out.0$loglik)
            Q.star[[i]][j]
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "logisregmix") {
    if (is.null(N)) 
      stop("Number of trials must be specified!")
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    logit <- function(x) 1/(1 + exp(-x))
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          H0.fit = glm(cbind(y, N - y) ~ x, family = binomial())
          Q0[i] = logLik(H0.fit)
          H1.fit = try(logisregmixEM(y = y, x = x, N = N, 
                                     k = (i + 1), ...), silent = TRUE)
          if (class(H1.fit) == "try-error") {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = coef(H0.fit)
            xbeta = cbind(1, x) %*% beta
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rbinom(length(y), size = N, prob = logit(xbeta))
          xy.simout = glm(cbind(y.sim, N - y.sim) ~ x, 
                          family = binomial())
          em.out = try(logisregmixEM(y = y.sim, x = x, 
                                     N = N, k = (i + 1), ...), silent = TRUE)
          if (class(em.out) == "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - logLik(xy.simout))
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(logisregmixEM(y = y, x = x, N = N, 
                                     k = i, ...), silent = TRUE)
          H1.fit = try(logisregmixEM(y = y, x = x, N = N, 
                                     k = (i + 1), ...), silent = TRUE)
          if (class(H0.fit) == "try-error" || class(H1.fit) == 
              "try-error") {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = H0.fit$beta
            xbeta = cbind(1, x) %*% beta
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          y.sim = sapply(1:length(y), function(i) rbinom(1, 
                                                         size = N[i], prob = logit(xbeta)[, (wt[, 
                                                                                                i] == 1)]))
          em.out.0 = try(logisregmixEM(y = y.sim, x = x, 
                                       N = N, k = i, ...), silent = TRUE)
          em.out.1 = try(logisregmixEM(y = y.sim, x = x, 
                                       N = N, k = (i + 1), ...), silent = TRUE)
          if (class(em.out.0) == "try-error" || class(em.out.1) == 
              "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik - em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "poisregmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          H0.fit = glm(y ~ x, family = poisson())
          Q0[i] = logLik(H0.fit)
          H1.fit = try(poisregmixEM(y = y, x = x, k = (i + 
                                                         1), ...), silent = TRUE)
          if (class(H1.fit) == "try-error") {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = coef(H0.fit)
            xbeta = cbind(1, x) %*% beta
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rpois(length(y), lambda = exp(xbeta))
          xy.simout = glm(y.sim ~ x, family = poisson())
          em.out = try(poisregmixEM(y = y.sim, x = x, 
                                    k = (i + 1), ...), silent = TRUE)
          if (class(em.out) == "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - logLik(xy.simout))
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(poisregmixEM(y = y, x = x, k = i, 
                                    ...), silent = TRUE)
          H1.fit = try(poisregmixEM(y = y, x = x, k = (i + 
                                                         1), ...), silent = TRUE)
          if (class(H0.fit) == "try-error" || class(H1.fit) == 
              "try-error") {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = H0.fit$beta
            xbeta = cbind(1, x) %*% beta
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          y.sim = sapply(1:length(y), function(i) rpois(1, 
                                                        lambda = exp(xbeta)[, (wt[, i] == 1)]))
          em.out.0 = try(poisregmixEM(y = y.sim, x = x, 
                                      k = i, ...), silent = TRUE)
          em.out.1 = try(poisregmixEM(y = y.sim, x = x, 
                                      k = (i + 1), ...), silent = TRUE)
          if (class(em.out.0) == "try-error" || class(em.out.1) == 
              "try-error") {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik - em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (hist) {
    if (length(p) == 2) {
      par(mfrow = c(1, 2))
      for (i in 1:length(p)) {
        hist(Q.star[[i]], xlab = c("Bootstrap Likelihood", 
                                   "Ratio Statistic"), 
             main = paste(i, "versus", i + 1, "Components"))
        segments(obs.Q[i], 0, obs.Q[i], B, col = 3, lwd = 2)
      }
    }
    else {
      g = ceiling(sqrt(length(p)))
      par(mfrow = c(g, g))
      for (i in 1:length(p)) {
        hist(Q.star[[i]], xlab = c("Bootstrap Likelihood", 
                                   "Ratio Statistic"), 
             main = paste(i, "versus", i + 1, "Components"))
        segments(obs.Q[i], 0, obs.Q[i], B, col = 2, lwd = 2)
      }
    }
  }
  if (p[length(p)] < sig) {
    cat("Diagnostic: Select max.", length(p) + 1, "Component(s)", 
        "\n")
  }
  else {
    cat("Diagnostic: Select max.", length(p), "Component(s)", "\n")
  }
  list(p.values = p, log.lik = Q.star, obs.log.lik = obs.Q)
}


es.boot <- comp.estim(es, max.comp=5, B=length(es), mix.type="normalmix",
                      maxit=400, epsilon=1e-2)

boot.return <- data.frame(Observation=1:length(es), 
                          y1=es.boot$log.lik[[1]], 
                          y2=es.boot$log.lik[[2]], 
                          y3=es.boot$log.lik[[3]])

sample <- boot.return[sample(nrow(boot.return), 100, replace = F), ] #whole data set (n = 854) is difficult to visualize
#instead, plotting random subset of boot.return
p4 <- ggplot(sample, aes(Observation, y1)) +
  geom_line(color="gray45", size=1) +
  geom_line(data=sample, aes(Observation, y2), color="#FF8F00", size=1) +
  #geom_line(data=boot.return, aes(Iteration, y3), color="blue") 
  ylab("Log-likelihood (posterior estimates)") +
  ggtitle("B") +
  default

### Figure 6. Decision process and final log-likelihoods ###  
p3p4 <- grid.arrange(p3, p4, ncol=2)

### Comparison between parametric and semi-parametric methods

# Multimix
theta4 <- matrix(runif(56), ncol = 14)
theta3 <- theta4[1:3,]
theta2 <- theta4[1:2,]

mixmdl.es <- multmixEM(es, lambda = rep(1, 3)/3, theta = theta3)

# Semi-par vs par
sempar <- spEMsymloc(es, mu0 = c(-1, 4)); sempar$lambdahat
par <- normalmixEM(es); par$lambda

# Function to plot components
plotcomp <- function(x, mu, sigma, lamda) {
  lamda * dnorm(x, mu, sigma)
}

# Plot and compare
data.frame(x = es) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = .05, colour = "gray", 
                 fill = "gray") +
  geom_vline(xintercept = sempar$lambdahat[1], linetype = "dashed", color = palette[1]) +
  geom_vline(xintercept = sempar$lambdahat[2], linetype = "dashed", color = palette[4]) +
  stat_function(geom = "line", fun = dnorm, args = list(mean = mean(es)), 
                colour = "black", lwd = .6, linetype = "dashed") +
  stat_function(geom = "line", fun = plotcomp,
                args = list(par$mu[1], par$sigma[1], lam = par$lambda[1]),
                colour = palette[4], lwd = 1.5) +
  stat_function(geom = "line", fun = plotcomp,
                args = list(par$mu[2], par$sigma[2], lam = par$lambda[2]),
                colour = palette[1], lwd = 1.5) +
  xlab("Effect Size") +
  ylab("Density") +
  ggtitle("B") +
  default


### Reliability
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
reliability <- function(lam, sigma, mean = c(0, 1), n = 100, rep = iterations, seed = 111){
  sd <- sqrt(sigma)
  set.seed(seed)
  results.lambda <- replicate(rep, normalmixEM(rnormmix(n, lam, mean, sd))$lambda)
  set.seed(seed)
  results.loglik <- replicate(rep, normalmixEM(rnormmix(n, lam, mean, sd))$loglik)
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
write.csv(data, file = "table1_results_100x1000.csv")
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


### Generalizability
metadata <- read.csv("Heterogeneity_Psych_Bull_1990_2013.csv")

# Remove missing values and 0
metadata.noNA <- metadata[-c(which(is.na(metadata$tau.2))), ] 
metadata.noNA.no0 <- metadata.noNA %>%
  filter(!tau.2 == 0) %>%
  filter(!tau == 0)

# Plot
ggplot(metadata.noNA.no0, aes(tau.2)) + geom_histogram(binwidth = .01)
ggplot(metadata.noNA.no0, aes(tau)) + geom_histogram(binwidth = .01)

# Count & describe
summary(metadata.noNA$tau)
length(which(metadata.noNA$tau >= .5))

# Mixture components (note that this is on heterogeneity estimates tau, not es)
comp.estim(metadata.noNA$tau,  
           mix.type="normalmix")

#######################################################################################
# When Power is Not Enough: Mixture Model Estimation in Replications and Meta-analyses
# David Moreau, Auckland, 2016
#######################################################################################


# ABOUT
# This script presents the EM model and its application. 
# More information at: https://github.com/davidmoreau/MixtureModel
# http://tinyheero.github.io/2016/01/03/gmm-em.html#checking-for-convergence
# http://rstudio-pubs-static.s3.amazonaws.com/1001_3177e85f5e4840be840c84452780db52.html

### SETUP ###

# Colors
palette <- c("#58D9BF", "#46ABE3", "#4485C5", "#3276BE", "#2A5766") #(from lighter to darker)

# Load packages (+ install if required)
packages <- c("psych", "ggplot2", "gridExtra", "compute.es", "metafor", "pwr", "scales", "mixtools", "dplyr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(packages, suppressPackageStartupMessages(require), warn.conflicts=F, quietly=T, character.only=T)

# Set default ggplot parameters
default <-  theme_bw() + 
  theme(plot.title=element_text(family="Arial", face="bold", size=20)) +
  theme(axis.line.x = element_line(colour = "black", size=.5, linetype = 1),
        axis.line.y = element_line(colour = "black", size=.5, linetype = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position="none")

###########################################################################################################################



################################
### LIMITS OF POWER ANALYSIS ###
################################

# Graph power and alpha for a given effect size
m1 <- 0  #mu H0
sd1 <- 1.5 #sigma H0
m2 <- 2.2 #mu HA
sd2 <- 1.5 #sigma HA
z_crit <- qnorm(1-(0.05/2), m1, sd1)

min1 <- m1-sd1*4
max1 <- m1+sd1*4
min2 <- m2-sd2*4
max2 <- m2+sd2*4 

# Sequence
x <- seq(min(min1,min2), max(max1, max2), .01)
# Normal distance 1
y1 <- dnorm(x, m1, sd1)
df1 <- data.frame("x" = x, "y" = y1)
#  Normal distance 2
y2 <- dnorm(x, m2, sd2)
df2 <- data.frame("x" = x, "y" = y2)

# Polygon (alpha)
y.poly <- pmin(y1,y2)
poly1 <- data.frame(x=x, y=y.poly)
poly1 <- poly1[poly1$x >= z_crit, ] 
poly1<-rbind(poly1, c(z_crit, 0))  # add lower-left corner

# Polygon (beta)
poly2 <- df2
poly2 <- poly2[poly2$x <= z_crit,] 
poly2<-rbind(poly2, c(z_crit, 0))  # add lower-left corner

# Polygon (1-beta)
poly3 <- df2
poly3 <- poly3[poly3$x >= z_crit,] 
poly3 <-rbind(poly3, c(z_crit, 0))  # add lower-left corner

# Merge polygons 
poly1$id <- 3 # alpha, give it the highest number to make it the top layer
poly2$id <- 2 # beta
poly3$id <- 1 # power; 1 - beta
poly <- rbind(poly1, poly2, poly3)
poly$id <- factor(poly$id,  labels=c("power","beta","alpha"))

# Plot power
p1 <- ggplot(poly, aes(x,y, fill=id, group=id)) +
  geom_polygon(show.legend=FALSE, alpha=I(8/10)) +
  geom_line(data=df1, aes(x,y, color="H0", group=NULL, fill=NULL), size=1.5, show.legend=FALSE) + 
  geom_line(data=df2, aes(color="HA", group=NULL, fill=NULL),size=1.5, show.legend=FALSE) +
  geom_vline(xintercept = z_crit, size=1, linetype="longdash") +
  scale_color_manual("Group", 
                     values= c("HA" = "#FF8F00","H0" = "#4F94CD")) +
  scale_fill_manual("test", values= c("alpha" = palette[1],"beta" = "#FFE4B5","power"=palette[2])) +
  annotate("text", label="beta", x=1.4, y=0.08, parse=T, size=4) +
  annotate("text", label="frac(alpha,2)", x=3.2, y=0.01, parse=T, size=4) +
  annotate("text", label="1-beta", x=3.8, y=0.08, parse=T, size=4) +
  annotate("text", label="H[0]", x=m1, y=0.28, parse=T, size=7) +
  annotate("text", label="H[1]", x=m2, y=0.28, parse=T, size=7) +
  ggtitle("") +
  # remove some elements
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=22)); p1


## PPV

# Function to calculate PPV
ppv <- function(obs, power, alpha=.05, prevalence){
  effect <- prevalence * obs
  no.effect <- obs - effect
  true.pos <- power * effect
  false.pos <- no.effect * alpha
  ppv <- true.pos / (false.pos + true.pos)
  return(ppv)
}


# Df for PPV plot
data <- data.frame(vector <- seq(.0005, 1, .001), 
                    ppv.line.2 <- ppv(1000, .2, .05, vector),
                    ppv.line.4 <- ppv(1000, .4, .05, vector),
                    ppv.line.6 <- ppv(1000, .6, .05, vector),
                    ppv.line.8 <- ppv(1000, .8, .05, vector),
                    ppv.line.1 <- ppv(1000, 1, .05, vector))

# Plot PPV
p2 <- ggplot(data, aes(vector, ppv.line.2), color=palette[1], size=1.5) +
  geom_line(color=palette[1], size=1.5) +
  geom_line(data=data, aes(vector, ppv.line.4), color=palette[2], size=1.5) +
  geom_line(data=data, aes(vector, ppv.line.6), color=palette[3], size=1.5) +
  geom_line(data=data, aes(vector, ppv.line.8), color=palette[4], size=1.5) +
  geom_line(data=data, aes(vector, ppv.line.1), color=palette[5], size=1.5) +
  xlab("ω") + ylab("PPV") +
  ggtitle("") +
  default


### FIGURE 1
bothplots <- grid.arrange(p1, p2, ncol=2); bothplots

###########################################################################################################################


###################
### SIMULATIONS ###
###################


set.seed(113)

# Set parameter values
n <- 100
lam <- 10 #lambda <- rep(14, 2)/2
mu <- c(0, 2)
sigma <- rep(1, 2)

# Generate data from single Gaussian
unimod <- rnorm(n, 
                mean=(mu[1] + mu[2])/2, 
                sd=(sigma[1] + sigma[2]/2))

# Generate data from two Gaussian
bimod <- c(rnorm(n/2+15, mean=mu[1], sd=sigma[1]), rnorm(n/2-15, mean=mu[2], sd=sigma[2]))

# Could also be generated via mixed Gaussian
multimod <- rnormmix(n, 
                     lam, #Vector of mixture probabilities, with length equal to m, the desired number of components (subpopulations). 
                     #This is assumed to sum to 1; if not, it is normalized
                     mu, 
                     sigma)

# Save in a data frame
data <- data.frame(unimod, bimod)

# Mixture model
mixmdl <- normalmixEM(bimod, k=2) ; summary(mixmdl) ; plot(mixmdl)
#mixmdl <- data.frame(mixmdl$x) #change this line to mixmdl <- data.frame(mixmdl$x) for graph

# Data frame for bimod
post <- data.frame(iteration = 1:length(mixmdl$posterior[,1]), 
                   value1 = mixmdl$posterior[,1], 
                   value2 = mixmdl$posterior[,2])
post <- post %>%
  mutate(cumsum1 = cumsum(value1), cumsum2 = cumsum(value2)) %>%
  mutate(converg1 = cumsum1 / iteration, converg2 = cumsum2 / iteration)
esplot <- c(-2,0,2,4,6)

### PLOT DENSITY DISTRIBUTIONS FOR BIMODAL (k = 2)

# Create function to plot mixture components
plotcomp <- function(x, mu, sigma, lamda) {
  lamda * dnorm(x, mu, sigma)
}

### FIGURE 2
data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = .2, colour = "gray", 
                 fill = "gray") +
  #geom_density(aes(x, ..density..), color = palette[1], size=1) +
  stat_function(geom = "line", fun = dnorm, args = list(mean = 1), colour = "black", lwd = .6, linetype = "dashed") +
  stat_function(geom = "line", fun = plotcomp,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = palette[1], lwd = 1.5) +
  stat_function(geom = "line", fun = plotcomp,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = palette[2], lwd = 1.5) +
  xlab("Effect Size") +
  ylab("Density") +
  xlim(-3, 6) +
  scale_x_continuous(breaks = esplot, labels = c(-1, 0, 1, 2, 3)) +
  default

### FIGURE 4

posterior.plot <- ggplot(data=post, aes(x=iteration, y=converg1)) +
  geom_line(color=palette[1], size=1.5) +
  geom_line(data=post, aes(x=iteration, y=converg2), color=palette[2], size=1.5) +
  geom_hline(yintercept=.5, color="gray50", linetype="solid", size=1) +
  geom_hline(yintercept=mixmdl$lambda[1], color="gray50", linetype="dashed", size=1) +
  geom_hline(yintercept=mixmdl$lambda[2], color="gray50", linetype="dashed", size=1) +
  xlab("Iteration") + ylab("Estimated φ") +
  ggtitle("") +
  default
posterior.plot

mean(mixmdl$lambda[1])
mean(mixmdl$lambda[2])

### ADDITIONAL PLOTS ##################################


# Combine plots
#par(mfrow=c(2,1)) 

# Base plot 
plot(hist(unimod, breaks=30), col="grey",border="grey",freq=FALSE,
     xlab="Simulated g value", main="Single Gaussian")
lines(density(unimod),lty=2, col="#FF8F00", lwd=2)
#sapply(1:2, plot.normal.components, mixture=es.2k)

### FIGURE 5
# Base plot 
plot(hist(bimod, breaks=30), col="grey",border="grey",freq=FALSE,
     xlab="Simulated g value", main="Mixture", ylim=c(0,1.09))
lines(density(bimod),lty=2, col="#FF8F00", lwd=2)
sapply(1:2, plot.normal.components, mixture=es.2k, col="#4F94CD", lwd=2)


p1 <- ggplot(data, aes(bimod)) + 
  geom_histogram(aes(y=..density..), color="gray", fill="gray", binwidth=.3) +
  geom_density(color="#FF8F00", size=1) +
  #ylim(0.01, 0.25) +
  geom_vline(aes(xintercept=mean(bimod)),
             color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=mean(bimod)/2),
             color="#4F94CD", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=mean(bimod)+(mean(bimod)/2)),
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
plot(hist(multimod, breaks=30), col="grey",border="grey",freq=FALSE,
     xlab="Simulated d value", main="")
lines(density(multimod),lty=2)


# Base plot 
plot(hist(unimod, breaks=30), col="grey", border="grey", freq=FALSE,
     xlab="Simulated d value", main="")
lines(density(unimod), lty=2)


###########################################################################################################################


################
### ANALYSES ###
################


# Measurement-error Gaussian mixture helpers
get_observation_weights <- function(se, mode = "measurement_error") {
  n <- length(se)
  if (mode %in% c("measurement_error", "none")) {
    return(rep(1 / n, n))
  }
  
  raw <- switch(
    mode,
    inverse_variance = 1 / pmax(se^2, 1e-12),
    inverse_se = 1 / pmax(se, 1e-12),
    stop("Unsupported weighting mode.")
  )
  
  raw / sum(raw)
}

mixture_density <- function(x, pi, mu, tau, se = NULL, mode = "measurement_error") {
  k <- length(pi)
  out <- numeric(length(x))
  
  if (is.null(se)) {
    se <- rep(0, length(x))
  }
  
  if (length(se) == 1L) {
    se <- rep(se, length(x))
  }
  
  for (j in seq_len(k)) {
    sd_j <- if (mode == "measurement_error") {
      sqrt(tau[j]^2 + se^2)
    } else {
      rep(tau[j], length(x))
    }
    
    out <- out + pi[j] * dnorm(x, mean = mu[j], sd = sd_j)
  }
  
  pmax(out, .Machine$double.xmin)
}

.unpack_theta <- function(theta, k = 2) {
  if (k == 1) {
    pi <- 1
    mu <- theta[1]
    tau <- exp(theta[2])
    return(list(pi = pi, mu = mu, tau = tau))
  }
  
  eta <- theta[1]
  pi1 <- 1 / (1 + exp(-eta))
  pi <- c(pi1, 1 - pi1)
  mu <- theta[2:(1 + k)]
  tau <- exp(theta[(2 + k):(1 + 2 * k)])
  
  list(pi = pi, mu = mu, tau = tau)
}

negloglik_mixture <- function(theta, y, se, k = 2, mode = "measurement_error") {
  pars <- .unpack_theta(theta, k = k)
  weights <- get_observation_weights(se, mode)
  dens <- mixture_density(
    x = y,
    pi = pars$pi,
    mu = pars$mu,
    tau = pars$tau,
    se = se,
    mode = mode
  )
  
  -sum(weights * log(dens))
}

fit_precision_mixture <- function(y, se, k = 2, n_starts = 15, seed = 123,
                                  mode = "measurement_error") {
  keep <- is.finite(y) & is.finite(se) & se > 0
  y <- y[keep]
  se <- se[keep]
  
  if (length(y) < k + 2) {
    stop("Not enough usable rows after cleaning to fit the mixture.")
  }
  
  set.seed(seed)
  
  y_sd <- stats::sd(y)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 0.1
  
  y_q <- stats::quantile(y, probs = c(0.33, 0.67), na.rm = TRUE)
  starts <- vector("list", n_starts)
  
  for (s in seq_len(n_starts)) {
    if (k == 1) {
      starts[[s]] <- c(
        rnorm(1, mean = mean(y), sd = y_sd / 5),
        log(pmax(y_sd * runif(1, 0.5, 1.5), 1e-3))
      )
    } else {
      pi1 <- runif(1, 0.2, 0.8)
      mu_start <- c(
        rnorm(1, mean = y_q[1], sd = y_sd / 4),
        rnorm(1, mean = y_q[2], sd = y_sd / 4)
      )
      tau_start <- pmax(c(y_sd / 3, y_sd / 2) * runif(2, 0.5, 1.5), 1e-3)
      
      starts[[s]] <- c(
        qlogis(pi1),
        mu_start,
        log(tau_start)
      )
    }
  }
  
  fits <- lapply(starts, function(st) {
    optim(
      par = st,
      fn = negloglik_mixture,
      y = y,
      se = se,
      k = k,
      mode = mode,
      method = "BFGS",
      control = list(maxit = 3000, reltol = 1e-10)
    )
  })
  
  conv <- vapply(fits, function(z) identical(z$convergence, 0L), logical(1))
  if (!any(conv)) {
    fits_nm <- lapply(starts, function(st) {
      optim(
        par = st,
        fn = negloglik_mixture,
        y = y,
        se = se,
        k = k,
        mode = mode,
        method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-10)
      )
    })
    conv_nm <- vapply(fits_nm, function(z) identical(z$convergence, 0L), logical(1))
    fits <- c(fits, fits_nm)
    conv <- c(conv, conv_nm)
  }
  
  if (!any(conv)) {
    warning("No optimization run formally converged; using the best available fit.")
    best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  } else {
    fits_ok <- fits[conv]
    best <- fits_ok[[which.min(vapply(fits_ok, `[[`, numeric(1), "value"))]]
  }
  
  pars <- .unpack_theta(best$par, k = k)
  if (k > 1) {
    ord <- order(pars$mu)
    pars$pi <- pars$pi[ord]
    pars$mu <- pars$mu[ord]
    pars$tau <- pars$tau[ord]
  }
  
  posterior <- if (k == 1) {
    matrix(1, nrow = length(y), ncol = 1)
  } else {
    post <- vapply(seq_along(y), function(i) {
      obs_sd <- if (mode == "measurement_error") {
        sqrt(pars$tau^2 + se[i]^2)
      } else {
        pars$tau
      }
      numer <- pars$pi * dnorm(y[i], mean = pars$mu, sd = obs_sd)
      numer / sum(numer)
    }, numeric(k))
    t(post)
  }
  
  colnames(posterior) <- paste0("comp", seq_len(ncol(posterior)))
  
  structure(
    list(
      x = y,
      se = se,
      lambda = pars$pi,
      mu = pars$mu,
      sigma = pars$tau,
      posterior = posterior,
      loglik = -best$value,
      fit = best,
      weighting_mode = mode,
      weights = get_observation_weights(se, mode)
    ),
    class = "precisionmix"
  )
}


# Load meta-analysis data
pps <- read.csv("PPS_MelbyLervag.csv", header=T)
es.d <- pps$`Std diff in means`
es <- pps$`Hedges's g`
se <- pps$`Std Err`

keep <- is.finite(es) & is.finite(se) & se > 0
es <- es[keep]
se <- se[keep]
es.d <- es.d[keep]

# Mixture model
mixmdl.es <- fit_precision_mixture(es, se, k=2, n_starts=25, seed=123,
                                   mode="measurement_error")
mixmdl.es

# Data frame for es
post.es <- data.frame(iteration = 1:length(mixmdl.es$posterior[,1]), 
                   value1 = mixmdl.es$posterior[,1], 
                   value2 = mixmdl.es$posterior[,2])
post.es <- post.es %>%
  mutate(cumsum1 = cumsum(value1), cumsum2 = cumsum(value2)) %>%
  mutate(converg1 = cumsum1 / iteration, converg2 = cumsum2 / iteration)

posterior.plot.es <- ggplot(data=post.es, aes(x=iteration, y=converg1)) +
  geom_line(color=palette[1], size=1.5) +
  geom_line(data=post.es, aes(x=iteration, y=converg2), color=palette[4], size=1.5) +
  geom_hline(yintercept=.5, color="gray50", linetype="solid", size=1) +
  geom_hline(yintercept=mixmdl.es$lambda[1], color="gray50", linetype="dashed", size=1) +
  geom_hline(yintercept=mixmdl.es$lambda[2], color="gray50", linetype="dashed", size=1) +
  xlab("Iteration") + ylab("Estimated φ") +
  ggtitle("") +
  default
posterior.plot.es

mean(mixmdl.es$lambda[1])
mean(mixmdl.es$lambda[2])

### FIGURE 
data.frame(x = mixmdl.es$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = .05, colour = "gray", 
                 fill = "gray") +
  #geom_density(aes(x, ..density..), color="#FF8F00", size=1) +
  stat_function(geom = "line", fun = dnorm, args = list(mean = 1), colour = "black", lwd = .6, linetype = "dashed") +
  stat_function(geom = "line", fun = plotcomp,
                args = list(mixmdl.es$mu[1], mixmdl.es$sigma[1], lam = mixmdl.es$lambda[1]),
                colour = palette[4], lwd = 1.5) +
  stat_function(geom = "line", fun = plotcomp,
                args = list(mixmdl.es$mu[2], mixmdl.es$sigma[2], lam = mixmdl.es$lambda[2]),
                colour = palette[1], lwd = 1.5) +
  xlab("Effect Size") +
  ylab("Density") +
  default


# Histogram of Effect Sizes
plot(hist(es,breaks=101), col="gray", border="gray", freq=FALSE,
     xlab="Effect size (Cohen's d)", main="Effect Sizes")
lines(density(es), lty=2)

# Two-component Gaussian mixture with measurement-error weighting

es.2k <- fit_precision_mixture(es, se, k=2, n_starts=25, seed=456,
                               mode="measurement_error")

# Function to add Gaussian mixture components, with k vertically scaled
plot.normal.components <- function(mixture, component.number,...) {
  curve(mixture$lambda[component.number] *
          dnorm(x,mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), add=TRUE, ...)
}

### Plot
plot(hist(es,breaks=101),col="gray",border="gray",freq=FALSE,
     xlab="Effect size (Hedge's g)", main="")
lines(density(es),lty=2, col="#FF8F00", lwd=2)
sapply(1:2,plot.normal.components,mixture=es.2k, col="#4F94CD", lwd=2)

# Function to calculate the cumulative distribution function of a Gaussian mixture model
pnormmix <- function(x, mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k, pnorm.from.mix, x=x)
  return(rowSums(pnorms))
}


# Distinct values in the data
distinct.es <- sort(unique(es))
# Theoretical CDF evaluated at each distinct value
tcdfs <- pnormmix(distinct.es, mixture=es.2k)
# Empirical CDF evaluated at each distinct value
# ecdf(es) returns an object which is a _function_, suitable for application
# to new vectors
ecdfs <- ecdf(es)(distinct.es)
# Plot them against each other
plot(tcdfs,ecdfs,xlab="Theoretical CDF", ylab="Empirical CDF", xlim=c(0,1),
     ylim=c(0,1))
# Main diagonal for visual reference
abline(0,1)


# Probability density function for a Gaussian mixture
dnormalmix <- function(x, mixture, se = NULL, mode = "measurement_error", log=FALSE) {
  lambda <- mixture$lambda
  d <- mixture_density(x, pi = lambda, mu = mixture$mu, tau = mixture$sigma,
                       se = se, mode = mode)
  if (log) {
    d <- log(d)
  }
  return(d)
}

# Log likelihood function for a Gaussian mixture
loglike.normalmix <- function(x, mixture, se = NULL, mode = "measurement_error") {
  loglike <- dnormalmix(x, mixture, se = se, mode = mode, log=TRUE)
  return(sum(loglike))
}

# Evaluate various numbers of Gaussian components by data-set splitting
n <- length(es)
data.points <- 1:n
data.points <- sample(data.points) # Permute randomly
train <- data.points[1:floor(n/2)] # First random half is training
test <- data.points[-(1:floor(n/2))] # 2nd random half is testing
candidate.component.numbers <- 2:5
loglikes <- vector(length=1+length(candidate.component.numbers))
# k=1 needs special handling
es.k1 <- fit_precision_mixture(es[train], se[train], k=1, n_starts=15, seed=789,
                               mode="measurement_error")
loglikes[1] <- sum(log(mixture_density(es[test], pi=es.k1$lambda, mu=es.k1$mu,
                                       tau=es.k1$sigma, se=se[test],
                                       mode="measurement_error")))
for (k in candidate.component.numbers) {
  mixture <- fit_precision_mixture(es[train], se[train], k=k, n_starts=20,
                                   seed=789 + k, mode="measurement_error")
  loglikes[k] <- loglike.normalmix(es[test], mixture=mixture, se=se[test],
                                   mode="measurement_error")
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
es.k2 <- fit_precision_mixture(es, se, k=2, n_starts=25, seed=999,
                               mode="measurement_error")
plot(hist(es, breaks=101), col="grey", border="grey", freq=FALSE,
     xlab="Effect size (g)", main="Effect Sizes", xlim=c(-2.5, 2.5), ylim=c(0, 1.1))
lines(density(es),lty=2, col="#FF8F00")
sapply(1:2, plot.normal.components, mixture=es.k2)

# 3 components
es.k3 <- fit_precision_mixture(es, se, k=3, n_starts=25, seed=1001,
                               mode="measurement_error")
plot(hist(es, breaks=101), col="grey", border="grey", freq=FALSE,
     xlab="Effect size (g)", main="3 components solution", xlim=c(-2.5, 2.5), ylim=c(0, 1.1))
lines(density(es),lty=2, col="#FF8F00")
sapply(1:3, plot.normal.components, mixture=es.k3)


distinct.es <- sort(unique(es))
tcdfs <- pnormmix(distinct.es, mixture=es.2k)
ecdfs <- ecdf(es)(distinct.es)
plot(tcdfs,ecdfs,xlab="Theoretical CDF", ylab="Empirical CDF", xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)


# Function to estimate number of components
source("comp_estim.R")
es.comp <- comp.estim(es, max.comp=5, B=length(es), mix.type="normalmix",
                     maxit=400, epsilon=1e-2)

return.comp <- data.frame(Observation=1:length(es), 
                          y1=es.comp$log.lik[[1]], 
                          y2=es.comp$log.lik[[2]], 
                          y3=es.comp$log.lik[[3]])

sample <- return.comp[sample(nrow(return.comp), 100, replace = F), ] #whole data set (n = 854) is difficult to visualize
#instead, plotting random subset of boot.return
p4 <- ggplot(sample, aes(Observation, y1)) +
  geom_line(color="gray45", size=1) +
  geom_line(data=sample, aes(Observation, y2), color="#FF8F00", size=1) +
  #geom_line(data=return.comp, aes(Iteration, y3), color="blue") 
  ylab("Log-likelihood (posterior estimates)") +
  ggtitle("B") +
  default

### Figure 6. Decision process and final log-likelihoods ###  
p3p4 <- grid.arrange(p3, p4, ncol=2)

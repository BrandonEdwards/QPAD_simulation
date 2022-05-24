library(bSims)
library(detect)
library(Distance)
library(tidyverse)
library(jagsUI)

rm(list=ls())

setwd("~/iles_ECCC/Landbirds/QPAD_simulation/scripts")

# ---------------------------------
# Simulation parameters
# ---------------------------------

tau_median <- 1

nproj <- 10                                               # Data comes from 10 different projects

density <- runif(nproj,5,20) %>% sort()                   # density of males per ha in each project (controls sample size)
tau <- rlnorm(nproj,log(tau_median),0.3) %>% sort()       # EDR in 100m units, where probability of detection = exp(-(d/tau)^2)
phi <- rep(1,nproj)                                       # singing rate.  probability that event has occurred by time t = F(t) = 1-exp(-t*phi)

Y_list <- vector(mode = "list", length = nproj)
D_list <- vector(mode = "list", length = nproj)

# Loop through projects and simulate datasets
for (i in 1:nproj){
  
  # ---------------------------------
  # Transcription parameters differ among projects.  
  # Also, projects 
  # ---------------------------------
  tint <- c(1, 2, 3, 4, 5)
  
  # First 8 projects have narrow bins
  if (i <= 8) rint <- c(seq(0.25,2,0.25), Inf)
  
  # Last 2 projects have wide bins (but also have the most data)
  if (i > 8) rint <- c(1, 2, Inf)
  
  # ---------------------------------
  # Simulate landscape
  # ---------------------------------
  l1 <- bsims_init()
  p1 <- bsims_populate(l1, density = density[i]) # by default landscape is 1 km2 (100 ha)
  
  get_abundance(p1)
  get_density(p1) # Males per ha
  
  # ---------------------------------
  # Simulate bird behaviour (singing)
  # ---------------------------------
  
  e1 <- bsims_animate(p1, vocal_rate = phi[i], move_rate = 0)
  
  # ---------------------------------
  # Simulate detection process
  # ---------------------------------
  
  d1 <- bsims_detect(e1, tau=tau[i])
  
  # ------------------------------------
  # Transcribe
  # ------------------------------------
  
  x <- bsims_transcribe(d1, tint=tint, rint=rint)
  y <- get_table(x, "removal")
  
  # ------------------------------------
  # Prepare data for analysis
  # ------------------------------------
  
  # Distance model (q - detectability)
  Y = matrix(rowSums(y),1)
  D = matrix(rint,1)
  fit.q <- cmulti.fit(Y,D, type = "dis") # Fit using frequentist model
  tau_est <- exp(fit.q$coef)
  
  # ------------------------------------
  # Append data to full data matrices
  # ------------------------------------
  Y_list[[i]] <- rowSums(y)
  D_list[[i]] <- rint
  
}

# Maximum number of bins
nbins <- sapply(D_list, FUN = function(x)length(x))
maxbins <- max(nbins)

Ymat <- matrix(NA,nrow = nproj, ncol = maxbins)
Dmat <- matrix(NA,nrow = nproj, ncol = maxbins)

for (i in 1:nproj){
  Ymat[i,1:(nbins[i])] <- Y_list[[i]]
  Dmat[i,1:(nbins[i])] <- D_list[[i]]
}

# ---------------------------------
# Fit "naive" model that assumes tau is constant among projects
# ---------------------------------

# Model script
sink("dist_baseline.jags")
cat("

    model {

    # ------------------------------
    # Priors for EDR
    # ------------------------------
    
    tau ~ dunif(0,5)
    
    # ------------------------------
    # Calculate multinomial cell probabilities
    # ------------------------------
    
    for (i in 1:nproj){
    
      for (j in 1:(nbins[i])){
        cdf_rint[i,j] <- 1-exp(-(rint[i,j]/tau)^2)
      }
    
      p[i,1] <- cdf_rint[i,1] - 0
      
      for (j in 2:(nbins[i])){
        p[i,j] <- cdf_rint[i,j] - cdf_rint[i,j-1]
      }
    }
    
    # ------------------------------
    # Likelihood
    # ------------------------------
    
    for (i in 1:nproj){
      Y[i,1:(nbins[i])] ~ dmulti(p[i,1:(nbins[i])],N[i])
    }
    
    }
",fill = TRUE)
sink()

jags_data <- list(Y = Ymat,
                  rint = Dmat,
                  p = Dmat*NA,
                  N = rowSums(Ymat,na.rm = TRUE),
                  nproj = nproj,
                  nbins = nbins)

jags_data$p[is.na(Dmat)] <- 0

out1 <- jags(data = jags_data,
            model.file =  "dist_baseline.jags",
            parameters.to.save = c("tau"),
            inits = NULL,
            n.chains = 3,
            n.thin = 5,
            n.iter = 60000,
            n.burnin = 10000,
            parallel = TRUE)

# ---------------------------------
# Fit hierarchical model using JAGS
# ---------------------------------

# Model script
sink("dist_hierarchical.jags")
cat("

    model {

    # ------------------------------
    # Priors for EDR
    # ------------------------------
    
    tau_median ~ dunif(0,5)
    tau_sd ~ dunif(0,2)
    tau_precision <- pow(tau_sd,-2)
    
    for (i in 1:nproj){
      tau[i] ~ dlnorm(log(tau_median),tau_precision)
    }
    
    # ------------------------------
    # Calculate multinomial cell probabilities
    # ------------------------------
    
    for (i in 1:nproj){
    
      for (j in 1:(nbins[i])){
        cdf_rint[i,j] <- 1-exp(-(rint[i,j]/tau[i])^2)
      }
    
      p[i,1] <- cdf_rint[i,1] - 0
      
      for (j in 2:(nbins[i])){
        p[i,j] <- cdf_rint[i,j] - cdf_rint[i,j-1]
      }
    }
    
    # ------------------------------
    # Likelihood
    # ------------------------------
    
    for (i in 1:nproj){
      Y[i,1:(nbins[i])] ~ dmulti(p[i,1:(nbins[i])],N[i])
    }
    
    }
",fill = TRUE)
sink()

jags_data <- list(Y = Ymat,
                  rint = Dmat,
                  p = Dmat*NA,
                  N = rowSums(Ymat,na.rm = TRUE),
                  nproj = nproj,
                  nbins = nbins)

jags_data$p[is.na(Dmat)] <- 0

out2 <- jags(data = jags_data,
            model.file =  "dist_hierarchical.jags",
            parameters.to.save = c("tau_median","tau_sd","tau"),
            inits = NULL,
            n.chains = 3,
            n.thin = 5,
            n.iter = 60000,
            n.burnin = 10000,
            parallel = TRUE)

# ---------------------------------
# Compare resulting estimates of EDR from each model
# ---------------------------------

par(mfrow=c(2,1))
# Estimated EDR from naive model
hist(out1$sims.list$tau, breaks = seq(0,1000,0.01), xlim = c(0,3), border = "transparent", col = "cornflowerblue", main = "EDR estimate from naive model")
abline(v = tau_median)

# Estimated median EDR from hierarchical model
hist(out2$sims.list$tau_median, breaks = seq(0,1000,0.01), xlim = c(0,3), border = "transparent", col = "cornflowerblue", main = "Median EDR estimate from mixed model")
abline(v = tau_median)
par(mfrow=c(1,1))

# ---------------------------------
# Plot true vs estimated EDR
# ---------------------------------

lim <- range(tau,out2$mean$tau)
plot(out2$mean$tau ~ tau, pch = 19,ylim=lim,xlim=lim, ylab = "Estimated project EDR", xlab = "True project EDR")
abline(a = 0, b = 1)

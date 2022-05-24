library(bSims)
library(detect)
library(Distance)
library(ggplot2)

setwd("~/iles_ECCC/Landbirds/QPAD_simulation/scripts")

rm(list=ls())

# ---------------------------------
# Simulation parameters
# ---------------------------------
density = 30     # density of males per ha
phi <- 1         # singing rate.  probability that event has occurred by time t = F(t) = 1-exp(-t*phi)
tau <- 1         # EDR in 100m units, where probability of detection = exp(-(d/tau)^2)

tint_scenarios <- list(c(1, 2, 3, 4, 5))

rint_scenarios <- list(c(seq(0.25,1.5,0.25),Inf),
                       c(0.5, Inf),
                       c(1,Inf),
                       c(1.5,Inf))


simulation_results <- data.frame()

for (sim_rep in 1:100){
  
  # ---------------------------------
  # Simulate landscape
  # ---------------------------------
  l1 <- bsims_init()
  p1 <- bsims_populate(l1, density = density) # by default landscape is 1 km2 (100 ha)
  
  get_abundance(p1)
  get_density(p1) # Males per ha
  
  # ---------------------------------
  # Simulate bird behaviour (singing)
  # ---------------------------------
  
  e1 <- bsims_animate(p1, vocal_rate = phi, move_rate = 0)
  
  # ---------------------------------
  # Simulate detection process
  # ---------------------------------
  
  d1 <- bsims_detect(e1, tau=tau)
  
  
  for (tint in tint_scenarios){
    for (rint in rint_scenarios){
      
      # ------------------------------------
      # Expected proportions in each distance bin
      # ------------------------------------
      
      cdf_rint <- 1-exp(-(rint/tau)^2)
      p_rint <- diff(c(0,cdf_rint))
      
      # ------------------------------------
      # Transcribe
      # ------------------------------------
      
      x <- bsims_transcribe(d1, tint=tint, rint=rint)
      y <- get_table(x, "removal")
      
      # ------------------------------------
      # Analysis
      # ------------------------------------
      
      # Removal model (p - detectability)
      Y = matrix(colSums(y),1)
      D = matrix(tint,1)
      fit.p <- cmulti.fit(Y,D, type = "rem")
      
      # Distance model (q - detectability)
      Y = matrix(rowSums(y),1)
      D = matrix(rint,1)
      fit.q <- cmulti.fit(Y,D, type = "dis")
      
      tau_est <- exp(fit.q$coef)
      
      simulation_results <- rbind(simulation_results, 
                                  data.frame(rint = toString(rint),
                                             tint = toString(tint),
                                             sim_rep = sim_rep,
                                             tau_est = tau_est))
      
    }
  }

  bp <- ggplot(simulation_results, aes(x = rint, y = tau_est, fill = tint))+
    geom_boxplot()+ggtitle(max(simulation_results$sim_rep))
  
  print(bp)
  
}


# ---------------------------------------------------
# Bayesian distance sampling model
# ---------------------------------------------------

library(jagsUI)
# Model script
sink("dist.jags")
cat("

    model {

    # ------------------------------
    # Prior for EDR
    # ------------------------------
    
    # tau ~ dgamma(0.001,0.001)
    tau ~ dunif(0,5)
    
    # ------------------------------
    # Calculate multinomial cell probabilities
    # ------------------------------
    
    for (i in 1:nbins){
      cdf_rint[i] <- 1-exp(-(rint[i]/tau)^2)
    }
    
    p[1] <- cdf_rint[1]- 0
    for (i in 2:nbins){
      p[i] <- cdf_rint[i] - cdf_rint[i-1]
    }
    
    # ------------------------------
    # Likelihood
    # ------------------------------
    
    Y[1,1:nbins] ~ dmulti(p[],N)
    
    }
",fill = TRUE)
sink()

rint <- c(0.5,1,1.5,Inf)
tint <- seq(1,5)

# ------------------------------------
# Expected proportions in each distance bin
# ------------------------------------

cdf_rint <- 1-exp(-(rint/tau)^2)
diff(c(0,cdf_rint))

# ------------------------------------
# Transcribe
# ------------------------------------

x <- bsims_transcribe(d1, tint=tint, rint=rint)
y <- get_table(x, "removal")

# ------------------------------------
# Frequentist analysis
# ------------------------------------

# Removal model (p - detectability)
Y = matrix(colSums(y),1)
D = matrix(tint,1)
fit.p <- cmulti.fit(Y,D, type = "rem")

# Distance model (q - detectability)
Y = matrix(rowSums(y),1)
D = matrix(rint,1)
fit.q <- cmulti.fit(Y,D, type = "dis")
tau_est <- exp(fit.q$coef)

# ------------------------------------
# Bayesian analysis
# ------------------------------------

jags_data <- list(Y = Y,
                  N = sum(Y),
                  rint = rint,
                  nbins = length(rint))

out <- jags(data = jags_data,
             model.file =  "dist.jags",
             parameters.to.save = c("tau","cdf_rint","p_rint"),
             inits = NULL,
             n.chains = 3,
             n.thin = 5,
             n.iter = 60000,
             n.burnin = 10000,
             parallel = TRUE)

hist(out$sims.list$tau)

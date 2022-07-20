library(bSims)
library(detect)
library(Distance)
library(ggplot2)

dir.create("plots")

# ---------------------------------
# Simulation parameters
# ---------------------------------
density = 30     # density of males per ha
phi <- 1         # singing rate.  probability that event has occurred by time t = F(t) = 1-exp(-t*phi)
tau <- 1         # EDR in 100m units, where probability of detection = exp(-(d/tau)^2)
n_rep <- 10

tint_scenarios <- list(c(1, 2, 3, 4, 5))

rint_scenarios <- list(seq(0.1,1.5,0.10),
                       c(0.5, 1.0, 1.5),
                       c(1.0, 2.0),
                       c(0.5, Inf))


simulation_results <- data.frame()

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

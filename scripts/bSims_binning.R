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

rint_scenarios <- list(c(seq(0.25,1.5,0.25),Inf),
                       c(0.5, Inf),
                       c(1,Inf),
                       c(1.5,Inf))


simulation_results <- data.frame()

for (sim_rep in 1:n_rep){
  print(sim_rep)
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
    for (r_index in 1:length(rint_scenarios)){
      rint <- rint_scenarios[[r_index]]
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

  if ((sim_rep %% 10000) == 0)
  {
    bp <- ggplot(simulation_results, aes(x = rint, y = tau_est, fill = tint))+
      geom_boxplot()+ggtitle(max(simulation_results$sim_rep)) +
      ylim(0.7,1.5)+
      #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
      NULL
    
    png(filename = paste0("plots/",sim_rep, ".png"), width = 8, height = 6, res = 300, units = "in")
    print(bp)
    dev.off()    
  }

  
}


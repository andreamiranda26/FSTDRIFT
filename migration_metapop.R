library(scales)
#setwd("~/Desktop/")
setwd("C:/Users/andre/Box/New folder")

# Census size
total = 13000
N1 = total*0.5  #east
N2 = total*0.35 #west
N3 = total*0.15 #central

# Harmonic mean population sizes for each pair (set up for calculation later)
k = c(0.1, 0.3)  #Ne Scaling factors
Ne12 = (2 * (N1 * N2) / (N1 + N2)) 
Ne13 = (2 * (N1 * N3) / (N1 + N3)) 
Ne23 = (2 * (N2 * N3) / (N2 + N3)) 

# How long
times = seq(0, 50, by = 1)

# Migration rates (need to be checked, is it equal between all three?)
M12 = c(0, 0.12, 0.24)
M13 = c(0, 0.12, 0.24)
M23 = c(0, 0.12, 0.24)

# All outputs for each m and k get saved here
ALLOUTPUT = NULL

# Run calculations for each migration rate
for(m in 1:length(M12)){
  # All outputs for each k get saved here
  OUTPUT = NULL
  
  # Run calculations with all scaling factors
  for(kk in 1:length(k)){
    # Data frame to store FSTs
    FST_values = data.frame(time = times, FST12 = rep(NA, length(times)), FST13 = rep(NA, length(times)), FST23 = rep(NA, length(times)))
    FST_values$FST12[1] = 0
    FST_values$FST13[1] = 0.0004
    FST_values$FST23[1] = 0.0008
    
    # Function to update FST given migration rate, effective population size, and scaling factor
    update_FST = function(FST, mig, N_e,k) {
      FST_new = 1 - (1 - FST) * exp(-1/((2*N_e*k)*(1-mig)))
      return(FST_new)
    }
    
    # Run calculation
    for (i in 2:length(times)) {
      FST_values$FST12[i] = update_FST(FST_values$FST12[i-1], M12[m], Ne12,k[kk])
      FST_values$FST13[i] = update_FST(FST_values$FST13[i-1], M13[m], Ne13,k[kk])
      FST_values$FST23[i] = update_FST(FST_values$FST23[i-1], M23[m], Ne23,k[kk])
    }
    FST_values$k = rep(k[kk], nrow(FST_values))
    
    # Record outputs
    OUTPUT = rbind(OUTPUT, FST_values)
  }
  #save all output
  temp = cbind(OUTPUT, rep(M12[m], nrow(OUTPUT)))
  colnames(temp)[6] = "migrate"
  ALLOUTPUT = rbind(ALLOUTPUT, temp)
  
  # Save all data into data frame
  FST_values = as.data.frame(OUTPUT) #note that migration rate is defined in ALLOUTPUT but not here
  
  # Plot output - note that if you increase 'times' need to update here
  pdf(file=paste("FST_", k[1], "to", k[2], "_", M12[m], "mig.pdf", sep=""), width=5, height=5, bg="transparent")
  plot(-100,-100, xlab = 'Time (generations)', ylab = 'FST', main = '', xlim=c(0,50), ylim=c(0,0.1))
  abline(h=0.05, lwd=2, col="grey50", lty=2)
  polygon(x=c(FST_values$time[1:51], rev(FST_values$time[1:51])), y=c(FST_values$FST12[1:51], rev(FST_values$FST12[52:102])), col=alpha('chartreuse3', 0.3), border=F)
  polygon(x=c(FST_values$time[1:51], rev(FST_values$time[1:51])), y=c(FST_values$FST13[1:51], rev(FST_values$FST13[52:102])), col=alpha('darkorchid3', 0.3), border=F)
  polygon(x=c(FST_values$time[1:51], rev(FST_values$time[1:51])), y=c(FST_values$FST23[1:51], rev(FST_values$FST23[52:102])), col=alpha('darkorange3', 0.2), border=F)
  lines(FST_values$time[1:51], apply(cbind(FST_values$FST12[1:51], FST_values$FST12[52:102]), 1, mean), col = 'chartreuse4', lwd = 3)
  lines(FST_values$time[1:51], apply(cbind(FST_values$FST13[1:51], FST_values$FST13[52:102]), 1, mean), col = 'darkorchid3', lwd = 3)
  lines(FST_values$time[1:51], apply(cbind(FST_values$FST23[1:51], FST_values$FST23[52:102]), 1, mean), col = 'darkorange3', lwd = 3)
  dev.off()
}


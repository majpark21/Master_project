sim_phase_shifted <- function(n, noise, freq = 1, len = 100, return.wide = F){
  # A fucntion to simulate a population of n noisy sinusoidals that are phase shifted
  # n: number of sinusoids
  # freq: sampling rate (in time unit)
  # noise: standard deviation phaseshift
  # len: length of simulations
  # return.wide: returns data in wide or long format
  
  require(data.table)
  # Create a matrix of shifted times 
  time_matrix <- matrix(seq(0, len-1, freq), nrow = len, ncol = n)
  shifts <- rnorm(n, 0, noise)
  shifts <- matrix(shifts, nrow = len, ncol = n, byrow = T)
  time_matrix <- time_matrix + shifts

  # Replace each shifted time by it sine function
  sins <- sin(time_matrix)
  
  # Go to data.table
  sins <- as.data.table(sins)
  sins <- cbind(seq(0, len-1, freq), sins)
  colnames(sins)[1] <- "Time"
  if(return.wide){
    return(sins)  
  }
  
  # Format long data table
  sins <- melt(sins, id.vars = "Time")
  return(sins)
}


sim_phase_shifted_with_fixed_trend <- function(n, noise, slope, freq = 1, len = 100){
  # Add a trend, i.e. a linear increase or decrease, to simulations
  # See sim_phase_shifted for arguments. Slope indicates the slope of the trend (change of mean value per unit of time)
  
  sins <- sim_phase_shifted(n, noise, freq,len, return.wide = F)
  trend_vec <- unique(sins$Time)
  trend_vec <- trend_vec * slope
  sins[, value := value + trend_vec, by = .(variable)]
  return(sins)
}


sim_noisy_amplitude <- function(n, noise, freq = 1, len = 100, return.wide = F){
  # Create a matrix of times and noise
  time_matrix <- matrix(seq(0, len-1, freq), nrow = len, ncol = n)
  noise_matrix <- replicate(n, rnorm(len, 0, noise))
  
  # Replace each shifted time by it sine function and add white noise
  sins <- sin(time_matrix)
  sins <- sins + noise_matrix
  
  # Go to data.table
  sins <- as.data.table(sins)
  sins <- cbind(seq(0, len-1, freq), sins)
  colnames(sins)[1] <- "Time"
  if(return.wide){
    return(sins)  
  }
  
  # Format long data table
  sins <- melt(sins, id.vars = "Time")
  return(sins)
}


multi_sims <- function(type, noises, ...){
  # /!\ not optimized, growing data table
  # Generate multiple simulations of the indicated type.
  # noises: numeric vector with noise value, DO NOT PASS 0, as it is already used to initialize the output
  # ps: phase shifted
  # pst: phase shifted with linear trend
  # na: noisy amplitude
  
  if(!(type %in% c("ps", "pst", "na"))){
    stop("type must be one of c('ps', 'pst', 'na')")
  }
  
  # Initialize data.table with no noise
  if(type == "ps"){multi_sim <- sim_phase_shifted(noise = 0, ...)}
  else if(type == "pst"){multi_sim <- sim_phase_shifted_with_fixed_trend(noise = 0, ...)}
  else if(type == "na"){multi_sim <- sim_noisy_amplitude(noise = 0, ...)}
  
  multi_sim$noise <- 0
  for(noise in noises){
    if(type == "ps"){temp <- sim_phase_shifted(noise = noise, ...)}
    else if(type == "pst"){temp <- sim_phase_shifted_with_fixed_trend(noise = noise, ...)}
    else if(type == "na"){temp <- sim_noisy_amplitude(noise = noise, ...)}
    temp$noise <- noise
    multi_sim <- rbind(multi_sim, temp)
  }
  return(multi_sim)
}


plot_sim <- function(data, x = "Time", y = "value", group = "variable", use.facet = T, facet = "noise", alpha = 0.2){
  require(ggplot2)
  p <- ggplot(data, aes_string(x = x, y = y)) + geom_line(aes_string(group = group), alpha = alpha)
  p <- p + stat_summary(fun.y=mean, geom="line", colour = "blue", size = 1.5)
  if(use.facet){
    p <- p + facet_grid(as.formula(paste(facet, "~ .")))
  }
  p
}



noises <- seq(0.2, 3, 0.2)
multi_sim <-  multi_sims(type="pst", noises=noises, n = 50, slope = 0.1)
plot_sim(multi_sim)

DistMean <- dist_mean(data = multi_sim, condition = "noise", tcol = "Time", measure = "value", label = "variable")
p <- ggplot(DistMean, aes(x = as.factor(noise), y = euclid_to_mean)) + geom_boxplot(aes(group = noise)) + ggtitle("Euclidian distance to mean trajectory"); p

Clip <- multi_sim[, .(clip = wrap_clip(value, k = 5)), by = .(noise, variable)]
# Convert noise (condition) and variable (label) to integers to comply with overlap constrains
Clip[, ':=' (noise = as.integer(noise*10),
             variable = as.integer(gsub("V", "", as.character(variable), fixed = T)))]
Overlap <- overlap_clipping(data = Clip, condition = "noise", label = "variable", measure = "clip")
q <- ggplot(Overlap, aes(x = as.factor(noise), y = Overlap)) + geom_boxplot(aes(group  = noise)) + scale_x_discrete(labels = as.character(c(0,noises))) + ggtitle("Pairwise overlap of clipped trajectories") ; q


pdf("Plots_Large_Sim.pdf", width = 10)
plot_sim(multi_sim)
p
q
dev.off()


multi_sim[, ':=' (noise = as.integer(noise*10),
             variable = as.integer(gsub("V", "", as.character(variable), fixed = T)))]
Correlations_Pearson <- correlations_group_label(multi_sim, "noise", "variable","value", method = "pearson")
Correlations_Spearman <- correlations_group_label(multi_sim, "noise", "variable","value", method = "spearman")
Correlations_Kendall <- correlations_group_label(multi_sim, "noise", "variable","value", method = "kendall")


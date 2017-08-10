sim_phase_shifted <- function(n, noise, freq = 1, end = 100, return.wide = F){
  # A fucntion to simulate a population of n noisy sinusoidals that are phase shifted
  # n: number of sinusoids
  # freq: sampling rate (in time unit)
  # noise: standard deviation phaseshift
  # end: end time of simulations
  # return.wide: returns data in wide or long format
  
  require(data.table)
  # If sampling rate < 1, adjust end such that seq(from, to, by) re
  
  # Create a matrix of shifted times 
  tvec <- seq(0, end-1, by = freq)
  time_matrix <- matrix(tvec, nrow = length(tvec), ncol = n)
  shifts <- rnorm(n, 0, noise)
  shifts <- matrix(shifts, nrow = length(tvec), ncol = n, byrow = T)
  time_matrix <- time_matrix + shifts

  # Replace each shifted time by it sine function
  sins <- sin(time_matrix)
  
  # Go to data.table
  sins <- as.data.table(sins)
  sins <- cbind(seq(0, end-1, by = freq), sins)
  colnames(sins)[1] <- "Time"
  if(return.wide){
    return(sins)  
  }
  
  # Format long data table
  sins <- melt(sins, id.vars = "Time")
  return(sins)
}


sim_phase_shifted_with_fixed_trend <- function(n, noise, slope, freq = 1, end = 100){
  # Add a trend, i.e. a linear increase or decrease, to simulations
  # See sim_phase_shifted for arguments. Slope indicates the slope of the trend (change of mean value per unit of time)
  
  sins <- sim_phase_shifted(n, noise, freq, end, return.wide = F)
  trend_vec <- unique(sins$Time)
  trend_vec <- trend_vec * slope
  sins[, value := value + trend_vec, by = .(variable)]
  return(sins)
}


sim_noisy_amplitude <- function(n, noise, freq = 1, end = 100, return.wide = F){
  # Create a matrix of times and noise
  tvec <- seq(0, end-1, by = freq)
  time_matrix <- matrix(tvec, nrow = length(tvec), ncol = n)
  noise_matrix <- replicate(n, rnorm(length(tvec), 0, noise))
  
  # Replace each shifted time by it sine function and add white noise
  sins <- sin(time_matrix)
  sins <- sins + noise_matrix
  
  # Go to data.table
  sins <- as.data.table(sins)
  sins <- cbind(seq(0, end-1, by = freq), sins)
  colnames(sins)[1] <- "Time"
  if(return.wide){
    return(sins)  
  }
  
  # Format long data table
  sins <- melt(sins, id.vars = "Time")
  return(sins)
}


sim_expodecay_lagged_stim <- function(n, noise, interval.stim = 5, lambda = 0.2, freq = 1, end = 100, return.wide = F){
  # A fucntion to simulate a population of n trajectories of expo decay, that are taken up at regular time by stimulus.
  # The noise add a lag between the stimulus and the increase.
  # n: number of trajectories
  # noise: standard deviation for the lag
  # interval.stim: how long between each stimulation?
  # lambda: desintegration rate
  # freq: sampling rate (in time unit), how often do we compute the state of the system
  # end: end time of simulations
  # return.wide: returns data in wide or long format
  
  # Time vector
  tvec <- seq(0, end-1, by = freq)
  # Matrix with stimulation times
  stim_time <- seq(interval.stim, end-1 , interval.stim)
  stim_time_matrix <- matrix(stim_time, nrow = length(stim_time), ncol = n)
  
  # Randomize the stimulation times (represent random lag), forbid lag < 0
  noise_matrix <- abs(replicate(n, rnorm(n = length(stim_time), mean = 0, sd = noise)))
  stim_time_matrix <- stim_time_matrix + noise_matrix
  
  # Initialize trajectories with 0 everywhere, set to 1 at stimulus times
  trajs <- matrix(0, nrow = length(tvec), ncol = n)
  for(col in 1:ncol(stim_time_matrix)){
    for(row in 1:nrow(stim_time_matrix)){
      index <- which(tvec >= stim_time_matrix[row, col])[1]
      trajs[index, col] <- 1
    }
  }
  
  # Expo decay computed thanks to previous value
  decrease_factor <- exp(-lambda * freq)
  for(col in 1:ncol(trajs)){
    for(row in 2:nrow(trajs)){
      # If not at a stim time, decay
      if(trajs[row, col] != 1){trajs[row, col] <- trajs[row-1, col] * decrease_factor}
    }
  }
  
  # Go to data.table
  trajs <- as.data.table(trajs)
  trajs <- cbind(seq(0, end-1, by = freq), trajs)
  colnames(trajs)[1] <- "Time"
  if(return.wide){
    return(trajs)  
  }
  
  # Format long data table
  trajs <- melt(trajs, id.vars = "Time")
  return(trajs)
}


multi_sims <- function(type, noises, ...){
  # /!\ not optimized, growing data table
  # Generate multiple simulations of the indicated type.
  # noises: numeric vector with noise value, DO NOT PASS 0, as it is already used to initialize the output
  # ps: phase shifted
  # pst: phase shifted with linear trend
  # na: noisy amplitude
  # edls: expo decay lagged stimulation
  
  if(!(type %in% c("ps", "pst", "na", "edls"))){
    stop("type must be one of c('ps', 'pst', 'na', 'edls')")
  }
  
  # Initialize data.table with no noise
  if(type == "ps"){multi_sim <- sim_phase_shifted(noise = 0, ...)}
  else if(type == "pst"){multi_sim <- sim_phase_shifted_with_fixed_trend(noise = 0, ...)}
  else if(type == "na"){multi_sim <- sim_noisy_amplitude(noise = 0, ...)}
  else if(type == "edls"){multi_sim <- sim_expodecay_lagged_stim(noise = 0, ...)}
  
  multi_sim$noise <- 0
  for(noise in noises){
    if(type == "ps"){temp <- sim_phase_shifted(noise = noise, ...)}
    else if(type == "pst"){temp <- sim_phase_shifted_with_fixed_trend(noise = noise, ...)}
    else if(type == "na"){temp <- sim_noisy_amplitude(noise = noise, ...)}
    else if(type == "edls"){temp <- sim_expodecay_lagged_stim(noise = noise, ...)}
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

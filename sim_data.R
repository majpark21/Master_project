sim_phase_shifted <- function(n, noise, freq = 1 ,len = 100, return.wide = F){
  # A fucntion to simulate a population of n noisy sinusoidals, at frequency freq 
  # n: number of sinusoids
  # freq: sampling rate (in time unit)
  # noise: standard deviation phaseshift
  # len: length of simulations
  # return.matrix: returns data in wide or long format
  
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


plot_sim <- function(data, alpha = 0.2){
  require(ggplot2)
  p <- ggplot(data, aes(x = Time, y = value)) + geom_line(aes(group = variable), alpha = alpha)
  p <- p + stat_summary(fun.y=mean, geom="line", colour = "blue", size = 1.5)
  p <- p + facet_grid(noise ~ .)
  p
}



# vector of noise?
noises <- seq(0.2, 1.2, 0.2)
n <- 50
len <- 100
multi_sim <- sim_phase_shifted(n, 0, len=len)
multi_sim$noise <- 0
for(noise in noises){
  temp <- sim_phase_shifted(n, noise, len=len)
  temp$noise <- noise
  multi_sim <- rbind(multi_sim, temp)
}

plot_sim(multi_sim)

DistMean <- dist_mean(data = multi_sim, condition = "noise", tcol = "Time", measure = "value", label = "variable")
p <- ggplot(DistMean, aes(x = as.factor(noise), y = euclid_to_mean)) + geom_boxplot(aes(group = noise)); p

Clip <- multi_sim[, .(clip = wrap_clip(value, k = 5)), by = .(noise, variable)]
# Convert noise (condition) and variable (label) to integers to comply with overlap constrains
Clip[, ':=' (noise = as.integer(noise*10),
             variable = as.integer(gsub("V", "", as.character(variable), fixed = T)))]
Overlap <- overlap_clipping(data = Clip, condition = "noise", label = "variable", measure = "clip")
q <- ggplot(Overlap, aes(x = as.factor(noise), y = Overlap)) + geom_boxplot(aes(group  = noise)) + scale_x_discrete(labels = as.character(c(0,noises))) ; q





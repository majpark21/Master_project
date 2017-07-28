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
  p
}
plot_sim(sim2)


# TODO: Implement condition as noise, simulate for a vector of noise?
sim <- sim_phase_shifted(n = 50, noise = 0.1)
sim$condition <- 1
a <- dist_mean(sim, "condition", "Time", "value", "variable")


sim2 <- sim_phase_shifted(n = 50, noise = 1)
sim2$condition <- 1
b <- dist_mean(sim2, "condition", "Time", "value", "variable")

par(mfrow=c(1,2))
boxplot(a$euclid_to_mean, ylim = c(0,10), main="noise 1")
boxplot(b$euclid_to_mean, ylim = c(0,10), main="noise X")

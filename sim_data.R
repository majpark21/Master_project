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


temp <- sim_phase_shifted(5, 1, 1, 50, return.wide = T)


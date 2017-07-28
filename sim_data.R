sim_phase_shifted <- function(n, freq, noise, len = 100){
  # A fucntion to simulate a population of n noisy sinusoidals, at frequency freq 
  # n: number of sinusoids
  # freq: oscillation freq (in time-1)
  # noise: standard deviation phaseshift
  # len: length of simulations
  
  # Create a matrix of shifted times 
  time_matrix <- matrix(seq(0, len-1), nrow = len, ncol = n)
  shifts <- rnorm(n, 0, noise)
  shifts <- matrix(shifts, nrow = len, ncol = n, byrow = T)
  time_matrix <- time_matrix + shifts

  # Replace each shifted time by it sine function
  sins <- sin(time_matrix)
  return(sins)
}

temp <- sim_phase_shifted(5, 1, 1, 50)

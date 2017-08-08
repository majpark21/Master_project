amplitude_oscillations <- function(data, condition, label, measure, k_roll_mean = 5){  
  # Compute the distance for each trajectory to its rolling mean
  
  # data: a data table containing trajectories
  # condition: column name that fully define an experimental condition; MUST BE INTEGERS, NUMERIC, FACTOR OR CHARACTERS
  # measure: column name with clipped trajectories
  # label: column name with label of individual objects in each condition (cell label); LABELS MUST BE INTEGERS, NUMERIC, FACTOR OR CHARACTERS
  
  out <- data[, .(euclid_to_roll_mean = sqrt( sum( (get(measure) - rollex(get(measure), k = k_roll_mean) )^2 )) ), by = c(condition, label)]
  return(out)
}

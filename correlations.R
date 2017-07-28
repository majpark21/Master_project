correlations_group_label <- function(data, condition, label, measure, method){
  # Compute overlap between all pairs of clipped trajectories in each condition
  
  # data: a data table containing CLIPPED trajectories
  # condition: column name that fully define an experimental condition; MUST BE INTEGERS
  # measure: column name with clipped trajectories
  # label: column name with label of individual objects in each condition (cell label); LABELS MUST BE INTEGERS
  # method: method for cor
  
  require(data.table)
  setkeyv(data, c(condition, label))
  
  # Number of rows for data table initialization, sum of number of pairs in each condition
  nber_row <- 0
  for(i in unique(data[, get(condition)])){
    nber_row <- nber_row + choose(length(unique(data[.(i), get(label)])), 2)
  }
  
  # One row = one pair in one condition; set column types
  out <- data.table(matrix(ncol = 4, nrow = nber_row))
  colnames(out) <- c(condition, "Label1", "Label2", "Correlations")
  out <- out[, lapply(.SD, as.integer)]
  out[, Overlap := as.numeric(Overlap)]
  
  curr_row <- 1L
  # Loop condition
  for(i in unique(data[, get(condition)])){
    labels <- unique(data[.(i), get(label)])
    # Loop 1st label
    for(j in 1:(length(labels)-1)){
      # Loop 2nd label
      for(k in (j+1):length(labels)){
        set(out, curr_row, 1L, i)
        set(out, curr_row, 2L, labels[j])
        set(out, curr_row, 3L, labels[k])
        set(out, curr_row, 4L, cor(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)]), method=method)
        curr_row <- curr_row + 1L
      }
    }
  }
  return(out)
}



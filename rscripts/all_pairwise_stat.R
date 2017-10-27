all_pairwise_stats <- function(data, condition, label, measure, k_roll_mean = 5){
  # Compute all pairwise measures: overlap_clipping and correlations (Pearson, Spearman, Kendall)
  # The clipping of trajectories is performed as first step
  
  # data: a data table containing CLIPPED trajectories
  # condition: column name that fully define an experimental condition; MUST BE INTEGERS, NUMERIC, FACTOR OR CHARACTERS
  # measure: column name with clipped trajectories
  # label: column name with label of individual objects in each condition (cell label); LABELS MUST BE INTEGERS, NUMERIC, FACTOR OR CHARACTERS
  
  if(!(class(data[,get(condition)]) %in% c("integer", "numeric", "character", "factor"))){
    stop("Column 'condition' must be integer, numeric, factor or character.")
  }
  if(!(class(data[,get(label)]) %in% c("integer", "numeric", "character", "factor"))){
    stop("Column 'label' must be integer, numeric, factor or character.")
  }
  
  require(data.table)
  setkeyv(data, c(condition, label))
  
  # Perform clipping
  ClipData <- data[, .(clip_measure = wrap_clip(get(measure), k = k_roll_mean)), by = c(condition, label)]
  setkeyv(ClipData, c(condition, label))
  
  # Number of rows for data table initialization, sum of number of pairs in each condition
  nber_row <- 0
  for(i in unique(data[, get(condition)])){
    nber_row <- nber_row + choose(length(unique(data[.(i), get(label)])), 2)
  }

  
  # One row = one pair in one condition; set column types according to input types
  out <- data.table(matrix(ncol = 7, nrow = nber_row))
  colnames(out) <- c(condition, "Label1", "Label2", "Overlap", "Pearson", "Spearman", "Kendall")
  if(class(data[,get(condition)]) == "integer"){out[[condition]] <- as.integer(out[[condition]])}
  else if(class(data[,get(condition)]) == "numeric"){out[[condition]] <- as.numeric(out[[condition]])}
  else if(class(data[,get(condition)]) == "character"){out[[condition]] <- as.character(out[[condition]])}
  else if(class(data[,get(condition)]) == "factor"){out[[condition]] <- as.factor(out[[condition]])}
  
  
  if(class(data[,get(label)]) == "integer"){
    out[, Label1 := as.integer(Label1)]
    out[, Label2 := as.integer(Label2)]
  }
  else if(class(data[,get(label)]) == "numeric"){
    out[, Label1 := as.numeric(Label1)]
    out[, Label2 := as.numeric(Label2)]
  }
  else if(class(data[,get(label)]) == "character"){
    out[, Label1 := as.character(Label1)]
    out[, Label2 := as.character(Label2)]
  }
  else if(class(data[,get(label)]) == "factor"){
    out[, Label1 := as.factor(Label1)]
    out[, Label2 := as.factor(Label2)]
  }
  
  
  out[, Overlap := as.numeric(Overlap)]
  out[, Pearson := as.numeric(Pearson)]
  out[, Spearman := as.numeric(Spearman)]
  out[, Kendall := as.numeric(Kendall)]
  
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
        # /!\ use clip data for overlap
        set(out, curr_row, 4L, overlap(ClipData[.(i, labels[j]), clip_measure], ClipData[.(i, labels[k]), clip_measure]))
        set(out, curr_row, 5L, cor(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)], method ="pearson"))
        set(out, curr_row, 6L, cor(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)], method ="spearman"))
        set(out, curr_row, 7L, cor(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)], method ="kendall"))
        curr_row <- curr_row + 1L
      }
    }
  }
  return(out)
}


complete.time.series <- function(data, cond.col, lab.col, time.col, time.vector, meas.col, fill = NA){
  # Add a row with fill measurement for all missing measurement in the long format
  # time vector: over which time should ALL time series span? Missing times will be added to series where it's not present
  require(data.table)
  # 1) Perform left outer join with a table that contains all RealTime
  temp <- CJ(Condition=unique(data[[cond.col]]), Label=unique(data[[lab.col]]), RealTime=time.vector)
  temp <- merge(temp, data, by = c(cond.col, lab.col, time.col), all.x = T)
  # 2) Trim the rows that contains only NA (i.e. this combination of Label and Condition does not exist in the real data)
  out <- temp[, if(!all(is.na(get(meas.col)))) .SD, by = .(Condition, Label)]
  return(out)
}

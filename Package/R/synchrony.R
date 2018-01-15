all_pairwise_stats <- function(data, condition, label, measure, k_roll_mean = 5){
  # Compute all pairwise measures: DTW, overlap_clipping and correlations (Pearson, Spearman, Kendall)
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
  require(dtw)
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
  out <- data.table(matrix(ncol = 8, nrow = nber_row))
  colnames(out) <- c(condition, "Label1", "Label2", "Overlap", "Pearson", "Spearman", "Kendall", "DTW")
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
  out[, DTW := as.numeric(DTW)]

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
        set(out, curr_row, 8L, dtw(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)])$distance)
        curr_row <- curr_row + 1L
      }
    }
  }
  return(out)
}


# # Solution to try
# library(data.table)
# dat <- data.table(Region = rep(LETTERS[24:26], each=3),
#                   Name = rep(LETTERS[1:3], 3),
#                   Val1 = rep(rnorm(3), 3),
#                   Val2 = rep(rnorm(3), 3))
# dat2 <- merge(dat, dat, by="Region", allow.cartesian = T)[Name.x < Name.y]
# dat2[, Val1Ratio := Val1.x / Val1.y]
# dat2[, Val2Ratio := Val2.x / Val2.y]
# dat2[, Diff := Val1Ratio - Val2Ratio]

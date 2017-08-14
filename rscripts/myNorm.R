
# Returns original dt with an additional column with normalized quantity.
# The column to be normalised is given by 'in.meas.col'.
# The name of additional column is the same as in.meas.col but with ".norm" suffix added.
# Normalisation is based on part of the trajectory;
# this is defined by in.rt.min and max, and the column with time in.rt.col.
# Additional parameters:
# in.by.cols - character vector with 'by' columns to calculate normalisation per group
#              if NULL, no grouping is done
# in.robust - whether robust measures should be used (median instead of mean, mad instead of sd)
# in.type - type of normalization: z.score or mean (fi.e. old change w.r.t. mean)

myNorm = function(in.dt,
                  in.meas.col,
                  in.rt.col = 'RealTime',
                  in.rt.min = 10,
                  in.rt.max = 20,
                  in.by.cols = NULL,
                  in.robust = TRUE,
                  in.type = 'z.score') {
  loc.dt <-
    copy(in.dt) # copy so as not to alter original dt object w intermediate assignments
  
  if (is.null(in.by.cols)) {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = mad(get(in.meas.col), na.rm = TRUE))]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = sd(get(in.meas.col), na.rm = TRUE))]
    
    loc.dt = cbind(loc.dt, loc.dt.pre.aggr)
  }  else {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = mad(get(in.meas.col), na.rm = TRUE)), by = in.by.cols]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = sd(get(in.meas.col), na.rm = TRUE)), by = in.by.cols]
    
    loc.dt = merge(loc.dt, loc.dt.pre.aggr, by = in.by.cols)
  }
  
  
  if (in.type == 'z.score') {
    loc.dt[, meas.norm := (get(in.meas.col) - meas.md) / meas.mad]
  } else {
    loc.dt[, meas.norm := (get(in.meas.col) / meas.md)]
  }
  
  setnames(loc.dt, 'meas.norm', paste0(in.meas.col, '.norm'))
  
  loc.dt[, c('meas.md', 'meas.mad') := NULL]
  return(loc.dt)
}
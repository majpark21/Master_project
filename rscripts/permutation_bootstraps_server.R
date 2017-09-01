library(stringr)
library(data.table)
library(plyr)

CastCluster <- function(data, time.col, condition.col, label.col, measure.col, k.clust, na.fill, plot = T, return.quality = T, ...){
  # Cast to wide, cluster and get quality indices
  require(dtwclust)
  temp <- myCast(data, time.col, condition.col, label.col, measure.col, na.fill)
  
  # Make clustering, and get quality indexes
  clust <- tsclust(temp$casted.matrix, type = "partitional", k = k.clust, distance = "dtw_basic", centroid = "pam", seed = 42, ...)
  names(clust) <- paste0("k_", k.clust)
  quality <- sapply(clust, cvi, type = "internal")
  
  # Add a column with the clusters to the casted table
  cluster.table <- temp$casted
  for(k in 1:length(clust)){
    cluster.table <- cbind(clust[[k]]@cluster, cluster.table)
    colnames(cluster.table)[1] <- names(clust)[k]
  }
  
  # Plot
  if(plot){
    require(ggplot2)
    mquality <- melt(quality)
    names(mquality) <- c("Stat", "Nb.Clust", "value")
    plot(ggplot(mquality, aes(x=Nb.Clust, y = value)) + geom_col(aes(group = Nb.Clust, fill = Nb.Clust)) + facet_grid(Stat ~ ., scales = "free_y"))
  }
  
  # Output
  if(return.quality) return(list(cluster = clust, table = cluster.table, quality = quality))
  else return(list(out = clust, table = cluster.table))
}


myCast <- function(data, time.col, condition.col, label.col, measure.col, na.fill){
  # Only cast to wide matrix
  temp <- copy(data)
  # dcast can change the order of the rows depending on the orde rin which the keyed columns are passed, keep the casted table in addition to the matrix to make the link afterwards
  temp <- dcast(temp, get(condition.col) + get(label.col) ~ get(time.col), value.var = measure.col)
  temp2 <- as.matrix(temp[, c(-1, -2)]) # remove 2 first columns with labels
  temp2[which(is.na(temp2))] <- na.fill
  return(list(casted = temp, casted.matrix = temp2))
}


plot_cluster <- function(data, id.vars.col, cluster.col, type){
  # Plot clusters directly from output$table of CastCluster
  # id.vars.col: given in indices (include ALL clustering columns)
  # cluster.col: name of the column with clustering to plot
  library(ggplot2)
  ids <- colnames(data)[id.vars.col]
  melted <- melt(data, id.vars = ids)
  if(type=="trajectory"){
    ggplot(melted, aes(x = as.numeric(variable), y = value)) + geom_line(aes(group = label.col)) +
      facet_wrap(as.formula(paste("~",cluster.col))) + stat_summary(fun.y=mean, geom="line", colour = "blue", size = 1.5) + xlab("Time")
  } else if(type=="composition"){
    melted[, c(cluster.col):=as.factor(get(cluster.col))]
    ggplot(melted, aes_string(cluster.col)) + geom_bar(aes(fill=condition.col))
  }
}

sep.meas.along.time <- function(data1, data2, time.col, measure.col){
  timev <- unique(data1[, get(time.col)])
  if(!(identical(unique(data2[, get(time.col)]), timev))) stop("Time vectors must be identical between the two data")
  out <- separability.measures(data1[get(time.col)==timev[1], get(measure.col)], data2[get(time.col)==timev[1], get(measure.col)])
  for(t in timev[2:length(timev)]){
    out <- rbind(out, separability.measures(data1[RealTime==t, get(measure.col)], data2[RealTime==t, get(measure.col)]))
  }
  out <- cbind(timev, out)
  return(out)
}


one.permutation.auc <- function(x, y, metric){
  n <- nrow(x)
  m <- nrow(y)
  temp <- rbind(x, y)
  samp.traj <- sample(1:nrow(temp), size = n, replace = FALSE)
  x.resamp <- temp[samp.traj, ]
  y.resamp <- temp[setdiff(1:nrow(temp), samp.traj), ]
  
  seps <- sapply(1:ncol(x), function(j) separability.measures(x.resamp[, j], y.resamp[, j]))
  return(sum(unlist(seps[metric, ])))
}

permutation.auc <- function(x, y, n, metric = "jm"){
  # x,y: two matrices representing time series, row: trajectory; col: time
  # n: number of permutations
  # metric: one of "jm", "bh", "div", "tdiv", "ks"
  if(ncol(x) != ncol(y)) stop("x and y must have same number of columns")
  return(replicate(n, one.permutation.auc(x,y,metric)))
}


wrap_perm <- function(x, y, meas, n){
  a <- myCast(x, "RealTime", "Condition", "Label", meas, 1100)$casted.matrix
  b <- myCast(y, "RealTime", "Condition", "Label", meas, 1100)$casted.matrix
  return(permutation.auc(a, b, n))
}

one.bootstrap.auc.percol <- function(x, y, metric){
  samp.col <- sample(1:ncol(x), size = ncol(x), replace = TRUE)
  x.resamp <- x[, samp.col]
  y.resamp <- y[, samp.col]
  seps <- sapply(1:ncol(x), function(j) separability.measures(x.resamp[, j], y.resamp[, j]))
  return(sum(unlist(seps[metric, ])))
}

bootstrap.auc.percol <- function(x, y, B, metric = "jm"){
  # x,y: two matrices representing time series, row: trajectory; col: time
  # B: number of boostraps
  # metric: one of "jm", "bh", "div", "tdiv", "ks"
  if(ncol(x) != ncol(y)) stop("x and y must have same number of columns")
  return(replicate(B, one.bootstrap.auc.percol(x,y,metric)))
}

wrap_bootcol <- function(x, y, meas, n){
  a <- myCast(x, "RealTime", "Condition", "Label", meas, 1100)$casted.matrix
  b <- myCast(y, "RealTime", "Condition", "Label", meas, 1100)$casted.matrix
  return(bootstrap.auc.percol(a, b, n))
}


# Load data
Yanni <- fread("../input/sust_E_F_N.csv")
gf <- str_extract(Yanni$Stim_All_Ch, "[E,F,N]")
conc <- str_extract(Yanni$Stim_All_Ch, "(([0-9]+\\.[0-9]*)|([0-9]+))")
Yanni$Condition <- paste(gf, conc, sep = "-")
Yanni[, c("TrackObjects_Label_uni","Condition") := list(as.factor(TrackObjects_Label_uni), as.factor(Condition))]
setkey(Yanni, Condition, TrackObjects_Label_uni)
setnames(Yanni, c("Intensity_MeanIntensity_Ratio", "TrackObjects_Label_uni"), c("Ratio","Label"))
rm(gf, conc)
setcolorder(Yanni, c("Condition", "Label", "Ratio", "RealTime", "Metadata_Series", "Metadata_T", "TrackObjects_Label","Stim_All_Ch", "Stim_All_S"))
del.cols <- names(Yanni)[5:9]
Yanni[, (del.cols) := NULL]

set(Yanni, i = which(Yanni[, Ratio] > 1800), j=3L, value = 1100) # Ratio is the 3rd column


# Compute AUC
# Get all pairs of conditions
conditions <- combn(as.character(unique(Yanni[,Condition])), m = 2)
conditions <- conditions[,c(1:3,12,13,22, 39:41,46,47,52, 61:66)]
max.val <- sqrt(2) * 101
# Compute separabilities of conditions at each time point
sep.meas.raw <- apply(conditions, 2, function(x) sep.meas.along.time(Yanni[Condition==x[1]],  Yanni[Condition==x[2]], "RealTime", "Ratio" ))
names(sep.meas.raw) <- apply(conditions, 2, function(x) paste(x[1], x[2], sep = ","))
# Go to data table
for(i in 1:length(sep.meas.raw)){
  temp <- unlist(strsplit(names(sep.meas.raw)[i], ","))
  sep.meas.raw[[i]]$Cond1 <- temp[1]
  sep.meas.raw[[i]]$Cond2 <- temp[2]
}
sep.meas.raw <- as.data.table(rbind.fill(sep.meas.raw))
#sep.meas.raw[, c("Cond1", "Cond2") := list(as.factor(Cond1), as.factor(Cond2))]
auc.raw <- sep.meas.raw[, .(auc = sum(jm, na.rm = T)/max.val), by = c("Cond1", "Cond2")] # a few NAs, slight bias in the values
auc.raw[, comb.cond := as.factor(paste(Cond1, Cond2, sep = ";"))]
auc.raw[, GF := as.factor(paste0(str_sub(Cond1,1,1), "GF"))]


# Distance matrix
get.dist.matrix <- function(dat){
  data <- copy(dat)
  data[, Cond1 := factor(data$Cond1, levels = union(data$Cond1, data$Cond2))]
  data[, Cond2 := factor(data$Cond2, levels = union(data$Cond1, data$Cond2))]
  dist.mat <- dcast(data=data[,1:3], formula = Cond1 ~ Cond2, value.var = "auc", drop = F)
  dist.mat <- as.matrix(dist.mat[,-1])
  rownames(dist.mat) <- colnames(dist.mat)
  dist.mat[which(is.na(dist.mat))] <- 0
  dist.mat <- dist.mat + t(dist.mat)
  return(dist.mat)
}
dist.raw.EGF <- get.dist.matrix(auc.raw[GF=="EGF"])
dist.raw.FGF <- get.dist.matrix(auc.raw[GF=="FGF"])
dist.raw.NGF <- get.dist.matrix(auc.raw[GF=="NGF"])


# Permutation tests and bootstraps
nperm <- 25
set.seed(7)
auc.perm <- apply(conditions, 2, function(x) wrap_perm(Yanni[Condition==x[1]], Yanni[Condition==x[2]], "Ratio", nperm))
bootcol <- apply(conditions, 2, function(x) wrap_bootcol(Yanni[Condition==x[1]], Yanni[Condition==x[2]], "Ratio", nperm))

save(auc.perm, file = "../output_perm_boot/auc.perm.Robj")
save(bootcol, file = "../output_perm_boot/bootcol.Robj")
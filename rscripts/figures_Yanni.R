# Load data ----
source("rscripts/package.R")
del.cols <- names(Yanni)[5:9]
Yanni[, (del.cols) := NULL]
Yanni <- myNorm(in.dt = Yanni, in.meas.col = "Ratio", in.rt.min = 0, in.rt.max = 36, in.by.cols = c("Condition", "Label"), in.type = "fold.change")
Yanni <- Yanni[RealTime <= 200]
Yanni[, Ratio.norm.smooth := rollex(Ratio.norm, k = 5), by = c("Condition", "Label")]

# Clip the big outliers
set(Yanni, i = which(Yanni[, Ratio] > 1800), j=3L, value = 1100) # Ratio is the 3rd column
set(Yanni, i = which(Yanni[, Ratio.norm] > 1.4), j=5L, value = 1.1) # Ratio.norm is the 5th column
set(Yanni, i = which(Yanni[, Ratio.norm.smooth] > 1.4), j=6L, value = 1.1) # Ratio.norm.smooth is the 6th column

# Functions definitions ----

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

# PCA ----

library(ggbiplot)
cond <- "Condition"
lab <-  "Label"
tim <- "RealTime"

EGF <- Yanni[Condition %in% c("E-0.25", "E-2.5", "E-25", "E-250")]
NGF <- Yanni[Condition %in% c("N-0.25", "N-2.5", "N-25", "N-250")]
FGF <- Yanni[Condition %in% c("F-0.25", "F-2.5", "F-25", "F-250")]

# Ignore the na.fill argument
cast_EGF <- myCast(EGF, tim, cond, lab, "Ratio.norm.smooth", na.fill = 1.1)
pca_EGF <- prcomp(cast_EGF$casted.matrix, center = T)

pca_EGF_plot <- ggbiplot(pca_EGF, obs.scale = 1, var.scale = 1, 
              groups = unlist(cast_EGF$casted[,1]), ellipse = TRUE, 
              circle = F, var.axes = F)
pca_EGF_plot <- pca_EGF_plot + scale_color_discrete(name = '')
pca_EGF_plot <- pca_EGF_plot + theme(legend.direction = 'horizontal', 
                                     legend.position = 'top') +
                ggtitle("PCA EGF - Raw Ratio - Ellipse 0.68")
print(pca_EGF_plot)

cast_NGF <- myCast(NGF, tim, cond, lab, "Ratio.norm.smooth", na.fill = 1.1)
pca_NGF <- prcomp(cast_NGF$casted.matrix, center = T)

pca_NGF_plot <- ggbiplot(pca_NGF, obs.scale = 1, var.scale = 1, 
                         groups = unlist(cast_NGF$casted[,1]), ellipse = TRUE, 
                         circle = F, var.axes = F)
pca_NGF_plot <- pca_NGF_plot + scale_color_discrete(name = '')
pca_NGF_plot <- pca_NGF_plot + theme(legend.direction = 'horizontal', 
                                     legend.position = 'top') + 
                ggtitle("PCA NGF - Raw Ratio - Ellipse 0.68")
print(pca_NGF_plot)

cast_FGF <- myCast(FGF, tim, cond, lab, "Ratio.norm.smooth", na.fill = 1.1)
pca_FGF <- prcomp(cast_FGF$casted.matrix, center = T)

pca_FGF_plot <- ggbiplot(pca_FGF, obs.scale = 1, var.scale = 1, 
                         groups = unlist(cast_FGF$casted[,1]), ellipse = TRUE, 
                         circle = F, var.axes = F)
pca_FGF_plot <- pca_FGF_plot + scale_color_discrete(name = '')
pca_FGF_plot <- pca_FGF_plot + theme(legend.direction = 'horizontal', 
                                     legend.position = 'top') +
                ggtitle("PCA FGF - Raw Ratio - Ellipse 0.68")
print(pca_FGF_plot)

# AUC heatmaps ----

library(plyr)
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
sep.meas.raw[, c("Cond1", "Cond2") := list(as.factor(Cond1), as.factor(Cond2))]
auc.raw <- sep.meas.raw[, .(auc = sum(jm, na.rm = T)/max.val), by = c("Cond1", "Cond2")] # a few NAs, slight bias in the values
auc.raw[, comb.cond := as.factor(paste(Cond1, Cond2, sep = ";"))]
auc.raw[, GF := as.factor(paste0(str_sub(Cond1,1,1), "GF"))]

# Split the GF prior to?
auc.raw[, Cond1 := factor(auc.raw$Cond1, levels = union(auc.raw$Cond1, auc.raw$Cond2))]
auc.raw[, Cond2 := factor(auc.raw$Cond2, levels = union(auc.raw$Cond1, auc.raw$Cond2))]
dist.raw <- dcast(data=auc.raw[,1:3], formula = Cond1 ~ Cond2, value.var = "auc", drop = F)

# # Compute separabilities of conditions at each time point
# sep.meas.norm <- apply(conditions, 2, function(x) sep.meas.along.time(Yanni[Condition==x[1]],  Yanni[Condition==x[2]], "RealTime", "Ratio.norm" ))
# names(sep.meas.norm) <- apply(conditions, 2, function(x) paste(x[1], x[2], sep = ","))
# # Go to data table
# for(i in 1:length(sep.meas.norm)){
#   temp <- unlist(strsplit(names(sep.meas.norm)[i], ","))
#   sep.meas.norm[[i]]$Cond1 <- temp[1]
#   sep.meas.norm[[i]]$Cond2 <- temp[2]
# }
# sep.meas.norm <- as.data.table(rbind.fill(sep.meas.norm))
# sep.meas.norm[, c("Cond1", "Cond2") := list(as.factor(Cond1), as.factor(Cond2))]
# auc.norm <- sep.meas.norm[, .(auc = sum(jm, na.rm = T)/max.val), by = c("Cond1", "Cond2")] # a few NAs, slight bias in the values
# 
# # Compute separabilities of conditions at each time point
# sep.meas.smooth <- apply(conditions, 2, function(x) sep.meas.along.time(Yanni[Condition==x[1]],  Yanni[Condition==x[2]], "RealTime", "Ratio.norm.smooth" ))
# names(sep.meas.smooth) <- apply(conditions, 2, function(x) paste(x[1], x[2], sep = ","))
# # Go to data table
# for(i in 1:length(sep.meas.smooth)){
#   temp <- unlist(strsplit(names(sep.meas.smooth)[i], ","))
#   sep.meas.smooth[[i]]$Cond1 <- temp[1]
#   sep.meas.smooth[[i]]$Cond2 <- temp[2]
# }
# sep.meas.smooth <- as.data.table(rbind.fill(sep.meas.smooth))
# sep.meas.smooth[, c("Cond1", "Cond2") := list(as.factor(Cond1), as.factor(Cond2))]
# auc.smooth <- sep.meas.smooth[, .(auc = sum(jm, na.rm = T)/max.val), by = c("Cond1", "Cond2")] # a few NAs, slight bias in the values

ggplot(auc.raw, aes(x=Cond1, y=Cond2)) + 
  geom_raster(aes(fill=auc)) +
  facet_wrap(~GF, scales = "free")


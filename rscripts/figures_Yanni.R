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
    ggplot(melted, aes_string(cluster.col)) + geom_bar(aes(fill=condition.col), position = "fill")
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
max.val <- 2 * 101

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



# 1) Distance matrix

# Split the GF prior to plot

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


# 2) Add p-values coming from permutation tests (load data from 2000 iterations)
load("data/Yannick_sustained/auc.perm2000.Robj")
# nperm <- 25
# set.seed(7)
# auc.perm <- apply(conditions, 2, function(x) wrap_perm(Yanni[Condition==x[1]], Yanni[Condition==x[2]], "Ratio", nperm))
# auc.perm <- auc.perm/max.val
# colnames(auc.perm) <- apply(conditions, 2, paste, collapse=";")

p.values.perm <- auc.raw[, .(Cond1, Cond2, auc)]
p.values.perm$p.val <- NA
for(i in 1:nrow(p.values.perm)){
  p.values.perm$p.val[i] <- sum(p.values.perm$auc[i] <= auc.perm[,i]) / 2000
}
# Go to symmetrical data
p.values.perm[, GF := as.factor(paste0(str_sub(Cond1,1,1), "GF"))]

get.pval.matrix <- function(dat){
  data <- copy(dat)
  data[, Cond1 := factor(data$Cond1, levels = union(data$Cond1, data$Cond2))]
  data[, Cond2 := factor(data$Cond2, levels = union(data$Cond1, data$Cond2))]
  data <- dcast(data=data[,c(1,2,4)], formula = Cond1 ~ Cond2, value.var = "p.val", drop = F)
  data <- as.matrix(data[,-1])
  rownames(data) <- colnames(data)
  data[which(is.na(data))] <- 0
  data <- data + t(data)
  return(data)
}

# Correction for multiple tests (Family-wise, no false positive correction)
p.values.perm$p.val <- p.adjust(p = p.values.perm$p.val, method = "holm")
p.values.EGF <- get.pval.matrix(p.values.perm[GF=="EGF"])
p.values.FGF <- get.pval.matrix(p.values.perm[GF=="FGF"])
p.values.NGF <- get.pval.matrix(p.values.perm[GF=="NGF"])


# 3) Add confidence intervals coming from bootstraps (load data from 2000 iterations)
load("data/Yannick_sustained/bootcol2000.Robj")
# nperm <- 25
# set.seed(7)
# bootcol <- apply(conditions, 2, function(x) wrap_bootcol(Yanni[Condition==x[1]], Yanni[Condition==x[2]], "Ratio", nperm))
# bootcol <- bootcol/max.val
# colnames(bootcol) <- apply(conditions, 2, paste, collapse=";")
ci.bootcol <- auc.raw[, .(Cond1, Cond2, auc)]
ci.bootcol$low <- NA
ci.bootcol$high <- NA
for(i in 1:nrow(ci.bootcol)){
  ci.bootcol$low[i] <- quantile(bootcol[,i], 0.025)
  ci.bootcol$high[i] <- quantile(bootcol[,i], 0.975)
}
ci.bootcol[, GF := as.factor(paste0(str_sub(Cond1,1,1), "GF"))]
ci.bootcol[, CI := paste0("[ " , as.character(round(low,3)), ",\n   ", as.character(round(high,3)), " ]")]

get.ci.matrix <- function(dat){
  data <- copy(dat)
  data[, Cond1 := factor(data$Cond1, levels = union(data$Cond1, data$Cond2))]
  data[, Cond2 := factor(data$Cond2, levels = union(data$Cond1, data$Cond2))]
  data <- dcast(data=data[,c(1,2,7)], formula = Cond1 ~ Cond2, value.var = "CI", drop = F)
  data <- as.matrix(data[,-1])
  rownames(data) <- colnames(data)
  data[which(is.na(data))] <- ""
  data <- matrix(paste0(data, t(data)), nrow = 4)
  return(data)
}


ci.bootcol.EGF <- get.ci.matrix(ci.bootcol[GF=="EGF"])
ci.bootcol.FGF <- get.ci.matrix(ci.bootcol[GF=="FGF"])
ci.bootcol.NGF <- get.ci.matrix(ci.bootcol[GF=="NGF"])


# 4) Create annotation matrices
# Create matrices where the lower triangle is CI and upper is pval
get.pval.char.matrix <- function(mat){
  out <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      if(mat[i,j] <= 1e-4) out[i,j] <- "***"
      if(mat[i,j] <= 1e-3 & mat[i,j] >= 1e-4) out[i,j] <- "**"
      if(mat[i,j] <= 1e-2 & mat[i,j] >= 1e-3) out[i,j] <- "*"
      if(mat[i,j] <= 5e-2 & mat[i,j] >= 1e-2) out[i,j] <- "."
      if(mat[i,j] >= 5e-2) out[i,j] <- "0"
    }
  }
  return(out)
}

p.values.EGF.char <- get.pval.char.matrix(p.values.EGF)
p.values.FGF.char <- get.pval.char.matrix(p.values.FGF)
p.values.NGF.char <- get.pval.char.matrix(p.values.NGF)

annot.matrix.EGF <- matrix("", nrow=4, ncol=4)
annot.matrix.FGF <- matrix("", nrow=4, ncol=4)
annot.matrix.NGF <- matrix("", nrow=4, ncol=4)

annot.matrix.EGF[lower.tri(annot.matrix.EGF)] <- ci.bootcol.EGF[lower.tri(ci.bootcol.EGF)]
annot.matrix.EGF[upper.tri(annot.matrix.EGF)] <- p.values.EGF.char[upper.tri(ci.bootcol.EGF)]
annot.matrix.FGF[lower.tri(annot.matrix.FGF)] <- ci.bootcol.FGF[lower.tri(ci.bootcol.FGF)]
annot.matrix.FGF[upper.tri(annot.matrix.FGF)] <- p.values.FGF.char[upper.tri(ci.bootcol.FGF)]
annot.matrix.NGF[lower.tri(annot.matrix.NGF)] <- ci.bootcol.NGF[lower.tri(ci.bootcol.NGF)]
annot.matrix.NGF[upper.tri(annot.matrix.NGF)] <- p.values.NGF.char[upper.tri(ci.bootcol.NGF)]

# 5) Plot the heatmaps and dendrograms
library(RColorBrewer)
library(gplots)
mypal <- brewer.pal(9, "Reds")
mypal <- colorspace::diverge_hsv(15, v=1, s=0.75)
mybreaks <- quantile(auc.raw$auc, probs = seq(0,1, length.out = 10))
mybreaks <- quantile(auc.raw$auc, probs = seq(0,1, length.out = 16))
color.annot <- "black"

heatmap.2(dist.raw.EGF, Rowv = F, Colv = F, dendrogram = "none", col=mypal, breaks = mybreaks, key = T, trace= "none", main = "AUC separability - EGF", margins=c(7,7), cellnote = annot.matrix.EGF, notecex = 1.75, notecol = color.annot)
heatmap.2(dist.raw.FGF, Rowv = F, Colv = F, dendrogram = "none", col=mypal, breaks = mybreaks, key = T, trace= "none", main = "AUC separability - FGF", margins=c(7,7), cellnote = annot.matrix.FGF, notecex = 1.75, notecol = color.annot)
heatmap.2(dist.raw.NGF, Rowv = F, Colv = F, dendrogram = "none", col=mypal, breaks = mybreaks, key = T, trace= "none", main = "AUC separability - NGF", margins=c(7,7), cellnote = annot.matrix.NGF, notecex = 1.75, notecol = color.annot)

plot(as.dendrogram(hclust(as.dist(dist.raw.EGF))), ylim = c(0, 0.225), main="AUC as distance matrix - EGF")
plot(as.dendrogram(hclust(as.dist(dist.raw.FGF))), ylim = c(0, 0.225), main="AUC as distance matrix - FGF")
plot(as.dendrogram(hclust(as.dist(dist.raw.NGF))), ylim = c(0, 0.225), main="AUC as distance matrix - NGF")


# Pamk plot ----

EGF_short <- EGF[RealTime >= 25 & RealTime <= 100]
FGF_short <- FGF[RealTime >= 25 & RealTime <= 100]
NGF_short <- NGF[RealTime >= 25 & RealTime <= 100]

clust_EGF <- CastCluster(EGF_short, time.col = "RealTime", condition.col = "Condition", label.col = "Label", k.clust = 2:8, na.fill = 1.1, plot = F, measure.col = "Ratio.norm.smooth")
clust_FGF <- CastCluster(FGF_short, time.col = "RealTime", condition.col = "Condition", label.col = "Label", k.clust = 2:8, na.fill = 1.1, plot = F, measure.col = "Ratio.norm.smooth")
clust_NGF <- CastCluster(NGF_short, time.col = "RealTime", condition.col = "Condition", label.col = "Label", k.clust = 2:8, na.fill = 1.1, plot = F, measure.col = "Ratio.norm.smooth")

plot_cluster(clust_EGF$table, id.vars.col = 1:9, cluster.col = "k_2", type = "trajectory") + ggtitle("2 Cluster trajectories - EGF")
plot_cluster(clust_EGF$table, id.vars.col = 1:9, cluster.col = "k_2", type = "composition") + ggtitle("2 Cluster compostions - EGF")

plot_cluster(clust_FGF$table, id.vars.col = 1:9, cluster.col = "k_2", type = "trajectory") + ggtitle("2 Cluster trajectories - FGF")
plot_cluster(clust_FGF$table, id.vars.col = 1:9, cluster.col = "k_2", type = "composition") + ggtitle("2 Cluster compostions - FGF")

plot_cluster(clust_NGF$table, id.vars.col = 1:9, cluster.col = "k_2", type = "trajectory") + ggtitle("2 Cluster trajectories - NGF")
plot_cluster(clust_NGF$table, id.vars.col = 1:9, cluster.col = "k_2", type = "composition") + ggtitle("2 Cluster compostions - NGF")

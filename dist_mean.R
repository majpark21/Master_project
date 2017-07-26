# 1) Get mean trajectory for a condition
# 2) Euclidean distance between each trajectory and the mean one
# 3) Distance distribution
library(data.table)
library(ggplot2)

dist_mean <- function(data, condition, tcol, measure, label, return.mean = F){
  # data: a data table containing trajectories and conditions
  # condition: column name that fully define an experimental condition, thus a set of trajectories analyzed together
  # tcol: time column name
  # measure: column name with quantity of interest
  # label: column name with label of individual objects in each condition (cell label)
  
  require(data.table)
  setkeyv(data, condition)
  
  # 1) Get mean trajectory for each condition
  means <- data[, .(Average = mean(get(measure))), by = .(get(condition), get(tcol))]
  colnames(means)[1:2] <- c(condition, tcol)
  dat <- merge(data[, c(condition, label, tcol, measure), with = F], means, by = c(condition, tcol))
  
  # 2) Euclidean distance between each trajectory and the mean one
  distances <- dat[, .(sq.dist = (get(measure) - Average)^2) , by = .(get(condition), get(label))]
  colnames(distances)[1:2] <- c(condition, label)
  
  euclid <- distances[, sqrt(sum(sq.dist)), by = .(get(condition), get(label))]
  colnames(euclid) <- c(condition, label, "euclid.to.mean")
  
  if(return.mean){
    return(list(means = means, euclid = euclid))
  } else {
    return(euclid)
  }
}



Cora <- fread("C:/Users/pixel/Dropbox/Marc-Antoine/data/set1-Coralie/tCoursesSelected.csv")
Cora[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
setkey(Cora, Image_Metadata_Site)
#means <- Cora[, .(mean.traj = mean(Ratio)), by = .(Image_Metadata_Site, RealTime),]


euclids <- dist_mean(Cora, "Image_Metadata_Site", "RealTime", "Ratio", "objNuc_TrackObjects_Label", return.mean = T)
#setkey(euclids, Image_Metadata_Site)

#plot(Cora[Image_Metadata_Site==0 & objNuc_TrackObjects_Label==3, Ratio], type = "b")
#plot(Cora[Image_Metadata_Site==0 & objNuc_TrackObjects_Label==17, Ratio], type = "b")

p <- ggplot(euclids$euclid, aes(euclid.to.mean)) + geom_histogram() + facet_grid(. ~ Image_Metadata_Site) + scale_x_continuous(limits = c(0, 15))
p

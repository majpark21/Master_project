library(data.table)
library(ggplot2)
library(wavelets)
library(TSdist)

# Dataset
Cora <- fread("C:/Users/pixel/Dropbox/Marc-Antoine/data/set1-Coralie/tCoursesSelected.csv")
# "Dummy clusters" to extract examples
Clust <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/temp/clust_hierch_data_euclidean_ward.D2.csv")

Clust[, c("Image_Metadata_Site", "objNuc_TrackObjects_Label") := strsplit(id, "_", fixed = T)]

# Trajectories
tr1 <- Cora[Well == 1 & objNuc_TrackObjects_Label == 1, objCyto_Intensity_MeanIntensity_imErkCorrOrig]
tr2 <- Cora[Well == 1 & objNuc_TrackObjects_Label == 2, objCyto_Intensity_MeanIntensity_imErkCorrOrig]

# Wavelet Transforms
wv1 <- dwt(tr1, filter = "haar")
wv2 <- dwt(tr2, filter = "haar")

# Fourier Transforms
fr1 <- fft(tr1)
fr2 <- fft(tr2)

# Cross Correlation based distance
CCorDistance(tr1, tr2)
DTWDistance(tr1, tr2)

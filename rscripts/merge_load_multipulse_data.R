library(data.table)
library(stringr)
# One table per frequency of pulse
x.1 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Cora_multipulse/1min.csv")
x.2 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Cora_multipulse/5min.csv")
x.3 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Cora_multipulse/10min.csv")
x.4 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Cora_multipulse/20min.csv")

# Create a Condition column that has all experimental condition
# Format: "PXX-IYY(-UO)"  X: pulse frequency; I: light intensity; UO: inhibited (optional)
create.condition <- function(data, pulse){
  require(data.table)
  require(stringr)
  pulse <- rep(pulse, nrow(data))
  light <- paste0(rep("I", nrow(data)), str_extract(data$Stimulation_intensity, "^[0-9]+"))
  inhib <- ifelse(data$Stimulation_treatment=="", "", "UO")
  condition <- paste(pulse, light, inhib, sep = "-")
  return(str_replace(condition, "-$", "")) # remove dashes at the end when no inhibitor
}

x.1$Condition <- create.condition(x.1, "P1"); x.1[, Condition := as.factor(Condition)]
x.2$Condition <- create.condition(x.2, "P5"); x.2[, Condition := as.factor(Condition)]
x.3$Condition <- create.condition(x.3, "P10"); x.3[, Condition := as.factor(Condition)]
x.4$Condition <- create.condition(x.4, "P20"); x.4[, Condition := as.factor(Condition)]

x.1[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]
x.2[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]
x.3[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]
x.4[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]

x.1[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
x.2[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
x.3[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
x.4[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]

# Merge and keep only columns of interest
Cora <- rbindlist(list(x.1, x.2, x.3, x.4))
Cora <- Cora[, c("Condition", "Label", "RealTime", "Ratio")]

# Column types
# Change order of levels in Condition (especially useful for plotting)
cond.good.order <- c("P1-I10", "P1-I25", "P1-I50", "P1-I100", "P1-I100-UO",
                     "P5-I10", "P5-I25", "P5-I50", "P5-I100", "P5-I100-UO",
                     "P10-I10", "P10-I25", "P10-I50", "P10-I100", "P10-I100-UO",
                     "P20-I10", "P20-I25", "P20-I50", "P20-I100", "P20-I100-UO")
Cora[, Condition := factor(Cora$Condition, levels = cond.good.order)]
Cora[, Label := as.factor(Label)]
Cora[, RealTime := as.numeric(RealTime)]
setkey(Cora, Condition, Label)

rm(cond.good.order, x.1, x.2, x.3, x.4)

library(data.table)
library(stringr)
# One table per length of pulse
x.1 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Coralie_singlePulse_doseResponse/tCoursesSelected_0050ms.csv")
x.2 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Coralie_singlePulse_doseResponse/tCoursesSelected_0100ms.csv")
x.3 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Coralie_singlePulse_doseResponse/tCoursesSelected_0500ms.csv")
x.4 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Coralie_singlePulse_doseResponse/tCoursesSelected_1000ms.csv")
x.5 <- fread("C:/Users/pixel/Desktop/Lectures/Semester_3/Master_thesis/data/Coralie_singlePulse_doseResponse/tCoursesSelected_5000ms.csv")

# Create a Condition column that has all experimental condition
# Format: "PXX-IYY(-UO)"  X: pulse duration (in s); I: light intensity; UO: inhibited (optional)
create.condition <- function(data, pulse){
  require(data.table)
  require(stringr)
  pulse <- rep(pulse, nrow(data))
  light <- paste0(rep("I", nrow(data)), str_extract(data$Stimulation_intensity, "^[0-9]+"))
  inhib <- ifelse(is.na(data$Stimulation_treatment), "", "UO")
  condition <- paste(pulse, light, inhib, sep = "-")
  return(str_replace(condition, "-$", "")) # remove dashes at the end when no inhibitor
}

x.1$Condition <- create.condition(x.1, "P0.05"); x.1[, Condition := as.factor(Condition)]
x.2$Condition <- create.condition(x.2, "P0.1"); x.2[, Condition := as.factor(Condition)]
x.3$Condition <- create.condition(x.3, "P0.5"); x.3[, Condition := as.factor(Condition)]
x.4$Condition <- create.condition(x.4, "P1"); x.4[, Condition := as.factor(Condition)]
x.5$Condition <- create.condition(x.5, "P5"); x.5[, Condition := as.factor(Condition)]

x.1[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]
x.2[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]
x.3[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]
x.4[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]
x.5[, Label := paste(Image_Metadata_Site, objNuc_TrackObjects_Label, sep="_")]

x.1[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
x.2[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
x.3[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
x.4[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
x.5[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]


# Merge and keep only columns of interest
Cora <- rbindlist(list(x.1, x.2, x.3, x.4, x.5))
Cora <- Cora[, c("Condition", "Label", "RealTime", "Ratio")]
# Add the condition missing to the grid:
temp <- length(unique(Cora$RealTime))
missing <- data.table(Condition = rep("P0.1-I0",temp), Label = rep("0_0",temp), RealTime = 0:(temp-1), Ratio = rep(0.35, temp))

Cora <- rbind(Cora, missing)

# Column types
# Change order of levels in Condition (especially useful for plotting)
cond.good.order <- c("P0.05-I0", "P0.05-I5", "P0.05-I10", "P0.05-I25", "P0.05-I50", "P0.05-I100",
                     "P0.1-I0", "P0.1-I5", "P0.1-I10", "P0.1-I25", "P0.1-I50", "P0.1-I100",
                     "P0.5-I0", "P0.5-I5", "P0.5-I10", "P0.5-I25", "P0.5-I50", "P0.5-I100",
                     "P1-I0", "P1-I5", "P1-I10", "P1-I25", "P1-I50", "P1-I100",
                     "P5-I0", "P5-I5", "P5-I10", "P5-I25", "P5-I50", "P5-I100")
Cora[, Condition := factor(Cora$Condition, levels = cond.good.order)]
Cora[, Label := as.factor(Label)]
Cora[, RealTime := as.numeric(RealTime)]
setkey(Cora, Condition, Label)
Cora <- myNorm(in.dt = Cora, in.meas.col = "Ratio", in.rt.min = 0, in.rt.max = 10, in.by.cols = c("Condition", "Label"), in.type = "fold.change")



rm(cond.good.order, x.1, x.2, x.3, x.4, x.5, missing, temp)

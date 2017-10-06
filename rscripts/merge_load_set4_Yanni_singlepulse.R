library(data.table)
library(stringr)
library(ggplot2)
source("../rscripts/myNorm.R")


# Read data
x.1 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170216_60min_noBSA_FGF250_25_2-5_0-25/tCoursesSelected.csv")
x.2 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170221_20min_FGF250_25_2-5_0-25_coll_ow//tCoursesSelected.csv")
x.3 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170223_3min_FGF250_25_2-5_0-25_coll/tCoursesSelected.csv")
x.4 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170224_5min_FGF250_25_2-5_0-25_coll/tCoursesSelected.csv")
x.5 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170629_10min_FGF250_25_2-5_0-25_coll/tCoursesSelected.csv")
x.6 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170912_10min_EGF250_25_2-5_0-25/tCoursesSelected.csv")
x.7 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170919_10min_NGF250_25_2-5_0-25/tCoursesSelected.csv")

# read stimulation times
y.1 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170216_60min_noBSA_FGF250_25_2-5_0-25/stimPulses.csv")
y.2 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170221_20min_FGF250_25_2-5_0-25_coll_ow//stimPulses.csv")
y.3 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170223_3min_FGF250_25_2-5_0-25_coll/stimPulses.csv")
y.4 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170224_5min_FGF250_25_2-5_0-25_coll/stimPulses.csv")
y.5 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170629_10min_FGF250_25_2-5_0-25_coll/stimPulses.csv")
y.6 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170912_10min_EGF250_25_2-5_0-25/stimPulses.csv")
y.7 <- fread("../data/set4_Yannick_EGFNGFFGFsinglepulse_PC12/20170919_10min_NGF250_25_2-5_0-25/stimPulses.csv")

Yanni_sp <- rbind(x.1, x.2, x.3, x.4, x.5, x.6, x.7)
StimYanni_sp <- rbind(y.1, y.2, y.3, y.4, y.5, y.6, y.7)

# ------------
# Create condition columns
gf <- str_extract(Yanni_sp$Stim_All_Ch, "[E,F,N]")
conc <- str_extract(Yanni_sp$Stim_All_Ch, "([ENF]GF +)(([0-9]+\\.[0-9]*)|([0-9]+))") %>% str_extract("(([0-9]+\\.[0-9]*)|([0-9]+))")
pulse <- str_extract(Yanni_sp$Stim_All_Ch, "[0-9]+(min|\')") %>% str_extract("[0-9]+")
Yanni_sp$Condition <- paste(paste0(gf,conc), paste0("P", pulse), sep = "-")

gf <- str_extract(StimYanni_sp$Stim_All_Ch, "[E,F,N]")
conc <- str_extract(StimYanni_sp$Stim_All_Ch, "([ENF]GF +)(([0-9]+\\.[0-9]*)|([0-9]+))") %>% str_extract("(([0-9]+\\.[0-9]*)|([0-9]+))")
pulse <- str_extract(StimYanni_sp$Stim_All_Ch, "[0-9]+(min|\')") %>% str_extract("[0-9]+")
StimYanni_sp$Condition <- paste(paste0(gf,conc), paste0("P", pulse), sep = "-")

# ------------
# Modif columns
setnames(Yanni_sp, old = c("TrackObjects_Label_uni", "Intensity_MeanIntensity_Ratio"), new = c("Label", "Ratio"))
Yanni_sp <- Yanni_sp[,-c(1,2,3,7,8)]
setkey(Yanni_sp, Condition, Label)
setcolorder(Yanni_sp, c("Condition", "Label", "RealTime", "Ratio"))

setnames(StimYanni_sp, old = c("Intensity_MeanIntensity_pulse"), new = c("Pulse_intensity"))
StimYanni_sp <- StimYanni_sp[,-c(1,2,3,5,7,8)]
setkey(StimYanni_sp, Condition)
setcolorder(StimYanni_sp, c("Condition", "RealTime", "Pulse_intensity"))

# ------------
# Normalization
Yanni_sp <- myNorm(in.dt = Yanni_sp, in.meas.col = "Ratio", in.rt.min = 0, in.rt.max = 36, in.by.cols = c("Condition", "Label"), in.type = "fold.change")
StimYanni_sp <- myNorm(in.dt = StimYanni_sp, in.meas.col = "Pulse_intensity", in.rt.min = 0, in.rt.max = 36, in.by.cols = c("Condition"), in.type = "fold.change")

# Outlier clipping and column order
Yanni_sp[which(is.na(Yanni_sp$Ratio)), c("Ratio", "Ratio.norm") := list(1215, 1)]
temp <- which(Yanni_sp$Ratio.norm > 2.5)
Yanni_sp[temp, c("Ratio", "Ratio.norm") := list(1215, 1)]
cond.good.order <- c(unique(Yanni_sp$Condition)[1:4],
                     unique(Yanni_sp$Condition)[25:28],
                     unique(Yanni_sp$Condition)[seq(5, 24, 5)],
                     unique(Yanni_sp$Condition)[seq(7, 24, 5)],
                     unique(Yanni_sp$Condition)[seq(8, 24, 5)],
                     unique(Yanni_sp$Condition)[seq(6, 24, 5)],
                     unique(Yanni_sp$Condition)[seq(9, 24, 5)])
Yanni_sp[, Condition := factor(Yanni_sp$Condition, levels = cond.good.order)]

# --------------
# Clean
rm(x.1, x.2, x.3, x.4, x.5, x.6, x.7,
   y.1, y.2, y.3, y.4, y.5, y.6, y.7,
   gf, conc, pulse, cond.good.order, temp)

###########################
#    dist_mean     #
###########################

Cora <- fread("C:/Users/pixel/Dropbox/Marc-Antoine/data/set1-Coralie/tCoursesSelected.csv")
Cora[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
setkey(Cora, Image_Metadata_Site)
#means <- Cora[, .(mean.traj = mean(Ratio)), by = .(Image_Metadata_Site, RealTime),]


euclids <- dist_mean(Cora, "Image_Metadata_Site", "RealTime", "Ratio", "objNuc_TrackObjects_Label", return.mean = F)
#setkey(euclids, Image_Metadata_Site)
euclids[, .(Mean = mean(euclid_to_mean), Variance = var(euclid_to_mean), Min = min(euclid_to_mean), Max = max(euclid_to_mean)), Image_Metadata_Site]


p <- ggplot(euclids$euclid, aes(x=Image_Metadata_Site, y=euclid_to_mean)) + geom_boxplot(aes(group=Image_Metadata_Site)) + scale_y_continuous(limits = c(0,4.5))
p


###########################
#    Overlap clipping     #
###########################

library(data.table)
library(ggplot2)
Cora <- fread("C:/Users/pixel/Dropbox/Marc-Antoine/data/set1-Coralie/tCoursesSelected.csv")
Cora[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
setkey(Cora, Image_Metadata_Site, objNuc_TrackObjects_Label)

ClipRatio <- Cora[, .(clip_ratio = wrap_clip(Ratio, k = 5)), by = .(Image_Metadata_Site, objNuc_TrackObjects_Label)]
Overlap <- overlap_clipping(data = ClipRatio, condition = "Image_Metadata_Site", label = "objNuc_TrackObjects_Label", measure = "clip_ratio")

##### 
#Plot example of a trajectory, rolling mean and clipped trajectory
#plot(ClipRatio[Image_Metadata_Site==5 & objNuc_TrackObjects_Label==2, clip_ratio], type = 'l', col = 'blue')
#points(Cora[Image_Metadata_Site == 5 & objNuc_TrackObjects_Label == 2, Ratio])
#lines(Cora[Image_Metadata_Site == 5 & objNuc_TrackObjects_Label == 2, Ratio])
#lines(rollex(Cora[Image_Metadata_Site == 5 & objNuc_TrackObjects_Label == 2, Ratio]), col = 'red')

#ClipRatio$RealTime <- Cora$RealTime
#p <- ggplot(ClipRatio[.(7,2)], aes(x=RealTime, y=clip_ratio)) + geom_step(alpha = 1) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); p
#####


overlap(ClipRatio[.(7,1), clip_ratio], ClipRatio[.(7,42), clip_ratio])
cor(ClipRatio[.(7,3), clip_ratio],ClipRatio[.(7,42), clip_ratio], method = "k")



p <- ggplot(Overlap, aes(x = Image_Metadata_Site, y = Overlap)) + geom_boxplot(aes(group  = Image_Metadata_Site))
p


temp <- overlap_clipping(data = ClipRatio, condition = "Image_Metadata_Site", label = "objNuc_TrackObjects_Label", measure = "clip_ratio")








###########################
#         sim data        #
###########################

noises <- seq(0.2, 3, 0.2)
multi_sim <-  multi_sims(type="na", noises=noises, n = 50)
plot_sim(multi_sim)

DistMean <- dist_mean(data = multi_sim, condition = "noise", tcol = "Time", measure = "value", label = "variable")
p <- ggplot(DistMean, aes(x = as.factor(noise), y = euclid_to_mean)) + geom_boxplot(aes(group = noise)) + ggtitle("Euclidian distance to mean trajectory"); p

Clip <- multi_sim[, .(clip = wrap_clip(value, k = 5)), by = .(noise, variable)]
Overlap <- overlap_clipping(data = Clip, condition = "noise", label = "variable", measure = "clip")
q <- ggplot(Overlap, aes(x = as.factor(noise), y = Overlap)) + geom_boxplot(aes(group  = noise)) + scale_x_discrete(labels = as.character(c(0,noises))) + ggtitle("Pairwise overlap of clipped trajectories") ; q


pdf("Plots_Large_Sim.pdf", width = 10)
plot_sim(multi_sim)
p
q
dev.off()


multi_sim[, ':=' (noise = as.integer(noise*10),
                  variable = as.integer(gsub("V", "", as.character(variable), fixed = T)))]
Correlations_Pearson <- correlations_group_label(multi_sim, "noise", "variable","value", method = "pearson")
Correlations_Spearman <- correlations_group_label(multi_sim, "noise", "variable","value", method = "spearman")
Correlations_Kendall <- correlations_group_label(multi_sim, "noise", "variable","value", method = "kendall")






###########################
#  apply stats to data    #
###########################

library(data.table)
library(ggplot2)
library(plotly)

Cora <- fread("C:/Users/pixel/Dropbox/Marc-Antoine/data/set1-Coralie/tCoursesSelected.csv")
Cora[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
setkey(Cora, Image_Metadata_Site, objNuc_TrackObjects_Label)


cond <- "Image_Metadata_Site"
tcol <- "RealTime"
mea <- "Ratio"
lab <- "objNuc_TrackObjects_Label"

Cora_mean <- dist_mean(data = Cora, condition = cond, tcol = tcol, measure = mea, label = lab, return.mean = F)
Cora_amp <- amplitude_oscillations(data = Cora, condition = cond, measure = mea, label = lab, k_roll_mean = 5)
Cora_pw <- all_pairwise_stats(data = Cora, condition = cond, label = lab, measure = mea, k_roll_mean = 5)
Cora_pw_long <- melt(Cora_pw, id.vars = c("Image_Metadata_Site", "Label1", "Label2"))

# Mean plot - violin
p1 <- ggplot(Cora_mean[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = as.factor(Image_Metadata_Site), y = euclid_to_mean, text = paste("Label:", objNuc_TrackObjects_Label))) +
  geom_violin(aes(group = as.factor(Image_Metadata_Site))) +
  geom_jitter(alpha=1, aes(col = as.factor(Image_Metadata_Site)), width = 0.15) +
  theme(legend.position="none") +
  ggtitle("Euclidean Distance to Average Trajectory")
p2 <- ggplotly(p1)

# PW plot - Boxplot
q1 <- ggplot(Cora_pw_long[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = as.factor(Image_Metadata_Site), y = value, text = paste("Label1:", Label1, "; Label2:", Label2 ))) +
  geom_boxplot(aes(group = as.factor(Image_Metadata_Site))) + 
  geom_point(alpha = 0) +
  facet_grid(. ~ variable) +
  ggtitle("Pairwise metrics: Correlations and Overlap of Clipped Trajectories")
q2 <- ggplotly(q1)

# Mean plot - density
r1 <- ggplot(Cora_mean[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = euclid_to_mean)) +
  geom_density() +
  facet_grid(Image_Metadata_Site ~ .) + 
  ggtitle("Euclidean Distance to Average Trajectory")

r2 <- ggplot(Cora_mean[Image_Metadata_Site %in% c(1,3,5,7)], aes(x= euclid_to_mean, colour = as.factor(Image_Metadata_Site), fill = as.factor(Image_Metadata_Site))) +
  geom_density(size = 1.75, alpha = 0.1) +
  ggtitle("Euclidean Distance to Average Trajectory") +
  theme(legend.position="bottom")

r3 <- ggplot(Cora_mean[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = as.factor(Image_Metadata_Site), y = euclid_to_mean, text = paste("Label:", objNuc_TrackObjects_Label))) +
  geom_boxplot(aes(group = as.factor(Image_Metadata_Site))) +
  theme(legend.position="none") +
  ggtitle("Euclidean Distance to Average Trajectory")

# PW plot - Densities
s1 <- ggplot(Cora_pw_long[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = value)) +
  geom_density() + 
  facet_wrap(Image_Metadata_Site ~ variable) + 
  ggtitle("Pairwise metrics: Correlations and Overlap of Clipped Trajectories")




t1 <- ggplot(Cora_amp[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = as.factor(Image_Metadata_Site), y = euclid_to_roll_mean, text = paste("Label:", objNuc_TrackObjects_Label))) +
  geom_violin(aes(group = as.factor(Image_Metadata_Site))) +
  geom_jitter(alpha=1, aes(col = as.factor(Image_Metadata_Site)), width = 0.15) +
  theme(legend.position="none") +
  ggtitle("Euclidean Distance to individual Rolling Mean", "Indicates depth of oscillations")

t2 <- ggplot(Cora_amp[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = as.factor(Image_Metadata_Site), y = euclid_to_roll_mean, text = paste("Label:", objNuc_TrackObjects_Label))) +
  geom_boxplot(aes(group = as.factor(Image_Metadata_Site))) +
  theme(legend.position="none") +
  ggtitle("Euclidean Distance to individual Rolling Mean", "Indicates depth of oscillations")

t3 <- ggplot(Cora_amp[Image_Metadata_Site %in% c(1,3,5,7)], aes(x = euclid_to_roll_mean)) +
  geom_density() + 
  facet_grid(Image_Metadata_Site ~ .) + 
  ggtitle("Euclidean Distance to individual Rolling Mean", "Indicates depth of oscillations")

t4 <- ggplot(Cora_amp[Image_Metadata_Site %in% c(1,3,5,7)], aes(x= euclid_to_roll_mean, colour = as.factor(Image_Metadata_Site), fill = as.factor(Image_Metadata_Site))) +
  geom_density(size = 1.75, alpha = 0.1) +
  ggtitle("Euclidean Distance to individual Rolling Mean", "Indicates depth of oscillations") +
  theme(legend.position="bottom")


pdf("Coherence_stats_Cora.pdf", height = 14, width = 20)
grid.arrange(s1,q1, ncol=2)
grid.arrange(p1, r3, r1,r2, ncol = 2)
grid.arrange(t1, t2, t3, t4, ncol = 2)
dev.off()


###############################
#   Amplitude oscillation     #
###############################
temp = amplitude_oscillations(Cora, "Image_Metadata_Site", "objNuc_TrackObjects_Label", "Ratio", 5)

p <- ggplot(data = temp, aes(x = Image_Metadata_Site, y = euclid_to_roll_mean, text=objNuc_TrackObjects_Label)) + geom_boxplot(aes(group = Image_Metadata_Site)) + geom_point(alpha=0)
p <- ggplotly(p)

conds <- 0:7
ks.p.val <- cbind(t(combn(conds, 2)), apply(combn(conds, 2), 2, function(x) ks.test(temp[Image_Metadata_Site==x[1], euclid_to_roll_mean], temp[Image_Metadata_Site==x[2], euclid_to_roll_mean])$p.value ))
ks.p.val[,3] <- p.adjust(ks.p.val[,3], method = "holm")

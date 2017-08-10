get.max.conv <- function(ts1, ts2, type = "circular", plot = T){
  temp <- convolve(ts1, ts2, type = type)
  if(plot){plot(temp, type = "b", xlab = "Lag", ylab = "Convolution")}
  return(list(temp[which.max(temp)], which.max(temp) - 1))
}

stim <- rep(0, 90)
stim[seq(11,66,5)] <- 1
conv.stim <- copy(Cora)
conv.stim[, c("max.conv", "lag") := get.max.conv(Ratio, stim, plot = F), by = .(Image_Metadata_Site, objNuc_TrackObjects_Label)]

p <- ggplot(conv.stim[Image_Metadata_Site %in% seq(1,7,2)], aes(x = as.factor(Image_Metadata_Site), y = lag)) + geom_boxplot(aes(group=as.factor(Image_Metadata_Site)))
q <- ggplotly(p)

library(lmtest)
tr1 <- Cora[.(7,3), Ratio]

grangertest(tr1 ~ stim, order=10)
grangertest(stim ~ tr1, order=6)

ccf(tr1, stim)
plot(convolve(tr1, stim))

for(i in 0:7){
  print(convolve(Cora[.(i,4), Ratio], stim, type="filter"))
}


filt.stim <- Cora[, .(filt.stim = convolve(Ratio, stim, type="filter")), .(Image_Metadata_Site, objNuc_TrackObjects_Label)]
filt.stim[, c("Image_Metadata_Site", "objNuc_TrackObjects_Label") := list(as.factor(Image_Metadata_Site), as.factor(objNuc_TrackObjects_Label))]

library(plotly)
p <- ggplot(filt.stim, aes(x=Image_Metadata_Site, y=filt.stim)) + geom_boxplot(aes(group=Image_Metadata_Site))
q <- ggplotly(p)



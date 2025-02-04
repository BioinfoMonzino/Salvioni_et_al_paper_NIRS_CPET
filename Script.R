rm(list = ls())
graphics.off()
cat("\014")

setwd("U:/LAVORO.Analisi/ANALISI/ANALISI.BioAI/AgostoniP/004.NIRS.Hb.muscolo.e.Respiro.Periodico/202405.withStartMarkers/")


library(RColorBrewer)
library(Hmisc)
library(pheatmap)
library(corrplot)
library(calibrate)
library(rgl)
library(car)
library(cluster)
library(ggrepel)
library(ggridges)
library(caret)
library(MASS)
library(reshape2)
library(ggplot2)
library(MatrixGenerics)
library(pracma)
# library(data.table)
library(gsignal)
library(plotly)



source("./0.Scripts/custom_functions.R")
source("./0.Scripts/disegna_multiplot.R")

# custom functions
'%!in%' <- function(x,y)!('%in%'(x,y))

##########################################################################################################
 patients_id <- "s021"

###########################################################################################################
#### Import
df_delay_metadata <- read.delim("./Metadata.txt",stringsAsFactors = T, na.strings = "")

# df must have the Time in the first column named "Time_s"
filename_nirs <- paste0("./",patients_id,"/NIRS.txt")
filename_cpet <- paste0("./",patients_id,"/CPET.txt")

df.nirs <- read.delim(filename_nirs,na.strings = "")
df.cpet <- read.delim(filename_cpet,na.strings = "")

idxlast_var_nirs <- ncol(df.nirs)-1
idxlast_var_cpet <- ncol(df.cpet)-1


######################################################################################################
#### convert NIRS in seconds

df.nirs <- df.nirs[-1,] # removing the first 'highly noisy' measurement
n_sec <- nrow(df.nirs) %/% 10
df.nirs <- df.nirs[1:(n_sec*10),]

df.nirs.sec <- as.data.frame(matrix(nrow = n_sec,ncol = ncol(df.nirs)))
colnames(df.nirs.sec) <- colnames(df.nirs)

# set marker
idx_new_mark <- which(!is.na(df.nirs$Marker)) %/% 10
df.nirs.sec[idx_new_mark,"Marker"] <- "start"

# average measurements every 10 dec sec
k <- 1
for (i in 1:nrow(df.nirs.sec)){
  
  df.nirs_sub <- df.nirs[k:(k+9),2:idxlast_var_nirs]
  df.nirs.sec[i,2:idxlast_var_nirs] <- colMeans(df.nirs_sub)
  
  k <- k + 10
}

# signals have to start from Time = 0
df.nirs.sec[,1] <- 0:(nrow(df.nirs.sec)-1)
df.nirs <- df.nirs.sec
rm(df.nirs.sec)

df.cpet$Time_s <- df.cpet$Time_s-min(df.cpet$Time_s)

# add costant to avoid NIRS signal with negative values
if(min(df.nirs[,2:3]) < 0){
  
  mins_vect_nirs <- min(df.nirs[,2:3])
  df.nirs[,2:3] <- df.nirs[,2:3] + abs(min(mins_vect_nirs))
  # plots
  df.nirs.ggplot <- melt(df.nirs,measure.vars = colnames(df.nirs)[2:idxlast_var_nirs])
  df.cpet.ggplot <- melt(df.cpet,measure.vars = colnames(df.cpet)[2:idxlast_var_cpet])
  
  p1 <- ggplot(df.nirs.ggplot,aes(x=Time_s,y=value,color=variable))+
    geom_line()+
    geom_point() +
    scale_color_manual(values = c("red","blue"))+
    ggtitle("NIRS")
  
  p2 <- ggplot(df.cpet.ggplot,aes(x=Time_s,y=value,color=variable))+
    geom_line()+
    geom_point() +
    ggtitle("CPET")
  
  
multiplot(p1,p2,cols = 1) #Raw signals
  
}

###################################################################################################
### calculate total and diff Hb
df.nirs$tHb <- df.nirs$O2Hb + df.nirs$HHb
df.nirs$diffHb <- df.nirs$O2Hb - df.nirs$HHb
idxlast_var_nirs <- ncol(df.nirs)-1
idxlast_var_cpet <- ncol(df.cpet)-1


######################################################################################################
####  Signal interpolation 
df.nirs.interp <- check_n_interpolate(df.nirs[,c(1,2,3,5,6)],"linear")
df.cpet.interp <- check_n_interpolate(df.cpet[,c(1,2,3,4,5,6)],"linear")


#####################################################################################################
#### CPET and NIRS synchronization
df.nirs.interp$Time_s_sync <- 0

T_marker_nirs <- df.nirs$Time_s[which(!is.na(df.nirs$Marker))]
T_marker_cpet <- df.cpet$Time_s[which(!is.na(df.cpet$Marker))]
T_marker_nirs_new <- T_marker_cpet + df_delay_metadata$s_delay_NIRS_CPET[which(rownames(df_delay_metadata) %in% patients_id)]

new_Time_nirs <- (T_marker_nirs_new - T_marker_nirs):(max(df.nirs.interp$Time_s)+(T_marker_nirs_new- T_marker_nirs))
df.nirs.interp$Time_s_sync <- new_Time_nirs

Time_shared <- sort(union(df.cpet.interp$Time_s,df.nirs.interp$Time_s_sync))

df_CPET_NIRS.synced <- as.data.frame(matrix(nrow = length(Time_shared), ncol = 10))
colnames(df_CPET_NIRS.synced) <- c("Time_s","VE","VO2","VCO2","VT","HR", "O2Hb","HHb","tHb","diffHb")
df_CPET_NIRS.synced$Time_s <- Time_shared

t_cpet_inters <- intersect(df_CPET_NIRS.synced$Time_s,df.cpet.interp$Time_s)
t_nirs_inters <- intersect(df_CPET_NIRS.synced$Time_s,df.nirs.interp$Time_s_sync)

df_CPET_NIRS.synced[which(df_CPET_NIRS.synced$Time_s %in% t_cpet_inters),c("VE","VO2","VCO2", "VT","HR")] <- df.cpet.interp[which(df.cpet.interp$Time_s == t_cpet_inters),c("VE","VO2","VCO2","VT","HR")]
df_CPET_NIRS.synced[which(df_CPET_NIRS.synced$Time_s %in% t_nirs_inters),c("O2Hb","HHb","tHb","diffHb")] <- df.nirs.interp[which(df.nirs.interp$Time_s_sync == t_nirs_inters),c("O2Hb","HHb","tHb","diffHb")]

# plots
df_CPET_NIRS.synced.ggplot.nirs <- melt(df_CPET_NIRS.synced[,c(1,7,8)],measure.vars = c("O2Hb","HHb"))
df_CPET_NIRS.synced.ggplot.cpet <- melt(df_CPET_NIRS.synced[,c(1,2,3,4,5,6)],measure.vars = c("VE","VO2","VCO2", "VT","HR"))

p7 <- ggplot(df_CPET_NIRS.synced.ggplot.nirs,aes(x=Time_s,y=value,color=variable))+
  geom_line()+
  geom_point() +
  scale_color_manual(values = c("red","blue"))+
  ggtitle("NIRS")

p8 <- ggplot(df_CPET_NIRS.synced.ggplot.cpet,aes(x=Time_s,y=value,color=variable))+
  #scale_color_manual(values = c("green","violet"))+
  geom_line()+
  geom_point() +
  ggtitle("CPET")

#multiplot(p7,p8,cols = 1)

###################################################################################################
######### cutting signals, selecting oscillating phase
cut_time_range <- c(30,400) 
df_CPET_NIRS.synced.cut <- df_CPET_NIRS.synced[which(df_CPET_NIRS.synced$Time_s %in% cut_time_range[1]:cut_time_range[2]),]

df_CPET_NIRS.synced.cut.ggplot.nirs <- melt(df_CPET_NIRS.synced.cut[,c(1,7,8)],measure.vars = c("O2Hb","HHb"))
df_CPET_NIRS.synced.cut.ggplot.cpet <- melt(df_CPET_NIRS.synced.cut[,c(1,2,3,4,5,6)],measure.vars = c("VE","VO2","VCO2", "VT","HR"))

# plot
p9 <- ggplot(df_CPET_NIRS.synced.cut.ggplot.nirs,aes(x=Time_s,y=value,color=variable))+
  geom_line()+
  geom_point() +
  scale_color_manual(values = c("red","blue"))+
  ggtitle("NIRS: synced and cut signal")

p10 <- ggplot(df_CPET_NIRS.synced.cut.ggplot.cpet,aes(x=Time_s,y=value,color=variable))+
  #scale_color_manual(values = c("green","violet"))+
  geom_line()+
  geom_point() +
  ggtitle("CPET: synced and cut signal")

#multiplot(p9,p10,cols = 1)

#######################################################################################################
######## Smoothing
df_CPET_NIRS.synced.cut.smoth <- ave_med_filt_df(df_CPET_NIRS.synced.cut,"average",10)

# plots
df_CPET_NIRS.synced.cut.smoth.ggplot.nirs <- melt(df_CPET_NIRS.synced.cut.smoth[,c(1,7,8)],measure.vars = c("O2Hb","HHb"))
df_CPET_NIRS.synced.cut.smoth.ggplot.cpet <- melt(df_CPET_NIRS.synced.cut.smoth[,c(1,2,3,4,5,6)],measure.vars = c("VE","VO2","VCO2", "VT","HR"))

p11 <- ggplot(df_CPET_NIRS.synced.cut.smoth.ggplot.nirs,aes(x=Time_s,y=value,color=variable))+
  geom_line()+
  geom_point() +
  scale_color_manual(values = c("red","blue"))+
  ggtitle("NIRS: synced, cut and smoothed signal")

p12 <- ggplot(df_CPET_NIRS.synced.cut.smoth.ggplot.cpet,aes(x=Time_s,y=value,color=variable))+
  #scale_color_manual(values = c("green","violet"))+
  geom_line()+
  geom_point() +
  ggtitle("CPET: synced, cut and smoothed signal") + facet_wrap(~variable)

#multiplot(p11,p12,cols = 1)



#####################################################################################################
# Scaling
df_CPET_NIRS.synced.cut.smoth.scale <- as.data.frame(cbind(df_CPET_NIRS.synced.cut.smoth[,1,drop=F],scale(df_CPET_NIRS.synced.cut.smoth[,2:ncol(df_CPET_NIRS.synced.cut.smoth)])))

##### Plots
df_CPET_NIRS.synced.cut.smoth.scale.ggplot.nirs <- melt(df_CPET_NIRS.synced.cut.smoth.scale,measure.vars = c("O2Hb","HHb"))
#df_CPET_NIRS.synced.cut.smoth.scale.ggplot.cpet <- melt(df_CPET_NIRS.synced.cut.smoth.scale,measure.vars = c("VE","VO2","VCO2", "VT","HR"))
df_CPET_NIRS.synced.cut.smoth.scale.ggplot.cpet <- melt(df_CPET_NIRS.synced.cut.smoth.scale,measure.vars = c("VE","VO2","VCO2", "VT"))

p13 <- ggplot(df_CPET_NIRS.synced.cut.smoth.scale.ggplot.nirs,aes(x=Time_s,y=value,color=variable))+
  geom_line()+
  geom_point() +
  scale_color_manual(values = c("red","blue"))+
  ggtitle("NIRS: synced, cut, smoothed and scaled signal")

p14 <- ggplot(df_CPET_NIRS.synced.cut.smoth.scale.ggplot.cpet,aes(x=Time_s,y=value,color=variable))+
  #scale_color_manual(values = c("green","violet"))+
  geom_line()+
  geom_point() +
  ggtitle("CPET: synced, cut, smoothed and scaled signal")



multiplot(p14,p13,cols = 1) # Preprocessed Signals


############################################################################################################
####################             Cross- correlation           #####################
df_CPET.fase <- df_CPET_NIRS.synced.cut.smoth[,c(1,2,3,4,5,6)]

# scale
df_CPET.fase.scale <- as.data.frame(cbind(df_CPET.fase[,1,drop=F],
                                          scale(df_CPET.fase[,2:ncol(df_CPET.fase)])))

time_range_cross_corr <- 100 * floor((cut_time_range[2]-cut_time_range[1]) / 100) -50

signal_name <- "VO2"
#signal_name <- "VCO2"
#  signal_name <- "VE"
# signal_name <- "VT"
#signal_name <- "HR"


cross_corr <- cross_corr_signals_v2(df_CPET.fase.scale,
                                    ref_sig = signal_name,
                                    test_sig = signal_name,
                                    time_range = c(-time_range_cross_corr:time_range_cross_corr)
)

gg1 <- ggplot(cross_corr,aes(x=Lags, y=CrossCorrelation))+
  geom_line(color="orange4")+
  geom_point(size=3,color="orange")+
  geom_vline(xintercept = 0,linetype="dashed", color = "red") +
  geom_hline(yintercept = 0,linetype="dashed", color = "red")+
  ylim(-1,1) + 
  ggtitle(signal_name)
#ggplotly(gg1)
gg1

########## calculate phase degree
phase_in_seconds = 39
avg_cycle_duration=80
phase_degree = phase_in_seconds * 360 / avg_cycle_duration

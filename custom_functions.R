###########################################################################
check_n_interpolate <- function(df,interp_type=c("linear", "constant", "nearest", "spline", "cubic")){
  # interpola il segnale per avere un segnale ogni secondo
  
  if(max(df$Time_s) != (dim(df)[1]-1)){
    
    time.step.new <- 0:max(df$Time_s)
    df.new <- as.data.frame(matrix(nrow = length(time.step.new),ncol = dim(df)[2]))
    colnames(df.new) <- colnames(df)
    df.new[,1] <- time.step.new
    
    for (i in 2:dim(df)[2]){
      
      df.new[,i] <- interp1(x = df$Time_s,y = df[,i],xi = time.step.new,method=interp_type)
      
    }
    return(df.new)
  }else{
    return(df)
    print("The dataframe has not to be interploated")
  }
}

##############################################################################Ã
ave_med_filt_df <- function(df,type=c("median", "average"),n_window=5){
  
  time.step.new <- df$Time_s
  df.new <- as.data.frame(matrix(nrow = length(time.step.new),ncol = dim(df)[2]))
  colnames(df.new) <- colnames(df)
  df.new[,1] <- time.step.new  
  
  
  ##### filters
  #Average Filter
  average_filt <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=2)} 
  #Median Filter
  median_filt <- function(x,n=5){runmed(x,n)} 

  
  for (i in 2:dim(df)[2]){
    if(type == "median"){
      df.new[,i] <- median_filt(df[,i],n=n_window)
      
    }else if(type == "average"){
      df.new[,i] <- average_filt(df[,i],n=n_window)
    }
  }
  
  df.new <- df.new[complete.cases(df.new),]
  
  return(df.new)
}

################################################################################
cross_corr_signals <- function(df,ref_sig,test_sig,time_range=c(-10:10)){
  
  # df deve avere 2 colonne e basta
  
  
  ### cross correlation
  cor.cross <- function(x0, y0, i=0) {
    #
    # Sample autocorrelation at (integral) lag `i`:
    # Positive `i` compares future values of `x` to present values of `y`';
    # negative `i` compares past values of `x` to present values of `y`.
    #
    if (i < 0) {x<-y0; y<-x0; i<- -i}
    else {x<-x0; y<-y0}
    n <- length(x)
    cor(x[(i+1):n], y[1:(n-i)], use="complete.obs")
  }
  ####
  df_1 <- df[,ref_sig]
  df_2 <- df[,test_sig]
  
  data.cor <- round(sapply(time_range, function(i) cor.cross(df_1, df_2, i)),2)
  
  df_plot <- as.data.frame(cbind(time_range,data.cor))
  colnames(df_plot) <- c("Lags","CrossCorrelation")
  
  best_x_corr_pos <- min(time_range[which(data.cor ==  max(data.cor)[1])])
  best_y_corr_pos <- max(data.cor)[1]
  
  best_x_corr_neg <- max(time_range[which(data.cor ==  min(data.cor)[1])])
  best_y_corr_neg <- min(data.cor)[1]
  
  
  gg1 <- ggplot(df_plot,aes(x=Lags, y=CrossCorrelation))+
          geom_line(color="orange4")+
          geom_point(size=3,color="orange")+
          geom_point(x=best_x_corr_pos,y=best_y_corr_pos,color="purple",size=5) +
          geom_point(x=best_x_corr_neg,y=best_y_corr_neg,color="purple",size=5) +
          geom_point(x=0,y=df_plot[which(df_plot$Lags ==0),2],color="green",size=5) +
          geom_vline(xintercept = 0,linetype="dashed", color = "red") +
          geom_hline(yintercept = 0,linetype="dashed", color = "red")+
          ylim(-1,1)+
          ggtitle(paste0(ref_sig," vs. ",test_sig))
        
        
        
  print(gg1)
  # print(paste0("Current Signal correlation = ",df_plot[which(df_plot$Lags ==0),2] ))
  # print(paste0("The ",test_sig," signal should be shifted by ",best_x_corr_pos, " seconds to reach the highest correlation: ", best_y_corr_pos))
  # print(paste0("The ",test_sig," signal should be shifted by ",abs(best_x_corr_neg), " seconds to reach the lowest correlation ", best_y_corr_neg))
  # 
  
  # plots
  #si muove il test e resta fermo il ref 
  
  # df_tr_1 <- as.data.frame(matrix(nrow = abs(time_range[1]), ncol = ncol(df)))
  # df_tr_2 <- as.data.frame(matrix(nrow = abs(time_range[2]), ncol = ncol(df)))
  # colnames(df_tr_1) <- colnames(df_tr_2) <- colnames(df)
  # 
  # df_tr_1$Time_s <- (min(df$Time_s)-abs(time_range[1])):(min(df$Time_s)-1)
  # df_tr_2$Time_s <- (max(df$Time_s)+1): (max(df$Time_s)+abs(time_range[2]))
  # 
  # df_2 <- rbind(df_tr_1,df,df_tr_2)
  # 
  # df_2$test_phase <- NA
  # df_2$test_countphase <- NA
  # 
  # 
  # time_range_corr_pos <- (min(df$Time_s)+best_x_corr_pos):(max(df$Time_s)+best_x_corr_pos)
  # time_range_corr_neg <- (min(df$Time_s)+best_x_corr_neg):(max(df$Time_s)+best_x_corr_neg)
  # 
  # 
  # 
  # df_2$test_phase[which(df_2$Time_s %in% time_range_corr_pos)] <- df_2[!is.na(df_2[,test_sig]),test_sig]
  # df_2$test_countphase[which(df_2$Time_s %in% time_range_corr_neg)] <- df_2[!is.na(df_2[,test_sig]),test_sig]
  # 
  # 
  # 
  # df_ggplot <- melt(df,measure.vars = c(ref_sig,test_sig))
  # print(ggplot(df_ggplot,aes(x=Time_s,y=value,color=variable))+
  #   geom_line()+
  #   geom_point() +
  #   scale_color_manual(values = c("red","blue"))+
  #   ggtitle(paste0("No lag. R = ",df_plot[which(df_plot$Lags ==0),2] ))
  # )
  
  
  
  
  # 
  # df_ggplot <- melt(df_2,measure.vars = c(ref_sig,test_sig))
  # p0 <- ggplot(df_ggplot,aes(x=Time_s,y=value,color=variable))+
  #   geom_line()+
  #   geom_point() +
  #   scale_color_manual(values = c("red","blue"))+
  #   ggtitle(paste0("No lag. R = ",df_plot[which(df_plot$Lags ==0),2] ))
  # 
  # 
  # 
  # 
  # df_phase <- df_2[,c("Time_s",ref_sig,"test_phase")]
  # colnames(df_phase) <- c("Time_s",ref_sig,test_sig)
  # df_phase.ggplot <- melt(df_phase,measure.vars = c(ref_sig,test_sig))
  # p1 <- ggplot(df_phase.ggplot,aes(x=Time_s,y=value,color=variable))+
  #   geom_line()+
  #   geom_point() +
  #   scale_color_manual(values = c("red","blue"))+
  #   ggtitle(paste0("Highest correlation: lag = ",best_x_corr_pos," sec. R = ",best_y_corr_pos))
  # 
  # 
  # df_countphase <- df_2[,c("Time_s",ref_sig,"test_countphase")]
  # colnames(df_countphase) <- c("Time_s",ref_sig,test_sig)
  # 
  # df_count_phase.ggplot <- melt(df_countphase,measure.vars = c(ref_sig,test_sig))
  # p2 <- ggplot(df_count_phase.ggplot,aes(x=Time_s,y=value,color=variable))+
  #   geom_line()+
  #   geom_point() +
  #   scale_color_manual(values = c("red","blue"))+
  #   ggtitle(paste0("Lowest correlation: lag= ",best_x_corr_neg," sec. R = ",best_y_corr_neg))
  # 
  # multiplot(p0,p1,p2,cols = 1)
  # 
  
  return(df_plot)
  
  
}
###############################################################################
cross_corr_signals_v2 <- function(df,ref_sig,test_sig,time_range=c(-10:10)){
  
  # df deve avere 2 colonne e basta
  
  
  ### cross correlation
  cor.cross <- function(x0, y0, i=0) {
    #
    # Sample autocorrelation at (integral) lag `i`:
    # Positive `i` compares future values of `x` to present values of `y`';
    # negative `i` compares past values of `x` to present values of `y`.
    #
    if (i < 0) {x<-y0; y<-x0; i<- -i}
    else {x<-x0; y<-y0}
    n <- length(x)
    cor(x[(i+1):n], y[1:(n-i)], use="complete.obs")
  }
  ####
  df_1 <- df[,ref_sig]
  df_2 <- df[,test_sig]
  
  data.cor <- round(sapply(time_range, function(i) cor.cross(df_1, df_2, i)),2)
  
  df_plot <- as.data.frame(cbind(time_range,data.cor))
  colnames(df_plot) <- c("Lags","CrossCorrelation")
  
  best_x_corr_pos <- min(time_range[which(data.cor ==  max(data.cor)[1])])
  best_y_corr_pos <- max(data.cor)[1]
  
  best_x_corr_neg <- max(time_range[which(data.cor ==  min(data.cor)[1])])
  best_y_corr_neg <- min(data.cor)[1]
  
  
  gg1 <- ggplot(df_plot,aes(x=Lags, y=CrossCorrelation))+
    geom_line(color="orange4")+
    geom_point(size=3,color="orange")+
    # geom_point(x=best_x_corr_pos,y=best_y_corr_pos,color="purple",size=5) +
    # geom_point(x=best_x_corr_neg,y=best_y_corr_neg,color="purple",size=5) +
    # geom_point(x=0,y=df_plot[which(df_plot$Lags ==0),2],color="green",size=5) +
    geom_vline(xintercept = 0,linetype="dashed", color = "red") +
    geom_hline(yintercept = 0,linetype="dashed", color = "red")+
    ylim(-1,1)+
    ggtitle(paste0(ref_sig," vs. ",test_sig))
  
  
  return(df_plot)
}


###########################################################################
# calcola derivata
derivata_2pt_time_series <- function(time_series,add_last=TRUE){
  derivata <- c()
  x <- 1:length(time_series)
  for (i in 1:(length(time_series)-1)){
    derivata[i] <- (time_series[i+1]-time_series[i])/(x[i+1]-x[i])
    if(isTRUE(add_last)){
      derivata <- c(derivata,derivata[length(time_series)-1])
    }
  }
  return(derivata)
}
############################################################################

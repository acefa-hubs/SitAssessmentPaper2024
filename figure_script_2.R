setwd("C:/Users/EALESO/R Projects/Winter 2024 paper/ForGitHub")


############################################################################################################################################
## Loading required packages

library(ggplot2)
library(rstan)
library(patchwork)
library(cowplot)
library(gridExtra)
library(tidyverse)

############################################################################################################################################
## Loading required functions

source('R/ps_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')


############################################################################################################################################
## Loading and formatting the data

################################################################################
## Loading in the data 

origin_date <- as.Date("2024-09-12")
date_column <- "notification_date"

# COVID data
df_cov <- read.csv("data/cov_dummy.csv")

# Influenza data
df_inf <- read.csv("data/inf_dummy.csv")

# RSV data
df_rsv <- read.csv("data/rsv_dummy.csv")


################################################################################
## Formatting the data and extracting correct COVID data

df_inf <- format_inf_type(df_inf, df_infT, date_column)
df_cov <- df_cov[df_cov$test_type=="PCR",]

################################################################################
##  Set limits on dates to consider

# Covid date limits
max_date <- origin_date
min_date <- as.Date("2023-01-01")
df_cov <- df_cov[df_cov[,date_column]<=max_date & df_cov[,date_column]>=min_date,]

#influenza date limits
max_date <- origin_date
min_date <- as.Date("2022-01-01")
df_inf <- df_inf[df_inf[,date_column]<=max_date & df_inf[,date_column]>=min_date,]


# RSV date limits
max_date <- origin_date
min_date <- as.Date("2022-01-01")
df_rsv <- df_rsv[df_rsv[,date_column]<=max_date & df_rsv[,date_column]>=min_date,]

################################################################################
## Ensuring data is in correct order

df_cov[,date_column] <- as.Date(df_cov[,date_column])
df_inf[,date_column] <- as.Date(df_inf[,date_column])
df_rsv[,date_column] <- as.Date(df_rsv[,date_column])

df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1
df_inf$time_index <- as.numeric(df_inf[,date_column]) - min(as.numeric(df_inf[,date_column]))+1
df_rsv$time_index <- as.numeric(df_rsv[,date_column]) - min(as.numeric(df_rsv[,date_column]))+1


df_cov <- df_cov[order(df_cov$time_index),]
df_inf <- df_inf[order(df_inf$time_index),]
df_rsv <- df_rsv[order(df_rsv$time_index),]


#########################################################################################################################################
## Reading in the stan models and getting modelled estimates
#########################################################################################################################################

# Get dates that will match name of saved model outputs
max_dates_considered <- max(df_rsv$notification_date) - seq(0, 7*25, by=7)

################################################################################
## COVID Stan models
for(i in 1:length(max_dates_considered)){
  print(i)
  max_date <- max_dates_considered[i]
  min_date <- max_dates_considered[i]-365
  df_tmp <- df_cov[df_cov[,date_column]<=max_date & df_cov[,date_column]>min_date,]
  
  df_tmp$time_index <- as.numeric(df_tmp[,date_column]) - min(as.numeric(df_tmp[,date_column]))+1
  
  cov_fit <- readRDS(paste('fitted_stan_models/',max_date, '-cov_fit.rds', sep=""))
  
  
  tmp_mod_inc <- ps_single_incidence(cov_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column])
  tmp_mod_gr <- ps_single_growth_rate(cov_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column])
  
  tmp_mod_inc$max_date <- max_date
  tmp_mod_gr$max_date <- max_date
  
  tmp_mod_inc$col <- i%%2
  tmp_mod_gr$col <- i%%2
  
  if(i==1){
    cov_mod_inc <- tmp_mod_inc
    cov_mod_gr <- tmp_mod_gr
    
  } else{
    cov_mod_inc <- rbind(cov_mod_inc, tmp_mod_inc)
    cov_mod_gr <- rbind(cov_mod_gr, tmp_mod_gr)
  }
  
}


################################################################################
## RSV Stan models
for(i in 1:length(max_dates_considered)){
  print(i)
  max_date <- max_dates_considered[i]
  min_date <- max_dates_considered[i]-365
  df_tmp <- df_rsv[df_rsv[,date_column]<=max_date & df_rsv[,date_column]>min_date,]
  
  df_tmp$time_index <- as.numeric(df_tmp[,date_column]) - min(as.numeric(df_tmp[,date_column]))+1
  
  rsv_fit <- readRDS(paste('fitted_stan_models/',max_date, '-rsv_fit.rds', sep=""))
  
  
  tmp_mod_inc <- ps_single_incidence(rsv_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column])
  tmp_mod_gr <- ps_single_growth_rate(rsv_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column])
  
  tmp_mod_inc$max_date <- max_date
  tmp_mod_gr$max_date <- max_date
  
  tmp_mod_inc$col <- i%%2
  tmp_mod_gr$col <- i%%2
  
  if(i==1){
    rsv_mod_inc <- tmp_mod_inc
    rsv_mod_gr <- tmp_mod_gr
    
  } else{
    rsv_mod_inc <- rbind(rsv_mod_inc, tmp_mod_inc)
    rsv_mod_gr <- rbind(rsv_mod_gr, tmp_mod_gr)
  }
  
}


################################################################################
## Influenza Stan models
for(i in 1:length(max_dates_considered)){
  print(i)
  max_date <- max_dates_considered[i]
  min_date <- max_dates_considered[i]-365
  df_tmp <- df_inf[df_inf[,date_column]<=max_date & df_inf[,date_column]>min_date,]
  
  df_tmp$time_index <- as.numeric(df_tmp[,date_column]) - min(as.numeric(df_tmp[,date_column]))+1
  
  inf_fit <- readRDS(paste('fitted_stan_models/',max_date, '-inf_fit.rds', sep=""))
  
  tmp_mod_inc <- ps_incidence(inf_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column],
                              num_path = 3,
                              pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))
  
  tmp_mod_gr <- ps_growth_rate(inf_fit, df_tmp$time_index, num_days=nrow(df_tmp), time_labels = df_tmp[,date_column],
                               num_path = 3,
                               pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))
  
  tmp_mod_inc$max_date <- max_date
  tmp_mod_gr$max_date <- max_date
  
  tmp_mod_inc$col <- i%%2
  tmp_mod_gr$col <- i%%2
  
  if(i==1){
    inf_mod_inc <- tmp_mod_inc
    inf_mod_gr <- tmp_mod_gr
    
  } else{
    inf_mod_inc <- rbind(inf_mod_inc, tmp_mod_inc)
    inf_mod_gr <- rbind(inf_mod_gr, tmp_mod_gr)
  }
  
}

################################################################################
# Some formatting so that each model is only plotted for a single week (and other estimates masked)

##################################################
# COVID
cov_mod_inc$min_date <- cov_mod_inc$max_date-7
cov_mod_inc$mask <- 0.0
cov_mod_inc[cov_mod_inc$time>=cov_mod_inc$min_date &cov_mod_inc$time<=cov_mod_inc$max_date,]$mask <- 1.0

cov_mod_gr$min_date <- cov_mod_gr$max_date-7
cov_mod_gr$mask <- 0.0
cov_mod_gr[cov_mod_gr$time>=cov_mod_gr$min_date &cov_mod_gr$time<=cov_mod_gr$max_date,]$mask <- 1.0

cov_mod_gr$col <- as.factor(cov_mod_gr$col)
cov_mod_inc$col <- as.factor(cov_mod_inc$col)

##################################################
# RSV
rsv_mod_inc$min_date <- rsv_mod_inc$max_date-7
rsv_mod_inc$mask <- 0.0
rsv_mod_inc[rsv_mod_inc$time>=rsv_mod_inc$min_date &rsv_mod_inc$time<=rsv_mod_inc$max_date,]$mask <- 1.0

rsv_mod_gr$min_date <- rsv_mod_gr$max_date-7
rsv_mod_gr$mask <- 0.0
rsv_mod_gr[rsv_mod_gr$time>=rsv_mod_gr$min_date &rsv_mod_gr$time<=rsv_mod_gr$max_date,]$mask <- 1.0

rsv_mod_gr$col <- as.factor(rsv_mod_gr$col)
rsv_mod_inc$col <- as.factor(rsv_mod_inc$col)

##################################################
# Influenza
inf_mod_inc$min_date <- inf_mod_inc$max_date-7
inf_mod_inc$mask <- 0.0
inf_mod_inc[inf_mod_inc$time>=inf_mod_inc$min_date &inf_mod_inc$time<=inf_mod_inc$max_date,]$mask <- 1.0

inf_mod_gr$min_date <- inf_mod_gr$max_date-7
inf_mod_gr$mask <- 0.0
inf_mod_gr[inf_mod_gr$time>=inf_mod_gr$min_date &inf_mod_gr$time<=inf_mod_gr$max_date,]$mask <- 1.0

inf_mod_gr$col <- as.factor(inf_mod_gr$col)
inf_mod_inc$col <- as.factor(inf_mod_inc$col)

#################################################################################################################################################
## Making the figures
#################################################################################################################################################

################################################################################
# Figure 4
################################################################################

# Set first date for plot
first_date <- as.Date("2024-09-10") - 25*7

# Set colours for figure
cols <- c("blue4","red4")

############################################
## COVID panels
############################################

# Modelled cases
cov1 <- ggplot(cov_mod_inc[cov_mod_inc$mask==1,])+
  geom_line(aes(x=time, y=y, group=max_date, color=col))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=max_date, fill=col), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=max_date, fill=col), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_cov, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=cov_mod_inc[cov_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=ub_95),linetype="dashed")+
  geom_line(data=cov_mod_inc[cov_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=lb_95),linetype="dashed")+
  ylab("Cases")+
  xlab("Date")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,400))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(legend.position = "none")

# Modelled growth rate
cov2<- ggplot(cov_mod_gr[cov_mod_gr$mask==1,])+
  geom_line(aes(x=time, y=y, group=max_date, color=col))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=max_date, fill=col), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=max_date, fill=col), alpha=0.2)+
  geom_line(data=cov_mod_gr[cov_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=ub_95),linetype="dashed")+
  geom_line(data=cov_mod_gr[cov_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=lb_95),linetype="dashed")+
  theme_bw(base_size = 14)+
  ylab("Growth rate")+
  xlab("Date")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  geom_hline(yintercept = 0, linetype= "dotted")+
  coord_cartesian(xlim=c(first_date, origin_date-10))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(legend.position = "none")


############################################
## rsv panels
############################################

# Modelled cases
rsv1<-ggplot(rsv_mod_inc[rsv_mod_inc$mask==1,])+
  geom_line(aes(x=time, y=y, group=max_date, color=col))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=max_date, fill=col), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=max_date, fill=col), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_rsv, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=rsv_mod_inc[rsv_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=ub_95),linetype="dashed")+
  geom_line(data=rsv_mod_inc[rsv_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=lb_95),linetype="dashed")+
  ylab("Cases")+
  xlab("Date")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,400))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(legend.position = "none")

# Modelled growth rate
rsv2<-ggplot(rsv_mod_gr[rsv_mod_gr$mask==1,])+
  geom_line(aes(x=time, y=y, group=max_date, col=col))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=max_date, fill=col), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=max_date, fill=col), alpha=0.2)+
  geom_line(data=rsv_mod_gr[rsv_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=ub_95),linetype="dashed")+
  geom_line(data=rsv_mod_gr[rsv_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=lb_95),linetype="dashed")+
  theme_bw(base_size = 14)+
  ylab("Growth rate")+
  xlab("Date")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  geom_hline(yintercept = 0, linetype= "dotted")+
  coord_cartesian(xlim=c(first_date, origin_date-10))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
    theme(legend.position = "none")

# Reformatting the subplots a little
cov1 <- cov1+labs(tag="A")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))+
  annotate("label",label="SARS-CoV-2", y=Inf, x = as.Date("2024-03-16"), fill= "white", color="black", vjust=1.2, size=5, hjust=0)

cov2 <- cov2+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

rsv1 <- rsv1+labs(tag="B")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))+
  annotate("label",label="RSV", y=Inf, x = as.Date("2024-03-16"), fill= "white", color="black", vjust=1.2, size=5, hjust=0)
rsv2 <- rsv2

# Combine subplots and save
cov1 + cov2 + rsv1 +rsv2 + plot_layout(nrow=4, heights=c(2,1,2,1))

ggsave('figure/Figure4.png', width=8, height=10)
ggsave('figure/Figure4.pdf', width=8, height=10)


############################################################################################
# Figure 5
############################################################################################

# Set first date for plot
first_date <- as.Date("2024-09-10") - 25*7

# Set colours for plot
cols <- RColorBrewer::brewer.pal(8, "Paired")
cols <- c(cols[1:2], cols[5:8], "black","grey40")

# Plot to steal legend from
leg <- ggplot(inf_mod_inc[inf_mod_inc$mask==1 ,])+
  geom_line(aes(x=time, y=y, group=interaction(pathogen, max_date), color=pathogen ))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=interaction(pathogen, max_date), fill=pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  scale_color_manual("Influenza subtype",values=cols[c(2,4,6,8)])+
  scale_fill_manual("Influenza subtype",values=cols[c(2,4,6,8)])+
  theme(legend.position = c(0.2,0.7),
        legend.background = element_rect(color="black"))

# Extracting legend from above plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(leg)


# Some reformatting for facetting in later plot
inf_mod_inc$pathogen2 <- factor(inf_mod_inc$pathogen, levels = c("Total", "Influenza B", "Influenza A H1N1", "Influenza A H3N2"))
inf_mod_gr$pathogen2 <- factor(inf_mod_gr$pathogen, levels = c("Total", "Influenza B", "Influenza A H1N1", "Influenza A H3N2"))

# Modelled cases by subtype
inf1<-ggplot(inf_mod_inc[inf_mod_inc$mask==1 ,])+
  geom_line(aes(x=time, y=y, group=interaction(pathogen, max_date), color=interaction(col,pathogen) ))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=interaction(pathogen, max_date), fill=interaction(col,pathogen)), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=interaction(pathogen, max_date), fill=interaction(col,pathogen)), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_inf, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=inf_mod_inc[inf_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=ub_95, color=interaction(col,pathogen)),linetype="dashed")+
  geom_line(data=inf_mod_inc[inf_mod_inc$max_date==max_dates_considered[1],], aes(x=time, y=lb_95, color=interaction(col,pathogen)),linetype="dashed")+
  ylab("Cases")+
  xlab("Date")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  #scale_color_brewer(palette = "Paired")+
  #scale_fill_brewer(palette = "Paired")+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,850))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(legend.position = "none")

# Modelled growth rate by subtype
inf2<-ggplot(inf_mod_gr[inf_mod_gr$mask==1,])+
  geom_line(aes(x=time, y=y, group=interaction(pathogen, max_date), col=interaction(col,pathogen)))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50,group=interaction(pathogen, max_date), fill=interaction(col,pathogen)), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95,group=interaction(pathogen, max_date), fill=interaction(col,pathogen)), alpha=0.2)+
  geom_line(data=inf_mod_gr[inf_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=ub_95, color=interaction(col,pathogen)),linetype="dashed")+
  geom_line(data=inf_mod_gr[inf_mod_gr$max_date==max_dates_considered[1],], aes(x=time, y=lb_95, color=interaction(col,pathogen)),linetype="dashed")+
  theme_bw(base_size = 14)+
  ylab("Growth rate")+
  xlab("Date")+
  facet_wrap(.~pathogen2, scales= "free_y", nrow=4)+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  geom_hline(yintercept = 0.0, linetype="dotted")+
  coord_cartesian(xlim=c(first_date, origin_date-10))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(strip.text = element_blank(),
        legend.position = "none")


# Create data frame with weekly number of H3N2 samples vs H1N1 samples
df <- data.frame()
for(i in 1:length(max_dates_considered) ){
  df_tmp <- df_inf[df_inf$notification_date<=max_dates_considered[i] & df_inf$notification_date>max_dates_considered[i]-7,]
  row_df <- data.frame(date = max_dates_considered[i],
                       A_H1N1 = sum(df_tmp$A_H1N1),
                       A_H3N2 = sum(df_tmp$A_H3N2),
                       A_unk = sum(df_tmp$A_unk))
  df <- rbind(df, row_df)
}

df <- pivot_longer(df[c("date","A_H3N2", "A_H1N1")], cols=c("A_H3N2", "A_H1N1"))

# Plot stacked bar of influenza A subtype data
inf3 <- ggplot(df)+
  geom_col(aes(x=date-3.5, y =value, fill=name),color="black", linewidth=0.2, width=7)+
  theme_bw(base_size = 14)+
  scale_fill_manual(values=cols[c(1,3)])+
  #coord_cartesian(xlim=c(min_date+30, max_date-30))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  ylab("Modelled cases")+
  xlab("Date")+
  scale_y_continuous("Subtypes determined")+
  coord_cartesian(xlim=c(first_date, origin_date-10))+
  theme(legend.position = "none",
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill=NA))


# Adjust subplots a little
inf1 <- inf1+labs(tag="A")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95),
         legend.position = "none")+
  annotation_custom(legend, x = first_date+30, xmax = first_date+30, 
                    ymin = 600, ymax = 600)


inf2 <- inf2 + labs(tag="B")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95),
         legend.position = "none")

inf3 <- inf3 + labs(tag="C")+
  theme(plot.tag.position = c(0.01,0.9),
        panel.grid = element_blank(),
        axis.title.y = element_blank())+
  annotate("label",label="Influenza A subtype data", y=Inf, x = as.Date("2024-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)

# Merge subplots and save
inf1 +inf2 + inf3 + plot_layout(nrow=3, heights=c(1,2,0.25))

ggsave('figure/Figure5.png', width=8, height=11)
ggsave('figure/Figure5.pdf', width=8, height=11)

############################################################################################################################################
## Loading required packages

library(ggplot2)
library(rstan)
library(patchwork)
library(RColorBrewer)
library(dplyr)
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
## Formatting the data

df_cov <- df_cov[df_cov$test_type=="PCR",]

df_cov[,date_column] <- as.Date(df_cov[,date_column])
df_inf[,date_column] <- as.Date(df_inf[,date_column])
df_rsv[,date_column] <- as.Date(df_rsv[,date_column])

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

df_cov$time_index <- as.numeric(df_cov[,date_column]) - min(as.numeric(df_cov[,date_column]))+1
df_inf$time_index <- as.numeric(df_inf[,date_column]) - min(as.numeric(df_inf[,date_column]))+1
df_rsv$time_index <- as.numeric(df_rsv[,date_column]) - min(as.numeric(df_rsv[,date_column]))+1

df_cov <- df_cov[order(df_cov$time_index),]
df_inf <- df_inf[order(df_inf$time_index),]
df_rsv <- df_rsv[order(df_rsv$time_index),]


############################################################################################################################################
## Reading in the fitted stan models (fit to all overall data)

cov_fit  <- readRDS(paste('fitted_stan_models/', 'cov_fit-overall.rds', sep=""))
inf_fit  <- readRDS(paste('fitted_stan_models/', 'inf_fit-overall.rds', sep=""))
rsv_fit  <- readRDS(paste('fitted_stan_models/', 'rsv_fit-overall.rds', sep=""))


############################################################################################################################################
## Extract posterior model estimates for RSV

# Smoothed cases including day of the week effect
rsv_mod_inc_dow <- ps_single_incidence_dow(rsv_fit,
                                           week_effect = 7,
                                           DOW = (df_rsv$time_index %% 7)+1,
                                           X=df_rsv$time_index, num_days=nrow(df_rsv), time_labels = df_rsv[,date_column])

# Smoothed cases
rsv_mod_inc <- ps_single_incidence(rsv_fit, df_rsv$time_index, num_days=nrow(df_rsv), time_labels = df_rsv[,date_column])

# Growth rate
rsv_mod_gr <- ps_single_growth_rate(rsv_fit, df_rsv$time_index, num_days=nrow(df_rsv), time_labels = df_rsv[,date_column])


############################################################################################################################################
## Extract posterior model estimates for COVID

# Smoothed cases including day of the week effect
cov_mod_inc_dow <- ps_single_incidence_dow(cov_fit,
                                           week_effect = 7,
                                           DOW = (df_cov$time_index %% 7)+1,
                                           X=df_cov$time_index, num_days=nrow(df_cov), time_labels = df_cov[,date_column])

# Smoothed cases
cov_mod_inc <- ps_single_incidence(cov_fit, df_cov$time_index, num_days=nrow(df_cov), time_labels = df_cov[,date_column])

# Growth rate
cov_mod_gr <- ps_single_growth_rate(cov_fit, df_cov$time_index, num_days=nrow(df_cov), time_labels = df_cov[,date_column])


###############################################################################################
## Extract posterior model estimates for influenza

# Smoothed cases
inf_mod_inc <- ps_incidence(inf_fit, df_inf$time_index, num_days=nrow(df_inf), time_labels = df_inf[,date_column],
                            num_path = 3,
                            pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))

# Smoothed cases including day of the week effect
inf_mod_inc_dow <- ps_incidence_dow(inf_fit, X=df_inf$time_index, num_days=nrow(df_inf), time_labels = df_inf[,date_column],
                                    week_effect = 7,
                                    DOW = (df_inf$time_index %% 7)+1,
                                    num_path = 3,
                                    pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))

# Growth rate
inf_mod_gr <- ps_growth_rate(inf_fit, df_inf$time_index, num_days=nrow(df_inf), time_labels = df_inf[,date_column],
                             num_path = 3,
                             pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))


# Proportion each subtype
inf_mod_prop <- ps_proportion(inf_fit, df_inf$time_index, num_days=nrow(df_inf), time_labels = df_inf[,date_column],
                              num_path = 3,
                              comb_num=list(c(1,2), c(3)),
                              comb_den=list(c(1,2,3), c(1,2,3)),
                              comb_names=c("Influenza A", "Influenza B"))

#########################################################################################################
# Plotting figures
#########################################################################################################
first_date <- as.Date("2024-09-10") - 27*7

###########################################################
## Figure 1
###########################################################

# COVID cases panel
cov1 <- ggplot(cov_mod_inc)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_cov, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=df_cov, aes(x=notification_date, y=cases) , linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,400))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")

# COVID growth rate panel
cov2 <- ggplot(cov_mod_gr)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(-0.052,0.052))+
  scale_y_continuous(breaks=c(-0.04,0.0,0.04))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")

# RSV cases panel
rsv1 <- ggplot(rsv_mod_inc)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_rsv, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=df_rsv, aes(x=notification_date, y=cases) , linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,400))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")

# RSV growth rate panel
rsv2 <- ggplot(rsv_mod_gr)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(-0.05,0.05))+
  scale_y_continuous(breaks=c(-0.04,0.0,0.04))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")

# Influenza cases panel
inf1 <- ggplot(inf_mod_inc[inf_mod_inc$pathogen=="Total",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_inf, aes(x=notification_date, y=cases), size=0.8)+
  geom_line(data=df_inf, aes(x=notification_date, y=cases) , linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,850))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")

# Influenza growth rate panel
inf2 <- ggplot(inf_mod_gr[inf_mod_gr$pathogen=="Total",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(-0.07,0.07))+
  #scale_y_continuous(breaks=c(-0.04,0.0,0.04))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")



## Formatting the panels and combining into single figure

cov1 <- cov1+labs(tag="A")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))+
  annotate("label",label="SARS-CoV-2", y=Inf, x = as.Date("2024-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)

cov2 <- cov2+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

rsv1 <- rsv1+labs(tag="B")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))+
  annotate("label",label="RSV", y=Inf, x = as.Date("2024-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)
rsv2 <- rsv2+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

inf1 <- inf1+labs(tag="C")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))+
  annotate("label",label="Influenza", y=Inf, x = as.Date("2024-03-15"), fill= "white", color="black", vjust=1.2, size=5, hjust=0.0)
inf2 <- inf2


cov1 + cov2 + rsv1 +rsv2 +inf1 + inf2+ plot_layout(nrow=6, heights=c(2,1,2,1,2,1))

ggsave('figure/Figure1.png', width=6, height=12)
ggsave('figure/Figure1.pdf', width=6, height=12)



###########################################################
## Figure 2
###########################################################

# Setting colours for use in figure
cols <- c("navy", "red4", "gold3", "black")

# Plot cases by subtype
inc <- ggplot(inf_mod_inc)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  #geom_point(data=df_inf, aes(x=notification_date, y=cases), size=0.8)+
  #geom_line(data=df_inf, aes(x=notification_date, y=cases) , linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  scale_color_manual("Influenza subtype",values=cols)+
  scale_fill_manual("Influenza subtype",values=cols)+
  coord_cartesian(xlim=c(first_date, origin_date-10), ylim=c(0,650))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(legend.position = c(0.2,0.7),
        legend.background=element_rect(color="black"))

# Plot growth rate by subtype
gr <- ggplot(inf_mod_gr[inf_mod_gr$pathogen!="Total" & inf_mod_gr$time>first_date-50,])+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  scale_color_manual("Influenza subtype",values=cols)+
  scale_fill_manual("Influenza subtype",values=cols)+
  facet_wrap(.~pathogen, nrow=3)+
  coord_cartesian(xlim=c(first_date, origin_date-10))+
  scale_y_continuous(breaks=c(-0.1,0.0,0.1))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme(strip.text = element_blank(),
        legend.position = "none")


# Format and plot panels as single figure
inc <- inc+labs(tag="A")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))
gr <- gr + labs(tag="B")+
  theme(plot.tag.position = c(0.01,0.95))


inc+gr +plot_layout(nrow=2, heights=c(1,1.5))

ggsave('figure/Figure2.png', width=6, height=8)
ggsave('figure/Figure2.pdf', width=6, height=8)



###########################################################
## Figure 3
###########################################################

# Set colours and alpha for figures
colours <- c(rev(RColorBrewer::brewer.pal(7,"Dark2")), "black")
alpha= c(rep(0.2, 7), 0.4)

# Format dates to all appear in same year (so multiple years can be compared)
cov_new <- cov_mod_inc
cov_new$year <- format(cov_new$time, format="%Y")
cov_new$time2 <- cov_new$time + (2024 - as.numeric(cov_new$year))*365

rsv_new <- rsv_mod_inc
rsv_new$year <- format(rsv_new$time, format="%Y")
rsv_new$time2 <- rsv_new$time + (2024 - as.numeric(rsv_new$year))*365

inf_new <- inf_mod_inc
inf_new$year <- format(inf_new$time, format="%Y")
inf_new$time2 <- inf_new$time + (2024 - as.numeric(inf_new$year))*365

cov_new$pathogen <- "SARS-CoV-2"
rsv_new$pathogen <- "RSV"
inf_new[inf_new$pathogen=="Total",]$pathogen <- "Influenza (all)"
path_new <- rbind(cov_new, rsv_new, inf_new)

path_new$pathogen <- factor(path_new$pathogen, levels = c("SARS-CoV-2", "RSV", "Influenza (all)",
                                                          "Influenza A H3N2", "Influenza A H1N1", "Influenza B") )

ggplot(path_new)+
  geom_line(aes(x=time2, y=y, group=year, color=year))+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time2, y=y, ymin=lb_95, ymax=ub_95, group=year, fill=year, alpha=year))+
  theme_bw(base_size = 14)+
  ylab("Cases")+
  xlab("Date")+
  facet_wrap(.~pathogen, scales = "free_y")+
  scale_color_manual("Season",values = colours, breaks=c("2015","2016","2017","2018","2019","2022","2023","2024"))+
  scale_fill_manual("Season",values = colours, breaks=c("2015","2016","2017","2018","2019","2022","2023","2024"))+
  scale_alpha_manual("Season",values=alpha, breaks=c("2015","2016","2017","2018","2019","2022","2023","2024"))+
  coord_cartesian(xlim=c(as.Date("2024-02-01"), origin_date-10))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b")+
  theme(legend.position = c(0.76,0.24),
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill=NA))


ggsave('figure/Figure3.png', width=10, height=6)
ggsave('figure/Figure3.pdf', width=10, height=6)


################################################################################
## Supplementary Figures
################################################################################

###########################################################
# RSV
###########################################################

# Set limits for the plot
min_date_plot <- as.Date("2022-02-01")+30
max_date_plot <- as.Date("2024-09-01")-30

# Plot smoothed cases with day of the week effect
rsv0 <- ggplot(rsv_mod_inc_dow)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_rsv, aes(x=notification_date, y=cases), size=0.2)+
  #geom_line(data=df_rsv, aes(x=notification_date, y=cases), linewidth=0.1)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

# Plot smoothed cases
rsv1 <- ggplot(rsv_mod_inc)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_rsv, aes(x=notification_date, y=cases), size=0.2)+
  #geom_line(data=df_rsv, aes(x=notification_date, y=cases), linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

# Plot growth rate
rsv2 <- ggplot(rsv_mod_gr)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  theme()

# Combine panels and save
rsv0 + rsv1 +rsv2 +plot_layout(nrow=3)
ggsave('figure/SFigure3.png', width=20, height=10)
ggsave('figure/SFigure3.pdf', width=20, height=10)



###########################################################
# COVID
###########################################################

# Set date limits for plot
min_date_plot <- as.Date("2023-01-01")+25
max_date_plot <- as.Date("2024-09-01")-15

# Plot smoothed cases with day of the week effect
cov0 <- ggplot(cov_mod_inc_dow)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_cov, aes(x=notification_date, y=cases), size=0.2)+
  #geom_line(data=df_cov, aes(x=notification_date, y=cases), linewidth=0.1)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

# Plot smoothed cases 
cov1 <- ggplot(cov_mod_inc)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_cov, aes(x=notification_date, y=cases), size=0.2)+
  #geom_line(data=df_cov, aes(x=notification_date, y=cases), linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

# Plot growth rate
cov2 <- ggplot(cov_mod_gr)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  theme()

# Combine panels and save
cov0 + cov1 +cov2 +plot_layout(nrow=3)
ggsave('figure/SFigure2.png', width=20, height=10)
ggsave('figure/SFigure2.pdf', width=20, height=10)


###########################################################
# Influenza
###########################################################

# Set date limits
min_date_plot <- as.Date("2022-01-01")+40
max_date_plot <- as.Date("2024-09-01")-30

# Plot smoothed cases with day of the week effect
inf0 <- ggplot(inf_mod_inc_dow[inf_mod_inc_dow$pathogen=="Total",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_inf, aes(x=notification_date, y=cases), size=0.2)+
  #geom_line(data=df_inf, aes(x=notification_date, y=cases), linewidth=0.1)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

# Plot smoothed cases
inf1 <- ggplot(inf_mod_inc[inf_mod_inc$pathogen=="Total",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_point(data=df_inf, aes(x=notification_date, y=cases), size=0.2)+
  #geom_line(data=df_inf, aes(x=notification_date, y=cases), linewidth=0.2)+
  ylab("Cases")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank())

# Plot growth rate
inf2 <- ggplot(inf_mod_gr[inf_mod_gr$pathogen=="Total",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  theme()

# Combine panels and save
inf0 + inf1 +inf2 +plot_layout(nrow=3)
ggsave('figure/SFigure4.png', width=20, height=10)
ggsave('figure/SFigure4.pdf', width=20, height=10)



###########################################################
# Influenza subtypes
###########################################################

# Set date limits for plot
min_date_plot <- as.Date("2022-01-01")+40
max_date_plot <- as.Date("2024-09-01")-30

# Plot smoothed cases by subtype
infS1 <- ggplot(inf_mod_inc[inf_mod_inc$pathogen!="Total",])+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill = pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill = pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  ylab("Cases")+
  xlab("Date")+
  scale_color_manual("Influenza subtype",values = cols)+
  scale_fill_manual("Influenza subtype",values = cols)+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         legend.position = c(0.75,0.75),
         legend.background = element_rect(color="black"))

# Plot growth rate by subtype
infS2 <- ggplot(inf_mod_gr[inf_mod_gr$pathogen!="Total",])+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylab("Growth rate")+
  xlab("Date")+
  scale_color_manual("Influenza subtype",values = cols)+
  scale_fill_manual("Influenza subtype",values = cols)+
  facet_wrap(.~pathogen, nrow=3, scales = "free_y")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  theme(strip.text=element_blank(),
        legend.position = "none")

# Combine panels and save
infS1 + infS2 +plot_layout(nrow=2, heights = c(1,2))
ggsave('figure/SFigure5.png', width=16, height=10)
ggsave('figure/SFigure5.pdf', width=16, height=10)


################################################################################
## Influenza subtype proportion
################################################################################

# Set date limits for plot
min_date_plot <- as.Date("2022-01-01")+40
max_date_plot <- as.Date("2024-09-01")-30

# Set colours for plot
cols <- c("navy", "red4", "gold3","green4", "black")

# Get modelled proporion A vs B
inf_mod_prop1 <- ps_proportion(inf_fit, df_inf$time_index, num_days=nrow(df_inf), time_labels = df_inf[,date_column],
                               num_path = 3,
                               comb_num=list(c(1,2), c(3)),
                               comb_den=list(c(1,2,3), c(1,2,3)),
                               comb_names=c("Influenza A", "Influenza B"))

# Get modelled proportion H3N2 vs H1N1
inf_mod_prop2 <- ps_proportion(inf_fit, df_inf$time_index, num_days=nrow(df_inf), time_labels = df_inf[,date_column],
                               num_path = 3,
                               comb_num=list(c(1), c(2)),
                               comb_den=list(c(1,2), c(1,2)),
                               comb_names=c("H3N2", "H1N1"))


# Get weekly total of cases by subtype 
infA_plot <- df_inf
infA_plot$week <- 3.5+infA_plot$notification_date - as.numeric(infA_plot$notification_date)%%7

# Calculate sum of cases/subtyping, grouped by week
infA_plot <- infA_plot %>%
  group_by(week) %>%
  summarize(Atot = sum(A_H1N1 + A_H3N2 + A_unk),
            cases = sum(cases),
            A_H3N2 = sum(A_H3N2),
            A_H1N1 = sum(A_H1N1),
            A_unk = sum(A_unk),
            B = sum(B),
            AB = sum(AB))

infA_plot$p <- infA_plot$Atot/(infA_plot$Atot + infA_plot$B)
infA_plot$H3 <- infA_plot$A_H3N2/(infA_plot$A_H1N1+infA_plot$A_H3N2)

# Loop through all weeks and get 95% confidence intervals in the raw proportions (A vs B, and H3N2 vs H1N1)
for(i in 1:nrow(infA_plot)){
  if(infA_plot$Atot[i] + infA_plot$B[i]==0){
    infA_plot$p[i] <- NA
    infA_plot$p_lb[i] <- NA
    infA_plot$p_ub[i] <- NA
    
    
  } else{
    prop <- prop.test(infA_plot$Atot[i],infA_plot$Atot[i]+infA_plot$B[i])
    infA_plot$p_lb[i] <- prop$conf.int[1]
    infA_plot$p_ub[i] <- prop$conf.int[2]
    
    
  }
  
  if( (infA_plot$A_H1N1[i]+infA_plot$A_H3N2[i]) ==0){
    infA_plot$H3[i] <- NA
    infA_plot$H3_lb[i] <- NA
    infA_plot$H3_ub[i] <- NA
  } else{
    prop <- prop.test(infA_plot$A_H3N2[i],infA_plot$A_H3N2[i]+infA_plot$A_H1N1[i])
    infA_plot$H3_lb[i] <- prop$conf.int[1]
    infA_plot$H3_ub[i] <- prop$conf.int[2]
  }
  
  
}



# Plot proportion influenza A (relative to B)
sub1 <- ggplot(inf_mod_prop1[inf_mod_prop1$pathogen=="Influenza A",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  geom_point(data=infA_plot, aes(x=week, y=p), size=0.2)+
  geom_errorbar(data=infA_plot, aes(x=week, y=p, ymin=p_lb, ymax=p_ub), width=0, size=0.2)+
  theme_bw(base_size = 14)+
  ylab("Proportion influenza A\n(relative to influenza B)")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  theme(strip.text=element_blank(),
        legend.position = "none")

# Plot proportion influenza A H3N2 (relative to A H1N1)
sub2 <- ggplot(inf_mod_prop2[inf_mod_prop2$pathogen=="H3N2",])+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  geom_point(data=infA_plot, aes(x=week, y=H3), size=0.2)+
  geom_errorbar(data=infA_plot, aes(x=week, y=H3, ymin=H3_lb, ymax=H3_ub), width=0, size=0.2)+
  theme_bw(base_size = 14)+
  ylab("Proportion H3N2\n(relative to H1N1)")+
  xlab("Date")+
  coord_cartesian(xlim=c(min_date_plot, max_date_plot))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  theme(strip.text=element_blank(),
        legend.position = "none")



# Reformat weekly data for plotting by subtype below
inf_new <- pivot_longer(infA_plot[c("week","A_H3N2", "A_H1N1", "B", "A_unk", "AB")], cols=c("A_H3N2", "A_H1N1", "B", "A_unk", "AB"))
inf_new$name <- factor(inf_new$name)
inf_new$name <- factor(inf_new$name, labels = c("Influenza A H1N1", "Influenza A H3N2", "Influenza A (not subtyped)","Influenza any (not typed)", "Influenza B"))

# Plot raw weekly case data by subtype
sub3 <- ggplot(inf_new)+
  geom_col(aes(x=week, y =value, fill = name))+
  scale_fill_brewer("Influenza subtype",palette = "Dark2")+
  theme_bw(base_size = 14)+
  coord_cartesian(xlim=c(min_date+30, max_date-30))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%b\n%Y")+
  ylab("Modelled cases")+
  xlab("Date")+
  scale_y_continuous("Cases")+
  theme(legend.position = c(0.4,0.6),
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill=NA))

# Plot raw weekly case data for H3N1 and H1N1 only
sub4 <- ggplot(inf_new[inf_new$name %in% c("Influenza A H1N1", "Influenza A H3N2"),])+
  geom_col(aes(x=week, y =value, fill = name))+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw(base_size = 14)+
  coord_cartesian(xlim=c(min_date+30, max_date-30))+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  ylab("Modelled cases")+
  xlab("Date")+
  scale_y_continuous("Cases")+
  theme(legend.position = "none",
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill=NA))


# Reformat all subplots
sub1 <- sub1+labs(tag="A")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))

sub2 <- sub2+labs(tag="C")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))

sub3 <- sub3+labs(tag="B")+
  theme( axis.text.x=element_blank(),
         axis.title.x = element_blank(),
         plot.tag.position = c(0.01,0.95))

sub4 <- sub4+labs(tag="D")+
  theme(plot.tag.position = c(0.01,0.95))

# Combine subplots and save
sub1 + sub3 + sub2 +sub4+ plot_layout(nrow=4)
ggsave('figure/SFigure1.png', width=16, height=12)
ggsave('figure/SFigure1.pdf', width=16, height=12)


################################################################################ 
## The below code is used to generate dummy data for each pathogen
################################################################################ 

# Dummy data is randomly sampled using the median of posterior model estimates 
# for the smoothed trends in total cases and pathogens (for influenza)

# Daily case numbers are randomly sampled from a negative binomial distribution
# Daily number of influenza types/subtypes are randomly sampled from binomial 
# distributions for A/B and H3N2/H1N1

# Function for generating the dummy data for RSV and COVID
generate_synthetic_data_single <- function(fit, df){
  mod_inc_dow <- ps_single_incidence_dow(fit,
                                             week_effect = 7,
                                             DOW = (df$time_index %% 7)+1,
                                             X=df$time_index, num_days=nrow(df), time_labels = df[,date_column])
  
  post <- rstan::extract(fit)
  
  mu = mod_inc_dow$y
  sigma = mod_inc_dow$y*(1+(1/mean(post$phi)) )
  df$cases <- rnbinom(length(mod_inc_dow$y), size=mu**2/(sigma - mu), prob=mu/sigma )
  return(df)
}

# Function for generating the dummy data for influenza
generate_synthetic_data_inf <- function(fit, df){
  mod_inc_dow <- ps_incidence_dow(fit,
                                  week_effect = 7,
                                  DOW = (df$time_index %% 7)+1,
                                  X=df$time_index, num_days=nrow(df), 
                                  time_labels = df[,date_column],
                                  num_path = 3,
                                  pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))
  
  post <- rstan::extract(fit)
  
  mu = mod_inc_dow[mod_inc_dow$pathogen=="Total",]$y
  sigma = mod_inc_dow[mod_inc_dow$pathogen=="Total",]$y*(1+(1/mean(post$phi)) )
  df$cases <- rnbinom(length(mod_inc_dow[mod_inc_dow$pathogen=="Total",]$y), size=mu**2/(sigma - mu), prob=mu/sigma )
  
  inf_mod_prop <- ps_proportion(fit, df$time_index, num_days=nrow(df), time_labels = df[,date_column],
                                num_path = 3,
                                comb_num=list(c(1), c(1,2), c(3)),
                                comb_den=list(c(1,2), c(1,2,3), c(1,2,3)),
                                comb_names=c("H3N2","Influenza A", "Influenza B"))
  total_samples <- rowSums(df[4:8])
  probA <- inf_mod_prop[inf_mod_prop$pathogen=="Influenza A",]$y
  numA <- rbinom(length(total_samples), size=total_samples,probA)
  
  totalH3H1 <- rowSums(df[5:6])
  probH3 <- inf_mod_prop[inf_mod_prop$pathogen=="H3N2",]$y
  numH3 <-rbinom(length(totalH3H1), size=totalH3H1,probH3)
  
  df$A_H3N2 <- numH3
  df$A_H1N1 <- totalH3H1  - numH3
  df$A_unk <- numA - totalH3H1
  df$B <- total_samples - numA
  df$AB <- 0
  
  return(df)
}


# Generate the synthetic data sets for each pathogen
cov_dummy <- generate_synthetic_data_single(fit=cov_fit, df=df_cov)
rsv_dummy <- generate_synthetic_data_single(fit=rsv_fit, df=df_rsv)
inf_dummy <- generate_synthetic_data_inf(fit=inf_fit, df=df_inf)

# Write the dummy data to csv files
write.csv(cov_dummy, "cov_dummy_example.csv")
write.csv(rsv_dummy, "rsv_dummy_example.csv")
write.csv(inf_dummy, "inf_dummy_example.csv")

# Visualise difference between data used to fit model and data randomly sampled from model fits below

plot(cov_dummy$notification_date, cov_dummy$cases, type='l')
points(df_cov$notification_date, df_cov$cases, type='l', col='red')

plot(rsv_dummy$notification_date, rsv_dummy$cases, type='l')
points(df_rsv$notification_date, df_rsv$cases, type='l', col='red')


plot(inf_dummy$notification_date, inf_dummy$cases, type='l')
points(df_inf$notification_date, df_inf$cases, type='l', col='red')

plot(inf_dummy$notification_date, inf_dummy$A_unk, type='l')
points(df_inf$notification_date, df_inf$A_unk, type='l', col='red')

plot(inf_dummy$notification_date, inf_dummy$B, type='l')
points(df_inf$notification_date, df_inf$B, type='l', col='red')

plot(inf_dummy$notification_date, inf_dummy$A_H3N2, type='l')
points(df_inf$notification_date, df_inf$A_H3N2, type='l', col='red')

plot(inf_dummy$notification_date, inf_dummy$A_H1N1, type='l')
points(df_inf$notification_date, df_inf$A_H1N1, type='l', col='red')



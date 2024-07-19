###********************************####
####<<METADATA>>####
###********************************####

#Manuscript title:
#Groundtruthing assessments of lab-based coral thermal tolerance with field-based photogrammetry


#Authors:
#Crystal J. McRae12*^, Nathaniel Hanna Holloway3^, Guanyan Keelung Chen1,4, Michael T. Connelly5, 
#Hung-Kai Chen1, Zong-Min Ye1, Kendall S. Chancellor3, Yu-Chi Chang1, Ming-Qi Jiang1, Kwok Wai Lam1, 
#Yu-Ting Qiu1, Tung-Yung Fan1,2*, Stuart A. Sandin3

#1 Department of Planning and Research, National Museum of Marine Biology and Aquarium, Checheng, Pingtung, Taiwan, 944401.
#2 Department of Marine Biotechnology and Resources, National Sun Yat-sen University, Kaohsiung, Taiwan, 804201. 
#3 Scripps Institution of Oceanography, University of California San Diego, La Jolla, CA 92093.
#4 Hawaii Institute of Marine Biology, University of Hawaii, Honolulu, Hawaii, United States of America, 96822.
#5 Department of Invertebrate Zoology, Smithsonian Institution, Washington, DC, United States of America, 20560.

#^ Contributed equally
#* Corresponding authors: Crystal J. McRae and Tung-Yung Fan

#Contact Information: 
#Address: No. 2, Houwan Road, Houwan Village, Checheng Township, Pingtung County, Taiwan 944401.
#Phone: (886) 08-8825001 ext. 2248
#Email: crystal.j.mcrae@gmail.com (CJM); tyfan@nmmba.gov.tw (TYF)



###********************************###

#Colony collection: March 16th 2022
#Collection sites: Haikou (HK), Haikou reef flat (HKS), Inlet (IT), and Outlet (OT)
#n=6 colonies / site
#Colony fragmentation: March 25th 2022
#Pre-experiment measurements (photos & Fv/Fm): March 30th 2022
#CBASS experiment start: Noon on March 31st 2022


####********************************###

#clear workspace
rm(list=ls())

#set language to English
Sys.setenv(LANG = "en")

####********************************###



#####********************************####
#####REEF SITE MONITORING####
#####********************************####

####< INSITU TEMPERATURE>####

####1) Libraries####

library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(patchwork)


####2) Data####

#In-situ temp data (HKS data removed after Nov 11th due to suspected air exposure due to re-positioning of logger peg)
one_year_insitu_temp_HKS_halfyear <- read_csv("one_year_insitu_temp_HKS_halfyear.csv")
CBASS_insitu_temp <- one_year_insitu_temp_HKS_halfyear 

#remove UTC
CBASS_insitu_temp$DateTime <-  sub(" UTC", "", CBASS_insitu_temp$DateTime)
CBASS_insitu_temp$timeline <-  sub(" UTC", "", CBASS_insitu_temp$timeline)

#class check & organize the data
summary(CBASS_insitu_temp)

#Separate month, day, hour, and minute
TEMP_DATA <-
  CBASS_insitu_temp %>%
  mutate(year= year(DateTime), 
         month = month(DateTime),
         day = day(DateTime), 
         hour = hour(DateTime), 
         minute = minute(DateTime))

#Define seasons
INSITU_TEMP_DATA<-
  TEMP_DATA %>%
  mutate(month_conversion = ifelse(year== 2022, -2, 10))%>%
  mutate(MONTHID = month + month_conversion)%>%
  mutate(
    season = case_when(
      MONTHID %in% 10:12 ~ "winter",
      MONTHID %in%  4:6  ~ "summer",
      MONTHID %in%  7:9  ~ "fall",
      TRUE ~ "spring"))


#create separate date column for plotting
INSITU_TEMP_DATA <- tidyr::separate(INSITU_TEMP_DATA, 'DateTime',
                                    into = c('date', 'time'),
                                    sep= ' ')


INSITU_TEMP_DATA$date <- as.Date (INSITU_TEMP_DATA$date,"%Y-%m-%d")
INSITU_TEMP_DATA$site <- as.factor (INSITU_TEMP_DATA$site)
INSITU_TEMP_DATA$month <- as.integer (INSITU_TEMP_DATA$month)
INSITU_TEMP_DATA$MONTHID <- as.integer (INSITU_TEMP_DATA$MONTHID)
INSITU_TEMP_DATA$season <- as.factor (INSITU_TEMP_DATA$season)
INSITU_TEMP_DATA$timeline <- as.Date (INSITU_TEMP_DATA$timeline,"%Y-%m-%d  %H:%M:%S")

str(INSITU_TEMP_DATA)
summary(INSITU_TEMP_DATA)
head(INSITU_TEMP_DATA)


#Monthly summaries

#remove NAs from dataframe
INSITU_TEMP_DATA <- na.omit(INSITU_TEMP_DATA)

month_sum <-
  INSITU_TEMP_DATA%>%
  group_by(month, site, year)%>%
  summarise(TEMP_mean = mean(temp), TEMP_sd =sd(temp))

month_var <-
  INSITU_TEMP_DATA%>%
  group_by(month, site, year)%>%
  summarise(TEMP_max = max(temp), TEMP_min =min(temp))


#Daily summaries
day_sum <-
  INSITU_TEMP_DATA%>%
  group_by(date, site, month, year, season, MONTHID)%>%
  summarise(TEMP_mean = mean(temp), TEMP_sd =sd(temp), TEMP_max=max(temp), TEMP_min = min(temp))%>%
  mutate(day_range = TEMP_max - TEMP_min)

#Average of daily values
day_sum_short<-
  day_sum%>%
  group_by(month, year, site, season, MONTHID)%>%
  summarise(day_range = mean(day_range), meanday_max =mean (TEMP_max), meanday_min = mean(TEMP_min), meanday = mean(TEMP_mean), meanday_sd = mean(TEMP_sd))

#Average range and SD
range_sum <-
  day_sum_short%>%
  group_by(site)%>%
  summarise(mean_range=mean(day_range), sd_range=sd(day_range))

range_sum

#A tibble: 4 x 3
#site  mean_range sd_range
#<fct>      <dbl>    <dbl>
#1 HK       1.39    0.203
#2 HKS      1.78    0.341
#3 IT       2.17    1.16 
#4 OT       3.17    1.14


####3) Plots####

#first make subsets for each site

#HK
HK_TEMP <-
  INSITU_TEMP_DATA%>%
  filter(site =="HK")

summary(HK_TEMP)
head(HK_TEMP)
tail(HK_TEMP)


#HKS
HKS_TEMP <-
  INSITU_TEMP_DATA%>%
  filter(site =="HKS")

summary(HKS_TEMP)
head(HKS_TEMP)
tail(HKS_TEMP)


#OT
OT_TEMP <-
  INSITU_TEMP_DATA%>%
  filter(site =="OT")

summary(OT_TEMP)
head(OT_TEMP)
tail(OT_TEMP)


#IT
IT_TEMP <-
  INSITU_TEMP_DATA%>%
  filter(site =="IT")

summary(IT_TEMP)
head(IT_TEMP)
tail(IT_TEMP)


#dashed line in plots shows the average MMM based on NOAA data from 1985-2023

####*HK####

HK_plot <- ggplot(data=HK_TEMP, 
                  aes(x=timeline, y=temp)) +
  geom_line(size=0.05, colour='seagreen')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=28.5874, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(20, 35), breaks=c(20,22, 24, 26, 28, 30, 32, 34))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Haikou reef", y="Temperature (°C)", x="Month")

HK_plot 


####*HKS####

HKS_plot <- ggplot(data=HKS_TEMP, 
                   aes(x=date, y=temp)) +
  geom_line(size=0.05, colour='lightgreen')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=28.5874, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(20, 35), breaks=c(20, 22, 24, 26, 28, 30, 32, 34))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Haikou reef flat", y="Temperature (°C)", x="Month")

HKS_plot 



####*OT####

OT_plot <- ggplot(data=OT_TEMP, 
                  aes(x=date, y=temp)) +
  geom_line(size=0.15, colour='red')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=28.5874, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(20, 35), breaks=c(20,22, 24, 26, 28, 30, 32, 34))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Outlet reef", y="Temperature (°C)", x="Month")

OT_plot 


####*IT####

IT_plot <- ggplot(data=IT_TEMP, 
                  aes(x=date, y=temp)) +
  geom_line(size=0.15, colour='steelblue')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=28.5874, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(20, 35), breaks=c(20, 22, 24, 26, 28, 30, 32, 34))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Inlet reef", y="Temperature (°C)", x="Month")

IT_plot 


raw_temp_combo_plot <- HK_plot / HKS_plot / IT_plot / OT_plot
raw_temp_combo_plot

####*daily mean####

daymean_plot <- 
  ggplot(data=day_sum, aes(x=date, y=TEMP_mean, colour=site)) +
  geom_smooth(aes(fill=site))+
  theme_bw()+
  scale_color_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle=0, margin = margin(t=5, r=50)))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month")+
  scale_y_continuous(limits=c(21.25, 32), breaks=c(22, 24, 26, 28, 30, 32))+
  labs(y="Temperature (?C)", x="", title = "Mean Temperature")+
  theme(text=element_text(size=14,  family="sans"))

daymean_plot 

tempA <- daymean_plot +   theme(legend.position = "none")

daymean_plot +   theme(legend.position = "top")




#####*daily range####

daymean_range_plot <- 
  ggplot(data=day_sum, aes(x=date, y=day_range, colour=site)) +
  geom_smooth(aes(fill=site))+
  theme_bw()+
  scale_color_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle=0, margin = margin(t=5, r=50)))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month")+
  scale_y_continuous(limits=c(0, 5), breaks=c(0,1,2,3,4, 5))+
  labs(title="Daily Temperature Range", y="Temperature (?C)", x="")+
  theme(text=element_text(size=14,  family="sans"))

daymean_range_plot 

tempB <- daymean_range_plot  + theme(legend.position = "none")

tempB


#####*daily_max####

daymean_max_plot <- 
  ggplot(data=day_sum, aes(x=date, y=TEMP_max, colour=site)) +
  geom_smooth(aes(fill=site))+
  theme_bw()+
  scale_color_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle=0, margin = margin(t=5, r=50)))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month")+
  scale_y_continuous(limits=c(21.25, 32), breaks=c(22, 24, 26, 28, 30, 32))+
  labs(y="Temperature (?C)", x="", title = "Maximum Temperature")+
  theme(text=element_text(size=14,  family="sans"))

daymean_max_plot

tempC <- daymean_max_plot  + theme(legend.position = "none")

tempC 


####*daily_min####

daymean_min_plot <- 
  ggplot(data=day_sum, aes(x=date, y=TEMP_min, colour=site)) +
  geom_smooth(aes(fill=site))+
  theme_bw()+
  scale_color_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red")) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle=0, margin = margin(t=5, r=50)))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month")+
  scale_y_continuous(limits=c(21.25, 32), breaks=c(22, 24, 26, 28, 30, 32))+
  labs(y="Temperature (?C)", x="", title = "Minimum Temperature")+
  theme(text=element_text(size=14,  family="sans"))

daymean_min_plot

tempD <- daymean_min_plot  + theme(legend.position = "none")

tempD 

#COMBO plot
combo <- tempA + tempB + tempC + tempD + 
  plot_layout(ncol = 2)

combo 


####4) Analysis####

####*by site: mean sd####

#HK
HK.meansd<-
  HK_TEMP%>%
  summarise(mean.HK = mean(temp), sd.HK = sd(temp))

HK.meansd
#mean.HK    sd.HK
# <dbl>     <dbl>
# 27.4      2.53

#HKS
HKS.meansd<-
  HKS_TEMP%>%
  summarise(mean.HKS = mean(temp), sd.HKS = sd(temp))

HKS.meansd
#mean.HKS   sd.HKS
#<dbl>      <dbl>
# 28.5      2.07


#OT
OT.meansd<-
  OT_TEMP%>%
  summarise(mean.OT = mean(temp), sd.OT = sd(temp))

OT.meansd
#mean.OT    sd.OT
#<dbl>      <dbl>
# 27.9      2.04

#IT
IT.meansd<-
  IT_TEMP%>%
  summarise(mean.IT = mean(temp), sd.IT = sd(temp))

IT.meansd
#mean.IT    sd.IT
#<dbl>      <dbl>
# 27.2      1.93


#Omit HKS from analysis due to missing data in some seasons

day_sum.noHKS <-
  day_sum%>%
  filter(site != "HKS")

summary(day_sum.noHKS)

####*daily mean mixed model####

meanday.mm <- 
  lmer(TEMP_mean~site*season + (1|MONTHID), data=day_sum.noHKS)

summary(meanday.mm)

#to get summary table of model
anova(meanday.mm)

#check assumptions
plot(fitted(meanday.mm),residuals(meanday.mm))
hist(residuals(meanday.mm ))
qqnorm(residuals(meanday.mm ))
vif(meanday.mm)


#ACF 
acf(residuals(meanday.mm))


#posthoc
emmeans(meanday.mm, pairwise~site|season)


#daily mean linear model
#just to take a look at ACF if month is not used as a random effect

my_lm_mean <- 
  lm(TEMP_mean~site*season, data=day_sum.noHKS)

summary(my_lm_mean)


#check assumptions
plot(fitted(my_lm_mean),residuals(my_lm_mean))
hist(residuals(my_lm_mean))
qqnorm(residuals(my_lm_mean))
vif(my_lm_mean)

#ACF 
acf(residuals(my_lm_mean))



####*daily range mixed model####

dayrange.mm <- 
  lmer(log(day_range)~site*season + (1|MONTHID), data=day_sum.noHKS)

summary(dayrange.mm)


#to get summary table of model
anova(dayrange.mm)

#check assumptions
plot(fitted(dayrange.mm),residuals(dayrange.mm))
hist(residuals(dayrange.mm))
qqnorm(residuals(dayrange.mm))
vif(dayrange.mm)

#ACF 
acf(residuals(dayrange.mm))


#posthoc
emmeans(dayrange.mm, pairwise~site|season)


#daily range linear model
#just to take a look at ACF if month is not used as a random effect

my_lm_range <- 
  lm(log(day_range)~site*season, data=day_sum.noHKS)

summary(my_lm_range)


#check assumptions
plot(fitted(my_lm_range),residuals(my_lm_range))
hist(residuals(my_lm_range))
qqnorm(residuals(my_lm_range))
vif(my_lm_range)

#ACF 
acf(residuals(my_lm_range))


####*daily max mixed model####

daymax.mm <- 
  lmer(TEMP_max~site*season + (1|MONTHID), data=day_sum.noHKS)

summary(daymax.mm)

#check assumptions
plot(fitted(daymax.mm),residuals(daymax.mm))
hist(residuals(daymax.mm))
qqnorm(residuals(daymax.mm))
vif(daymax.mm)

#ACF 
acf(residuals(daymax.mm))


#posthoc
emmeans(daymax.mm, pairwise~site|season)


#daily max linear model
#just to take a look at ACF if month is not used as a random effect

my_lm_max <- 
  lm(TEMP_max~site*season, data=day_sum.noHKS)

summary(my_lm_max)


#check assumptions
plot(fitted(my_lm_max),residuals(my_lm_max))
hist(residuals(my_lm_max))
qqnorm(residuals(my_lm_max))
vif(my_lm_max)

#ACF 
acf(residuals(my_lm_max))


####*daily min mixed model####

daymin.mm <- 
  lmer(TEMP_min~site*season + (1|MONTHID), data=day_sum.noHKS)

summary(daymin.mm)


#check assumptions
plot(fitted(daymin.mm),residuals(daymin.mm))
hist(residuals(daymin.mm))
qqnorm(residuals(daymin.mm))
vif(daymin.mm)

#ACF 
acf(residuals(daymin.mm))


#posthoc
emmeans(daymin.mm, pairwise~site|season)


#daily min linear model
#just to take a look at ACF if month is not used as a random effect

my_lm_min <- 
  lm(TEMP_min~site*season, data=day_sum)

summary(my_lm_min)


#check assumptions
plot(fitted(my_lm_min),residuals(my_lm_min))
hist(residuals(my_lm_min))
qqnorm(residuals(my_lm_min))
vif(my_lm_min)

#ACF 
acf(residuals(my_lm_min))


####____________________####  

####<DHW insitu>####

####1) Libraries####

library(tidyverse)
library(patchwork)
library(zoo)

####2) Data####

#Long-term data from NOAA
#https://coralreefwatch.noaa.gov/product/vs/data/southern_taiwan.txt
#Average maximum monthly mean for southern Taiwan from 1985-2023:
#28.5874

#raw NOAA data
NOAA_temp_southern_taiwan <- read_csv("NOAA_temp_southern_taiwan.csv")


#in-situ data with hotspots and DHW data (all hotspots >1C were divided by 7 
#to calculate cumulative DHW over rolling 84 day intervals)

DHW_insitu <- read_csv("DHW_insitu.csv")

#site subsets (>1C above MMM)

HK.DHW <-
  DHW_insitu%>%
  filter(site =="HK")

HKS.DHW <-
  DHW_insitu%>%
  filter(site =="HKS")

IT.DHW <-
  DHW_insitu%>%
  filter(site =="IT")

OT.DHW <-
  DHW_insitu%>%
  filter(site =="OT")


#in-situ data with hotspots and DHW data (all hotspots > MMM were divided by 7 
#to calculate cumulative DHW over rolling 63 day intervals)

DHW_insitu_NEW <- read_csv("DHW_insitu_NEW.csv")


#site subsets (>MMM)

HK.DHW_NEW <-
  DHW_insitu_NEW%>%
  filter(site =="HK")

HKS.DHW_NEW <-
  DHW_insitu_NEW%>%
  filter(site =="HKS")

IT.DHW_NEW <-
  DHW_insitu_NEW%>%
  filter(site =="IT")

OT.DHW_NEW <-
  DHW_insitu_NEW%>%
  filter(site =="OT")



####3) Plots####

####*HK####
HK.DHW_tempsum <-
  HK.DHW %>%
  mutate(cumsum = cumsum(DHW_data)) %>%
  mutate(cum_rolling84 = rollapplyr(DHW_data, width = 84, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

HK.DHW_tempsum

summary(HK.DHW_tempsum)

HK_DHW_plot <- 
  ggplot(data=HK.DHW_tempsum, aes(x=date, y=cum_rolling84)) +
  geom_line(size=1.5, colour='seagreen')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,21), breaks=c(0, 5, 10, 15, 20))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Haikou reef", y="Degree heating weeks (°C)", x="Month")

HK_DHW_plot


####*HK_NEW####
HK.DHW_tempsum_NEW <-
  HK.DHW_NEW %>%
  mutate(cumsum = cumsum(DHW_data_above_MMM)) %>%
  mutate(cum_rolling63 = rollapplyr(DHW_data_above_MMM, width = 63, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

HK.DHW_tempsum_NEW

summary(HK.DHW_tempsum_NEW)

HK_DHW_plot_NEW <- 
  ggplot(data=HK.DHW_tempsum_NEW, aes(x=date, y=cum_rolling63)) +
  geom_line(size=1.5, colour='seagreen')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,21), breaks=c(0, 5, 10, 15, 20))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Haikou reef", y="Degree heating weeks (°C)", x="Month")

HK_DHW_plot_NEW


####*HKS####

HKS.DHW_tempsum <-
  HKS.DHW %>%
  mutate(cumsum = cumsum(DHW_data)) %>%
  mutate(cum_rolling84 = rollapplyr(DHW_data, width = 84, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

HKS.DHW_tempsum

summary(HKS.DHW_tempsum)

HKS_DHW_plot <- 
  ggplot(data=HKS.DHW_tempsum, aes(x=date, y=cum_rolling84)) +
  geom_line(size=1.5, colour='lightgreen')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,21), breaks=c(0, 5, 10, 15, 20))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Haikou reef flat", y="Degree heating weeks  (°C)", x="Month")

HKS_DHW_plot

####*HKS_NEW####

HKS.DHW_tempsum_NEW <-
  HKS.DHW_NEW %>%
  mutate(cumsum = cumsum(DHW_data_above_MMM)) %>%
  mutate(cum_rolling63 = rollapplyr(DHW_data_above_MMM, width = 63, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

HKS.DHW_tempsum_NEW

summary(HKS.DHW_tempsum_NEW)

HKS_DHW_plot_NEW <- 
  ggplot(data=HKS.DHW_tempsum_NEW, aes(x=date, y=cum_rolling63)) +
  geom_line(size=1.5, colour='lightgreen')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,21), breaks=c(0, 5, 10, 15, 20))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Haikou reef flat", y="Degree heating weeks  (°C)", x="Month")

HKS_DHW_plot_NEW


####*IT####

IT.DHW_tempsum <-
  IT.DHW %>%
  mutate(cumsum = cumsum(DHW_data)) %>%
  mutate(cum_rolling84 = rollapplyr(DHW_data, width = 84, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

IT.DHW_tempsum

summary(IT.DHW_tempsum)

IT_DHW_plot <- 
  ggplot(data=IT.DHW_tempsum, aes(x=date, y=cum_rolling84)) +
  geom_line(size=1.5, colour='steelblue')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,21), breaks=c(0, 5, 10, 15, 20))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Inlet reef", y="Degree heating weeks  (°C)", x="Month")

IT_DHW_plot


####*IT_NEW####

IT.DHW_tempsum_NEW <-
  IT.DHW_NEW %>%
  mutate(cumsum = cumsum(DHW_data_above_MMM)) %>%
  mutate(cum_rolling63 = rollapplyr(DHW_data_above_MMM, width = 63, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

IT.DHW_tempsum_NEW

summary(IT.DHW_tempsum_NEW)

IT_DHW_plot_NEW <- 
  ggplot(data=IT.DHW_tempsum_NEW, aes(x=date, y=cum_rolling63)) +
  geom_line(size=1.5, colour='steelblue')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,21), breaks=c(0, 5, 10, 15, 20))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Inlet reef", y="Degree heating weeks  (°C)", x="Month")

IT_DHW_plot_NEW



####*OT####

OT.DHW_tempsum <-
  OT.DHW %>%
  mutate(cumsum = cumsum(DHW_data)) %>%
  mutate(cum_rolling84 = rollapplyr(DHW_data, width = 84, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

OT.DHW_tempsum

summary(OT.DHW_tempsum)

OT_DHW_plot <- 
  ggplot(data=OT.DHW_tempsum, aes(x=date, y=cum_rolling84)) +
  geom_line(size=1.5, colour='red')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,25), breaks=c(0, 5, 10, 15, 20, 25))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Outlet reef", y="Degree heating weeks  (°C)", x="Month")

OT_DHW_plot


####*OT_NEW####

OT.DHW_tempsum_NEW <-
  OT.DHW_NEW %>%
  mutate(cumsum = cumsum(DHW_data_above_MMM)) %>%
  mutate(cum_rolling63 = rollapplyr(DHW_data_above_MMM, width = 63, FUN = sum, partial = TRUE)) %>%
  drop_na(cumsum)

OT.DHW_tempsum_NEW

summary(OT.DHW_tempsum_NEW)

OT_DHW_plot_NEW <- 
  ggplot(data=OT.DHW_tempsum_NEW, aes(x=date, y=cum_rolling63)) +
  geom_line(size=1.5, colour='red')+
  theme_classic()+
  theme(text=element_text(size=12,  family="sans"))+
  geom_hline(yintercept=15, linetype="dashed", 
             color = "black", size=0.5)+
  scale_y_continuous(limits=c(0,25), breaks=c(0, 5, 10, 15, 20, 25))+
  scale_x_date(date_labels="%b",date_breaks  ="1 month", limits = as.Date(c("2022-03-16", "2023-02-21")))+
  labs(title="Outlet reef", y="Degree heating weeks  (°C)", x="Month")

OT_DHW_plot_NEW


####*combo####

#>1C above MMM
DHW_combo <- HK_DHW_plot / HKS_DHW_plot / IT_DHW_plot / OT_DHW_plot
DHW_combo


#TEMP > MMM
DHW_combo_NEW <- HK_DHW_plot_NEW / HKS_DHW_plot_NEW / IT_DHW_plot_NEW / OT_DHW_plot_NEW
DHW_combo_NEW

#DHW comparison combo plot
DHW_combo | DHW_combo_NEW


#combine > MMM with raw data
raw_DHW_COMBO_PLOT <- raw_temp_combo_plot | DHW_combo_NEW
raw_DHW_COMBO_PLOT


####____________________####   

####< INSITU LIGHT>####

####1) Libraries####

library(tidyverse)
library(patchwork)

####2) Data####

insitu_light_data <- read_csv("insitu_light_data.csv")

summary(insitu_light_data)

insitu_light_data$day_night <- as.factor(insitu_light_data$day_night)
insitu_light_data$DateTime <-  sub(" UTC", "", insitu_light_data$DateTime)
insitu_light_data$DateTime <- as.POSIXct(insitu_light_data$DateTime,format="%Y-%m-%d%H:%M:%S")
insitu_light_data$site <- as.factor(insitu_light_data$site)

summary(insitu_light_data)
head(insitu_light_data)

####3) Plots####

#summary table for mean (sd) light during the day vs. night
insitu_light_data.sum <-
  insitu_light_data%>%
  group_by(day_night, site)%>%
  summarise(mean_PAR = mean(PAR), sd_PAR = sd(PAR))

insitu_light_data.sum


#subsets
insitu_light_data.OT <-
  insitu_light_data%>%
  filter(site =="Outlet")

insitu_light_data.IT <-
  insitu_light_data%>%
  filter(site =="Inlet")

insitu_light_data.HK <-
  insitu_light_data%>%
  filter(site =="Haikou")

#OT
insitu_light_plot.OT <- 
  ggplot(insitu_light_data.OT, aes(x=DateTime, y=PAR)) + 
  geom_line(size=0.05, colour = "red")+
  theme_classic()+
  theme(text=element_text(size=14,  family="sans"))+
  geom_hline(yintercept=200, linetype="dashed", 
             color = "black", size=0.75)+
  scale_y_continuous(limits=c(0, 1200), breaks=c(0, 200, 400, 600, 800, 1000, 1200))+
  labs(title="A. Outlet reef", y="PAR", x="Date")

insitu_light_plot.OT


#IT
insitu_light_plot.IT <- 
  ggplot(insitu_light_data.IT, aes(x=DateTime, y=PAR)) + 
  geom_line(size=0.05, colour = "steelblue")+
  theme_classic()+
  theme(text=element_text(size=14,  family="sans"))+
  geom_hline(yintercept=200, linetype="dashed", 
             color = "black", size=0.75)+
  scale_y_continuous(limits=c(0, 1200), breaks=c(0, 200, 400, 600, 800, 1000, 1200))+
  labs(title="B. Inlet reef", y="PAR", x="Date")

insitu_light_plot.IT


#HK
insitu_light_plot.HK <- 
  ggplot(insitu_light_data.HK, aes(x=DateTime, y=PAR)) + 
  geom_line(size=0.05, colour="seagreen")+
  theme_classic()+
  theme(text=element_text(size=14,  family="sans"))+
  geom_hline(yintercept=200, linetype="dashed", 
             color = "black", size=0.75)+
  scale_y_continuous(limits=c(0, 1200), breaks=c(0, 200, 400, 600, 800, 1000, 1200))+
  labs(title="C. Haikou reef", y="PAR", x="Date")

insitu_light_plot.HK

light_plot <- insitu_light_plot.OT | insitu_light_plot.IT | insitu_light_plot.HK

light_plot





####____________________#### 

####<INSITU NUTRIENTS>####


####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)


####2) Data####

one_year_water_chem <- read_csv("one_year_water_chem.csv")
insitu_water_chem <- one_year_water_chem 

#class check

#WIDE format
summary(insitu_water_chem)

insitu_water_chem$site <- as.factor(insitu_water_chem$site)
insitu_water_chem$bottle_ID <- as.factor(insitu_water_chem$bottle_ID)
insitu_water_chem$collection_date <- as.Date(insitu_water_chem$collection_date)

str(insitu_water_chem)
summary(insitu_water_chem)
head(insitu_water_chem)


#convert to long format
insitu_water_chem_LONG <-
  pivot_longer(one_year_water_chem, -c(site, bottle_ID, collection_date), 
               values_to = "value", names_to = "parameter")



#LONG_format
summary(insitu_water_chem_LONG)

insitu_water_chem_LONG$site <- as.factor(insitu_water_chem_LONG$site)
insitu_water_chem_LONG$bottle_ID <- as.factor(insitu_water_chem_LONG$bottle_ID)
insitu_water_chem_LONG$collection_date <- as.Date(insitu_water_chem_LONG$collection_date)
insitu_water_chem_LONG$parameter <- as.factor(insitu_water_chem_LONG$parameter)

str(insitu_water_chem_LONG)
summary(insitu_water_chem_LONG)
head(insitu_water_chem_LONG)


####3) Plots####

#remove NAs
insitu_water_chem_LONG <- insitu_water_chem_LONG%>% filter(!is.na(value))

#site summary
water_chem_sum <-
  insitu_water_chem_LONG%>%
  group_by(site, parameter, collection_date)%>%
  summarise(mean = mean(value), n = n(), sd = sd(value), se = sd / sqrt(n))

water_chem_sum 


#site summary (for table)
insitu_water_chem_site_sum <-
  insitu_water_chem_LONG%>%
  group_by(site, parameter)%>%
  summarise(mean = mean(value), n = n(), sd = sd(value), se = sd / sqrt(n))

insitu_water_chem_site_sum


#water chemistry boxplot
water_chem_plot <- 
  ggplot(data=insitu_water_chem_LONG, aes(x=site, y=value, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, height=0) +
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red"))+
  labs(y = "", x="")+
  theme_bw()+
  facet_wrap(~parameter, scales="free", ncol = 2)

water_chem_plot 

water_chem_plot +  theme(text=element_text(size=14,  family="sans")) + theme(legend.position = "none")


####4) Analysis####

####*Alkalinity####
ALK_mm <- 
  lmer(ALK~site + (1|collection_date), data=insitu_water_chem)

summary(ALK_mm)

#check assumptions
plot(fitted(ALK_mm),residuals(ALK_mm))
hist(residuals(ALK_mm))
qqnorm(residuals(ALK_mm))

#ACF 
acf(residuals(ALK_mm))


#posthoc
emmeans(ALK_mm, pairwise~site)


####*Ammonia####
AMMO_mm <- 
  lmer(AMMO~site + (1|collection_date), data=insitu_water_chem)

summary(AMMO_mm)


#check assumptions
plot(fitted(AMMO_mm),residuals(AMMO_mm))
hist(residuals(AMMO_mm))
qqnorm(residuals(AMMO_mm))

#ACF 
acf(residuals(AMMO_mm))


#posthoc
emmeans(AMMO_mm, pairwise~site)



####*Calcium####
Ca_mm <- 
  lmer(Ca~site + (1|collection_date), data=insitu_water_chem)

summary(Ca_mm)


#check assumptions
plot(fitted(Ca_mm),residuals(Ca_mm))
hist(residuals(Ca_mm))
qqnorm(residuals(Ca_mm))

#ACF 
acf(residuals(Ca_mm))


#posthoc
emmeans(Ca_mm, pairwise~site)




####*Magnesium####
Mg_mm <- 
  lmer(Mg~site + (1|collection_date), data=insitu_water_chem)

summary(Mg_mm)


#check assumptions
plot(fitted(Mg_mm),residuals(Mg_mm))
hist(residuals(Mg_mm))
qqnorm(residuals(Mg_mm))

#ACF 
acf(residuals(Mg_mm))


#posthoc
emmeans(Mg_mm, pairwise~site)



####*Nitrate####


NITRATE_mm <- 
  lmer(NITRATE~site + (1|collection_date), data=insitu_water_chem)

summary(NITRATE_mm)


#check assumptions
plot(fitted(NITRATE_mm),residuals(NITRATE_mm))
hist(residuals(NITRATE_mm))
qqnorm(residuals(NITRATE_mm))

#ACF 
acf(residuals(NITRATE_mm))


#posthoc
emmeans(NITRATE_mm, pairwise~site)



####*Nitrite####

# all values are zero; no analysis needed


####*pH####
pH_mm <- 
  lmer(pH~site + (1|collection_date), data=insitu_water_chem)

summary(pH_mm)


#check assumptions
plot(fitted(pH_mm),residuals(pH_mm))
hist(residuals(pH_mm))
qqnorm(residuals(pH_mm))

#ACF 
acf(residuals(pH_mm))


#posthoc
emmeans(pH_mm, pairwise~site)



####*Phosphate####
PHOS_mm <- 
  lmer(PHOS~site + (1|collection_date), data=insitu_water_chem)

summary(PHOS_mm)


#check assumptions
plot(fitted(PHOS_mm),residuals(PHOS_mm))
hist(residuals(PHOS_mm))
qqnorm(residuals(PHOS_mm))

#ACF 
acf(residuals(PHOS_mm))


#posthoc
emmeans(PHOS_mm, pairwise~site)



####*********************************####
####HEAT-STRES EXPERIMENT####
####*********************************####

####<RECOVERY TANK>####

#colonies (and subsequent nubbins) were held in a recovery tank from March 16 2022 - March 31 2022 

####1) Libraries####

library(tidyverse)
library(patchwork)

####2) Data####

recovery_tank_temp_light <- read_csv("recovery_tank_temp_light.csv")

summary(recovery_tank_temp_light)

recovery_tank_temp_light$DateTime <-as.POSIXct(recovery_tank_temp_light$DateTime,format="%Y-%m-%d%H:%M:%S",tz=Sys.timezone())
recovery_tank_temp_light$day_night <-as.factor(recovery_tank_temp_light$day_night)

summary(recovery_tank_temp_light)
head(recovery_tank_temp_light)


####3) Plots####

####*Temperature####

#mean and sd temperature
mean(recovery_tank_temp_light$temp)
#26.18512

sd(recovery_tank_temp_light$temp)
#0.6921829

#recovery tank temp plot
recovery_tank_temp <-
  ggplot(recovery_tank_temp_light, aes(x=DateTime, y=temp)) +
  geom_line(size=.75)+
  scale_y_continuous(limits=c(24, 28), breaks=c(24, 25, 26, 27, 28))+
  scale_x_datetime(date_minor_breaks = "1 day", date_breaks = "2 day",date_labels = "%b%d") +
  theme_classic()+
  labs(title="A. Recovery Tank Temperature", y="Temperature (?C)", x="Date")+
  theme(text=element_text(size=14,  family="sans"))

recovery_tank_temp



####*Light####

#subset for day
recovery_tank_temp_light_DAY <-
  recovery_tank_temp_light%>%
  filter(day_night =="day")

mean(recovery_tank_temp_light_DAY$PAR)
#67.95676

sd(recovery_tank_temp_light_DAY$PAR)
#89.91996

#recovery tank light plot
recovery_light_plot <- 
  ggplot(recovery_tank_temp_light, aes(x=DateTime, y=PAR)) + 
  geom_line(size=0.05, colour = "black")+
  theme_classic()+
  scale_x_datetime(date_minor_breaks = "1 day", date_breaks = "2 day",date_labels = "%b%d") +
  theme(text=element_text(size=14,  family="sans"))+
  geom_hline(yintercept=200, linetype="dashed", 
             color = "black", size=0.75)+
  scale_y_continuous(limits=c(0, 800), breaks=c(0, 200, 400, 600, 800))+
  labs(title="B. Recovery Tank Light", y="PAR", x="Date")

recovery_light_plot

recovey_cond_plot <- recovery_tank_temp / recovery_light_plot

recovey_cond_plot


####____________________#### 

####<EXP TANK TEMP>####

####1) Libraries####

library(tidyverse)

####2) Data####

tank_temp <- read_csv("tank_temp.csv")

#change wide format to long format
tank_temp_long <- gather(tank_temp, tank_treatment, temp, T1_30:T15_27, factor_key=TRUE)
tank_temp_long 

#add a seperate treatment column to add plot colour
tank_temp_long <-
  tank_temp_long%>%
  mutate(Temperature = case_when(tank_treatment == 'T13_27'| tank_treatment =='T14_27'| tank_treatment == 'T15_27' ~ 27,
                                 tank_treatment == 'T1_30'| tank_treatment =='T5_30'| tank_treatment == 'T9_30' ~ 30,  
                                 tank_treatment == 'T2_33'| tank_treatment == 'T6_33'| tank_treatment == 'T10_33' ~ 33, 
                                 tank_treatment == 'T3_36'| tank_treatment =='T7_36'| tank_treatment == 'T11_36' ~ 36, 
                                 tank_treatment == 'T4_39'| tank_treatment =='T8_39'| tank_treatment == 'T12_39' ~ 39))


summary(tank_temp_long)

tank_temp_long$DateTime <- as.POSIXct(tank_temp_long$DateTime,tz=Sys.timezone())
tank_temp_long$DateTime <- as.POSIXct(tank_temp_long$DateTime, "%Y-%m-%d  %H:%M:%S",tz=Sys.timezone())
tank_temp_long$tank_treatment <- as.factor(tank_temp_long$tank_treatment)
tank_temp_long$Temperature <- as.factor(tank_temp_long$Temperature)


summary(tank_temp_long)
head(tank_temp_long)


####3) Plot####

tank_temp_plot <-
  ggplot(tank_temp_long, aes(x=DateTime, y=temp, group=tank_treatment, color=Temperature)) +
  geom_line(size=1)+
  scale_y_continuous(limits=c(26, 40), breaks=c(27, 30, 33, 36, 39))+
  scale_x_datetime(date_minor_breaks = "2 hour", date_breaks = "1 hour",date_labels = "%H") +
  scale_color_manual(values=c("seagreen", "blue", "yellow","orange", "red"))+
  theme_classic()+
  labs(title="", y="Temperature (?C)", x="Hour")+
  theme(text=element_text(size=13,  family="sans"))

tank_temp_plot 

tank_temp_plot + theme(legend.position = "bottom")

tank_temp_plot + facet_wrap(~tank_treatment) + theme_bw() + theme(legend.position = "top")



####____________________#### 

####<COLONY DATA>####

####1) Libraries####

library(tidyverse)

####2) Data####

colony_post_collection <- read_csv("colony_post_collection.csv")

#class check
summary(colony_post_collection)

colony_post_collection$colony_ID <- as.factor(colony_post_collection$colony_ID)
colony_post_collection$site <- as.factor(colony_post_collection$site)
colony_post_collection$collection_date <- as.Date(colony_post_collection$collection_date)

summary(colony_post_collection)
head(colony_post_collection)

####3) Table####

colony_sum <-
  colony_post_collection%>%
  group_by(site)%>%
  summarise(mean.size = mean(diameter_cm), sd.size = sd(diameter_cm), mean.fvfm =mean(FvFm_avg), sd.fvfm = sd (FvFm_avg))

colony_sum

# A tibble: 4 x 5
#  site             mean.size sd.size mean.fvfm sd.fvfm
#  <fct>                <dbl>   <dbl>     <dbl>   <dbl>
#1 Haikou                11.5   1.97      0.644  0.0272
#2 Haikou_reef.flat      12.5   2.43      0.656  0.0265
#3 Inlet                 14     0.632     0.626  0.0164
#4 Outlet                13.8   1.47      0.649  0.0319

mean(colony_post_collection$diameter_cm)
#12.95833
sd(colony_post_collection$diameter_cm)
#1.944427

mean(colony_post_collection$FvFm_avg)
#0.643875
sd(colony_post_collection$FvFm_avg)
#0.0270158


####4) Analysis####

size.anova <- aov(diameter_cm ~ site, data=colony_post_collection)
summary(size.anova)
plot(size.anova)
TukeyHSD(size.anova)

fvfm.anova <- aov(FvFm_avg ~ site, data=colony_post_collection)
summary(fvfm.anova)
plot(fvfm.anova)
TukeyHSD(fvfm.anova)


####____________________#### 

####<CHL A>####


####1) Libraries####

library(tidyverse)
library(emmeans)
library(car)

####2) Data####

chl_a_data <- read_csv("chl_a_data.csv")

summary(chl_a_data)

chl_a_data$coral_id <-as.factor(chl_a_data$coral_id)
chl_a_data$treatment <-as.factor(chl_a_data$treatment)
chl_a_data$site <-as.factor(chl_a_data$site)

summary(chl_a_data)
head(chl_a_data)

####3) Plots####

chl_a_plot <- 
  ggplot(chl_a_data, aes(x=treatment, y=chla_ug_cm2, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.1), alpha=0.4, size =1)+
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red"))+
  scale_y_continuous(limits=c(0, 5), breaks=c(0, 1, 2, 3, 4, 5))+
  labs(y = "Chlorophyll a (ug/cm^2)", x="Temperature (?C)")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_classic()



CHL <- chl_a_plot + theme(legend.position = "none") +  theme(text=element_text(size=16,  family="sans"))
CHL
chl_a_plot +  facet_wrap(~site) + theme_bw()



####4) Analysis####

plot(chl_a_data$chla_ug_cm2)
hist(chl_a_data$chla_ug_cm2)

#model 1 = no transformation
chl_a.lm <- lm(chla_ug_cm2 ~ site * treatment, data = chl_a_data)
Anova(chl_a.lm, type = "III")

#check assumptions
plot(fitted(chl_a.lm),residuals(chl_a.lm))
hist(residuals(chl_a.lm))
qqnorm(residuals(chl_a.lm))
vif(chl_a.lm)

#evidence of model assumption violation; use sqrt transformation:

#model 2 = sqrt transformation
chl_a.lm.sqrt <- lm(sqrt(chla_ug_cm2) ~ site * treatment, data = chl_a_data)
Anova(chl_a.lm.sqrt, type = "III")

#check assumptions
plot(fitted(chl_a.lm.sqrt),residuals(chl_a.lm.sqrt))
hist(residuals(chl_a.lm.sqrt))
qqnorm(residuals(chl_a.lm.sqrt))
vif(chl_a.lm.sqrt)

#assumptions are better met with sqrt transformation

#posthoc 
chl_a.lm.sqrt.posthoc1 <- emmeans(chl_a.lm.sqrt, pairwise ~ site|treatment, type="response")
summary(chl_a.lm.sqrt.posthoc1)

chl_a.lm.sqrt.posthoc2 <- emmeans(chl_a.lm.sqrt, pairwise ~ treatment|site, type="response")
summary(chl_a.lm.sqrt.posthoc2)





####____________________#### 

####<SYM DENSITY>####


####1) Libraries####

library(tidyverse)
library(emmeans)
library(car)

####2) Data####

sym_density_data <- read_csv("sym_density_data.csv")

summary(sym_density_data)

sym_density_data$coral_id <-as.factor(sym_density_data$coral_id)
sym_density_data$treatment <-as.factor(sym_density_data$treatment)
sym_density_data$site <-as.factor(sym_density_data$site)

summary(sym_density_data)
head(sym_density_data)

####3) Plots####

density_plot <- 
  ggplot(sym_density_data, aes(x=treatment, y=sym_density, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.1), alpha=0.4, size =1)+
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red"))+
  scale_y_continuous(limits=c(0, 22), breaks=c(0, 5, 10, 15, 20))+
  labs(y = "Symbiont density (10^5 cells/cm^2)", x="Temperature (?C)")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_classic()

density_plot 
DEN <- density_plot + theme(legend.position = "none") +  theme(text=element_text(size=16,  family="sans"))
DEN 
density_plot +  facet_wrap(~site) + theme_bw()


####4) Analysis####

plot(sym_density_data$sym_density)
hist(sym_density_data$sym_density)

#model 1 = no transformation
sym_density.lm <- lm(sym_density ~ site * treatment, data = sym_density_data)
Anova(sym_density.lm, type = "III")

#check assumptions
plot(fitted(sym_density.lm),residuals(sym_density.lm))
hist(residuals(sym_density.lm))
qqnorm(residuals(sym_density.lm))
vif(sym_density.lm)

#evidence of model assumption violation; use a sqrt transformation


#model 2 = sqrt transformation
sym_density.lm.sqrt <- lm(sqrt(sym_density) ~ site * treatment, data = sym_density_data)
Anova(sym_density.lm.sqrt, type = "III")

#check assumptions
plot(fitted(sym_density.lm.sqrt),residuals(sym_density.lm.sqrt))
hist(residuals(sym_density.lm.sqrt))
qqnorm(residuals(sym_density.lm.sqrt))
vif(sym_density.lm.sqrt)

#assumptions are better met with sqrt transformation


#posthoc 
sym_density.posthoc1 <- emmeans(sym_density.lm.sqrt, pairwise ~ site|treatment, type="response")
summary(sym_density.posthoc1)

sym_density.posthoc2 <- emmeans(sym_density.lm.sqrt, pairwise ~ treatment|site, type="response")
summary(sym_density.posthoc2)



####____________________#### 

####<PROTEIN>####


####1) Libraries####

library(tidyverse)
library(emmeans)
library(car)

####2) Data####

protein_data <- read_csv("protein_data.csv")

summary(protein_data)

protein_data$coral_id <-as.factor(protein_data$coral_id)
protein_data$treatment <-as.factor(protein_data$treatment)
protein_data$site <-as.factor(protein_data$site)

summary(protein_data)
head(protein_data)

####3) Plots####

protein_plot <- 
  ggplot(protein_data, aes(x=treatment, y=protein_mg_cm2, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.1), alpha=0.4, size =1)+
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red"))+
  scale_y_continuous(limits=c(0, 0.35), breaks=c(0, .1, .2, .3))+
  labs(y = "Host protein (mg/cm^2)", x="Temperature (?C)")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_classic()

protein_plot 
PRO <- protein_plot  + theme(legend.position = "none") +  theme(text=element_text(size=16,  family="sans"))
PRO
protein_plot  +  facet_wrap(~site) + theme_bw()


####4) Analysis####

plot(protein_data$protein_mg_cm2)
hist(protein_data$protein_mg_cm2)

#model 1 = no transformation
protein.lm <- lm(protein_mg_cm2 ~ site * treatment, data = protein_data)
Anova(protein.lm, type = "III")

#check assumptions
plot(fitted(protein.lm),residuals(protein.lm))
hist(residuals(protein.lm))
qqnorm(residuals(protein.lm))
vif(protein.lm)

#evidence of model assumption violation; use a sqrt transformation


#model 2 = sqrt transformation
protein.lm.sqrt <- lm(sqrt(protein_mg_cm2) ~ site * treatment, data = protein_data)
Anova(protein.lm.sqrt, type = "III")

#check assumptions
plot(fitted(protein.lm.sqrt),residuals(protein.lm.sqrt))
hist(residuals(protein.lm.sqrt))
qqnorm(residuals(protein.lm.sqrt))
vif(protein.lm.sqrt)

#assumptions are better met with sqrt transformation


#posthoc 
protein.posthoc1 <- emmeans(protein.lm.sqrt, pairwise ~ site|treatment, type="response")
summary(protein.posthoc1)

protein.posthoc2 <- emmeans(protein.lm.sqrt, pairwise ~ treatment|site, type="response")
summary(protein.posthoc2)



####____________________#### 


####<CBASS BLEACH>####


####1) Libraries####

library(tidyverse)
library(emmeans)
library(car)
library(chisq.posthoc.test)
library(corrplot)


####2) Data####

CBASS_bleach <- read_csv("CBASS_bleach.csv")

summary(CBASS_bleach)

CBASS_bleach$tank <- as.factor(CBASS_bleach$tank)
CBASS_bleach$temp <- as.factor(CBASS_bleach$temp)
CBASS_bleach$treatment <- as.factor(CBASS_bleach$treatment)
CBASS_bleach$site <- as.factor(CBASS_bleach$site)
CBASS_bleach$colony_id <- as.factor(CBASS_bleach$colony_id)
CBASS_bleach$nubbin_id <- as.factor(CBASS_bleach$nubbin_id)
CBASS_bleach$score_1 <- as.factor(CBASS_bleach$score_1)
CBASS_bleach$score_2 <- as.factor(CBASS_bleach$score_2)
CBASS_bleach$score_3 <- as.factor(CBASS_bleach$score_3)
CBASS_bleach$median_score <- as.factor(CBASS_bleach$median_score)

summary(CBASS_bleach)
head(CBASS_bleach)


bleach_percent <-
  CBASS_bleach%>%
  group_by(site, treatment, median_score)%>%
  summarise(n=n())

#if error occurs while running the above code it means that previously loaded package is masking the n=n() 
#function; close and open R again to run properly

bleach_percent <-
  bleach_percent%>%
  group_by(site, treatment)%>%
  mutate(sum = sum(n))  

bleach_percent <-
  bleach_percent%>%
  group_by(site, treatment)%>%
  mutate(percent = n/sum*100) 

HK_subset.bleach <-
  bleach_percent%>%
  filter(site=='HK')

OT_subset.bleach <-
  bleach_percent%>%
  filter(site=='OT')

IT_subset.bleach <-
  bleach_percent%>%
  filter(site=='IT')


####3) Plots####

bleach_percent$median_score <- factor(bleach_percent$median_score, levels=c("4", "3", "2", "1"))


cBASS_bleach_plot <-
  ggplot(bleach_percent, aes(fill=median_score, y=percent, x=treatment)) + 
  geom_bar(position="stack", stat="identity", colour="black")+
  scale_fill_manual(values=c("white", "bisque1","bisque3", "bisque4"))+
  theme_bw()+
  labs(y = "Percentage of colonies", x="Temperature ?C")+
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100))

cBASS_bleach_plot

cBASS_bleach_plot + facet_wrap(~site)  + theme(legend.position = "top")

cBASS_bleach_plot + facet_wrap(~site)  + theme(legend.position = "none") +  theme(text=element_text(size=16,  family="sans"))



####4) Analysis####

CBASS_bleach_check <- table(CBASS_bleach$site, CBASS_bleach$median_score)

CBASS_bleach_check

mosaicplot(CBASS_bleach_check)

prop.table(CBASS_bleach_check, margin =1)

addmargins(CBASS_bleach_check)

CBASS_chisq <- chisq.test(CBASS_bleach_check)

CBASS_chisq

corrplot(CBASS_chisq$residuals, is.cor = FALSE)

chisq.posthoc.test(CBASS_bleach_check)



####____________________#### 


####<FVFM>####

####1) Libraries####

library(tidyverse)
library(drc)
library(lmerTest)
library(emmeans)
library(Rmisc)
library(reshape2)
library(car)
library(patchwork)


####2) Data####

fvfm_data <- read_csv("fvfm_data.csv")
fvfm <-fvfm_data 

summary(fvfm)

fvfm$tank <-as.factor (fvfm$tank)
fvfm$treatment <-as.factor (fvfm$treatment)
fvfm$temp <-as.numeric(fvfm$temp)
fvfm$site <-as.factor (fvfm$site)
fvfm$colony_id <-as.factor (fvfm$colony_id)
fvfm$nubbin_id <-as.factor (fvfm$nubbin_id)

summary(fvfm)
head(fvfm)


####3) Plots####

#quick look at data (pre-experiment; March 30th)
fvfm_boxplot_pre <-
  ggplot(fvfm, aes(x=treatment, y=pre_fvfm, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.1), alpha=0.4, size =1)+
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red"))+
  scale_y_continuous(limits=c(0, 0.8), breaks=c(0,.2,.4,.6,.8))+
  labs(title="A. CBASS Fv/Fm: pre-experiment")+
  theme_classic()

PRE <- fvfm_boxplot_pre +  theme(text=element_text(size=16,  family="sans"))

fvfm_boxplot_pre + facet_wrap(~site)

#quick look at data (post-experiment; March 31st)
fvfm_boxplot_post <-
  ggplot(fvfm, aes(x=treatment, y=post_fvfm, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.1), alpha=0.4, size =1)+
  scale_fill_manual(values=c("seagreen", "lightgreen", "steelblue","red"))+
  labs(title="B. CBASS Fv/Fm: post-experiment")+
  theme_classic()

POST <- fvfm_boxplot_post +  theme(text=element_text(size=16,  family="sans"))

fvfm_boxplot_post + facet_wrap(~site)

FVFM_COMBO <- PRE / POST 
FVFM_COMBO 

#Use dose curve approach; code used is based on approach described by Evensen et al. 2022
#Evensen NR, Voolstra CR, Fine M, Perna G, Buitrago-Lopez C, Cardenas A, Banc-Prandi G, Rowe K, Barshis DJ.
#(2022) Empirically derived thermal thresholds of four coral species along the Red Sea using a portable and 
#standardized experimental approach. Coral Reefs 41(2): 239-52. [https://doi.org/10.1007/s00338-022-02233-y]


#site subsets

#Haikou
HK_fvfm <-
  fvfm %>%
  filter(site == "HK")

HK_fvfm


#Haikou .shallow
HKS_fvfm <-
  fvfm %>%
  filter(site == "HKS")

HKS_fvfm


#Outlet
OT_fvfm <-
  fvfm %>%
  filter(site == "OT")

OT_fvfm


#Inlet
IT_fvfm <-
  fvfm %>%
  filter(site == "IT")

IT_fvfm


####*HK plots####

#Run population-level fit
Haikou_pop <- drm(post_fvfm ~ temp, data = HK_fvfm, fct = LL.3())
summary(Haikou_pop)
plot(Haikou_pop)

#Run individual colony fits
Haikou_indiv <- drm(post_fvfm ~ temp, data = HK_fvfm, curveid=colony_id, fct = LL.3())
summary(Haikou_indiv)
plot(Haikou_indiv)

#extract coeffs by colony, then compute 95% CIs
Haikou_genocoeffs_50<-data.frame(ED(Haikou_indiv, c(50)))
Haikou_coeff_mean<-mean(Haikou_genocoeffs_50$Estimate)
Haikou_coeff_mean

Haikou_summary<-data.frame(CI(Haikou_genocoeffs_50$Estimate, ci=0.95))
Haikou_coeff_lower<-Haikou_summary[3,]
Haikou_coeff_upper<-Haikou_summary[1,]



####*HKS plots####

#Run population-level fit
Haikou.shallow_pop <- drm(post_fvfm ~ temp, data = HKS_fvfm, fct = LL.3())
summary(Haikou.shallow_pop)
plot(Haikou.shallow_pop)

#Run individual colony fits
Haikou.shallow_indiv <- drm(post_fvfm ~ temp, data = HKS_fvfm, curveid=colony_id, fct = LL.3())
summary(Haikou.shallow_indiv)
plot(Haikou.shallow_indiv)

#extract coeffs by colony, then compute 95% CIs
Haikou.shallow_genocoeffs_50<-data.frame(ED(Haikou.shallow_indiv, c(50)))
Haikou.shallow_coeff_mean<-mean(Haikou.shallow_genocoeffs_50$Estimate)
Haikou.shallow_coeff_mean

Haikou.shallow_summary<-data.frame(CI(Haikou.shallow_genocoeffs_50$Estimate, ci=0.95))
Haikou.shallow_coeff_lower<-Haikou.shallow_summary[3,]
Haikou.shallow_coeff_upper<-Haikou.shallow_summary[1,]


####*OT plots####

#Run population-level fit
Outlet_pop <- drm(post_fvfm ~ temp, data = OT_fvfm, fct = LL.3())
summary(Outlet_pop)
plot(Outlet_pop)

#Run individual colony fits
Outlet_indiv <- drm(post_fvfm ~ temp, data = OT_fvfm, curveid=colony_id, fct = LL.3())
summary(Outlet_indiv)
plot(Outlet_indiv)

#extract coeffs by colony, then compute 95% CIs
Outlet_genocoeffs_50<-data.frame(ED(Outlet_indiv, c(50)))
Outlet_coeff_mean<-mean(Outlet_genocoeffs_50$Estimate)
Outlet_coeff_mean

Outlet_summary<-data.frame(CI(Outlet_genocoeffs_50$Estimate, ci=0.95))
Outlet_coeff_lower<-Outlet_summary[3,]
Outlet_coeff_upper<-Outlet_summary[1,]




####*IT plots####

#Run population-level fit
Inlet_pop <- drm(post_fvfm ~ temp, data = IT_fvfm, fct = LL.3())
summary(Inlet_pop)
plot(Inlet_pop)

#Run individual colony fits
Inlet_indiv <- drm(post_fvfm ~ temp, data = IT_fvfm, curveid=colony_id, fct = LL.3())
summary(Inlet_indiv)
plot(Inlet_indiv)

#extract coeffs by colony, then compute 95% CIs
Inlet_genocoeffs_50<-data.frame(ED(Inlet_indiv, c(50)))
Inlet_coeff_mean<-mean(Inlet_genocoeffs_50$Estimate)
Inlet_coeff_mean

Inlet_summary<-data.frame(CI(Inlet_genocoeffs_50$Estimate, ci=0.95))
Inlet_coeff_lower<-Inlet_summary[3,]
Inlet_coeff_upper<-Inlet_summary[1,]


####*Calculate ED50s####
#combine colony-ED50s into data frame for statistical analysis

ED50s <- data.frame(cbind(Haikou_genocoeffs_50[,1],Haikou.shallow_genocoeffs_50[,1],Outlet_genocoeffs_50[,1],
                          Inlet_genocoeffs_50[,1]))

ED50s <- ED50s %>% 
  dplyr::rename(Haikou= X1,
                Haikou.shallow=X2,
                Outlet=X3,
                Inlet=X4)


ED50s$colony_id<-as.factor(1:nrow(ED50s))
str(ED50s)

ED50s_long<-melt(ED50s, id="colony_id")

ED50s_long<-ED50s_long %>% 
  dplyr::rename(site= variable,
                ED50=value)


####*ED50 plot####

PA_coeff_means<-data.frame(Haikou_coeff_mean, Haikou.shallow_coeff_mean, Outlet_coeff_mean, Inlet_coeff_mean)
PA_coeff_lowers<-data.frame(Haikou_coeff_lower,  Haikou.shallow_coeff_lower, Outlet_coeff_lower, Inlet_coeff_lower)
PA_coeff_uppers<-data.frame(Haikou_coeff_upper,  Haikou.shallow_coeff_upper, Outlet_coeff_upper, Inlet_coeff_upper)

Haikou_preddata = data.frame(temp = seq(27,39, length.out = 400))
Haikou_pred = as.data.frame(predict(Haikou_pop, newdata = Haikou_preddata, interval = 'confidence'))
Haikou_preddata = data.frame(Haikou_preddata, post_fvfm = Haikou_pred$Prediction, Lower = Haikou_pred$Lower, Upper = Haikou_pred$Upper)

Haikou.shallow_preddata = data.frame(temp = seq(27,39, length.out = 400))
Haikou.shallow_pred = as.data.frame(predict(Haikou.shallow_pop, newdata = Haikou.shallow_preddata, interval = 'confidence'))
Haikou.shallow_preddata = data.frame(Haikou.shallow_preddata, post_fvfm = Haikou.shallow_pred$Prediction, Lower = Haikou.shallow_pred$Lower, Upper = Haikou.shallow_pred$Upper)

Outlet_preddata = data.frame(temp = seq(27,39, length.out = 400))
Outlet_pred = as.data.frame(predict(Outlet_pop, newdata = Outlet_preddata, interval = 'confidence'))
Outlet_preddata = data.frame(Outlet_preddata, post_fvfm = Outlet_pred$Prediction, Lower = Outlet_pred$Lower, Upper = Outlet_pred$Upper)

Inlet_preddata = data.frame(temp = seq(27,39, length.out = 400))
Inlet_pred = as.data.frame(predict(Inlet_pop, newdata = Inlet_preddata, interval = 'confidence'))
Inlet_preddata = data.frame(Inlet_preddata, post_fvfm = Inlet_pred$Prediction, Lower = Inlet_pred$Lower, Upper = Inlet_pred$Upper)


combo_plot <- ggplot() +
  geom_jitter(data = fvfm, aes(x = temp, y = post_fvfm, color = site), size = 1, width = 0.25) +
  scale_x_continuous(limits=c(26,40), breaks=c(27, 30, 33, 36, 39)) +
  scale_y_continuous(limits=c(-0.01, 0.85), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = Haikou_preddata, aes(x = temp, y = post_fvfm), color = 'seagreen', show.legend = FALSE) +
  geom_ribbon(data = Haikou_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'seagreen', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = Haikou_coeff_mean), color = 'seagreen', show.legend = FALSE, size =0.75) +
  
  geom_line(data = Haikou.shallow_preddata, aes(x = temp, y = post_fvfm), color = 'lightgreen', show.legend = FALSE) +
  geom_ribbon(data = Haikou.shallow_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'lightgreen', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = Haikou.shallow_coeff_mean), color = 'lightgreen', show.legend = FALSE, size =0.75) +
  
  geom_line(data = Outlet_preddata, aes(x = temp, y = post_fvfm), color = 'red', show.legend = FALSE) +
  geom_ribbon(data = Outlet_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'red', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = Outlet_coeff_mean), color = 'red', show.legend = FALSE, size =0.75) +
  
  geom_line(data = Inlet_preddata, aes(x = temp, y = post_fvfm), color = 'steelblue', show.legend = FALSE) +
  geom_ribbon(data = Inlet_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'steelblue', linetype=2, alpha = 0.2) +
  geom_vline(data = PA_coeff_means, aes(xintercept = Inlet_coeff_mean), color = 'steelblue', show.legend = FALSE, size =0.75) +
  
  ggtitle("") +
  scale_color_manual(values=c("seagreen", "lightgreen", "steelblue","red")) +
  ylab("Maximum quantum yield (Fv/Fm)") +
  xlab("Temperature (°C)") +
  theme_classic()

combo_plot +
  theme(text=element_text(size=12,  family="sans"))  

combo_plot +
  theme(text=element_text(size=16,  family="sans")) +
  theme(legend.position = "none") 

combo_plot +
  geom_text(data = PA_coeff_means, aes(label=round(Haikou_coeff_mean, digits = 2)), x = 31, y = 0.3, show.legend = FALSE, color = 'seagreen') +
  geom_text(data = PA_coeff_means, aes(label=round(Haikou.shallow_coeff_mean, digits = 2)), x = 31, y = 0.25, show.legend = FALSE, color = 'lightgreen') +
  geom_text(data = PA_coeff_means, aes(label=round(Outlet_coeff_mean, digits = 2)), x = 31, y = 0.20, show.legend = FALSE, color = 'red') +
  geom_text(data = PA_coeff_means, aes(label=round(Inlet_coeff_mean, digits = 2)), x = 31, y = 0.15, show.legend = FALSE, color = 'steelblue')


####4) Analysis####

#compare fvfm among sites before the experiment
pre_fvfm.lm <- lm(pre_fvfm ~ site, data = fvfm)
Anova(pre_fvfm.lm, type = "III")
emmeans(pre_fvfm.lm, pairwise ~ site, type="response")

mean(fvfm$pre_fvfm)
#0.7154167

sd(fvfm$pre_fvfm)
#0.04094393


#compare ED50s (ANOVA approach)

ED50_lm <- lm(ED50 ~ site, data = ED50s_long)
Anova(ED50_lm, type = "III")

#check assumptions
plot(fitted(ED50_lm),residuals(ED50_lm))
hist(residuals(ED50_lm))
qqnorm(residuals(ED50_lm))

#posthoc 
ED50.posthoc1 <- emmeans(ED50_lm, pairwise ~ site, type="response")
summary(ED50.posthoc1)



####mixed model assessment of Fv/Fm post experiment

fvfm.mm1 <- lmer(post_fvfm~treatment*site + (1|tank), data=fvfm)
summary(fvfm.mm1)
anova(fvfm.mm1, type = "III")


#check assumptions
plot(fitted(fvfm.mm1),residuals(fvfm.mm1))
hist(residuals(fvfm.mm1))
qqnorm(residuals(fvfm.mm1))

fvfm.mm1.posthoc1 <- emmeans(fvfm.mm1, pairwise~ site|treatment)
fvfm.mm1.posthoc1

fvfm.mm1.posthoc2 <- emmeans(fvfm.mm1, pairwise~ treatment|site)
fvfm.mm1.posthoc2


####*********************************####
####LARGE-AREA IMAGING####
####*********************************####

####<GROWTH>####

####1) Libraries####

library(tidyverse)
library(emmeans)
library(car)


####2) Data####

#long format (for models & boxplot)
growth_data_long <- read_csv("growth_data_long.csv") #no fused colonies included in this df

summary(growth_data_long)

growth_data_long$site <- as.factor(growth_data_long$site)
growth_data_long$plot_id <- as.factor(growth_data_long$plot_id)
growth_data_long$timepoint <- as.factor(growth_data_long$timepoint)

summary(growth_data_long)
head(growth_data_long)



#wide format for scatter plot and grew/shrunk counts
growth_data_wide <- read_csv("growth_data_wide.csv")

summary(growth_data_wide)

growth_data_wide$site <- as.factor(growth_data_wide$site)
growth_data_wide$plot_id <- as.factor(growth_data_wide$plot_id)
growth_data_wide$genet_id <- as.factor(growth_data_wide$genet_id)
growth_data_wide$status <- as.factor(growth_data_wide$status)

summary(growth_data_wide)
head(growth_data_wide)


#remove colonies that had fused together
growth_data_wide <-
  growth_data_wide%>%
  filter(include =="OK")

summary(growth_data_wide)

#quick look at change over time
growth_data_wide_sum <-
  growth_data_wide%>%
  group_by(site)%>%
  summarise(mean.change = mean(change), sd.change = sd(change))

growth_data_wide_sum

#if an error occurs while running the above code (e.g., not grouped properly by site) 
#it means that a previously loaded package is interfering 
#close and open R again to run properly


#site subsets

#Haikou reef
growth_data_wide.HK <-
  growth_data_wide%>%
  filter(site =="Haikou")

summary(growth_data_wide.HK)
#grew  :54
#shrank:8

54+8
#62

54/62
#0.8709677
#0.87

#Inlet reef
growth_data_wide.IT <-
  growth_data_wide%>%
  filter(site =="Inlet")

summary(growth_data_wide.IT)
#grew  :119
#shrank:46

119+46
#165

119/165
#0.7212121
#0.72

#Outlet reef
growth_data_wide.OT <-
  growth_data_wide%>%
  filter(site =="Outlet")

summary(growth_data_wide.OT)
#grew  :217
#shrank:71

217+71
#288

217/288
#0.7534722
#0.75



####3) Plots####

#boxplot

growth_data_long$timepoint <- factor(growth_data_long$timepoint, levels=c("spring_2022", "fall_2022"))
summary(growth_data_long)

size_plot <-
  ggplot(growth_data_long, aes(x=timepoint, y=size, fill=site, group=timepoint)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.5), alpha=0.4, size =1)+
  scale_fill_manual(values=c("seagreen", "steelblue","red"))+
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  labs(title="")+
  theme_bw()

size_plot + facet_wrap(~site) + theme(text=element_text(size=14,  family="sans")) 


#scatter plot
growth_plot <-
  ggplot(growth_data_wide, aes(T1_spring, T2_fall)) +
  geom_point(aes(color=status), size=2) +
  scale_colour_manual(values=c("black", "gray"))+
  geom_abline(slope=1, intercept=0)+
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  labs(y = "Fall 2022 (planar area cm2)", x="Spring 2022 (planar area cm2)")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_bw()

growth_plot 

growth_plot + facet_wrap(~site) + theme(legend.position = "none") + theme(text=element_text(size=16,  family="sans"))


#scatter plot (HK)
growth_plot.HK <-
  ggplot(growth_data_wide.HK, aes(T1_spring, T2_fall)) +
  geom_point(aes(color=status), size=2) +
  scale_colour_manual(values=c("black", "gray"))+
  geom_abline(slope=1, intercept=0)+
  labs(y = "Fall 2022", x="Spring 2022", title = "Haikou")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_bw()

growth_plot.HK

#scatter plot (OT)
growth_plot.OT <-
  ggplot(growth_data_wide.OT, aes(T1_spring, T2_fall)) +
  geom_point(aes(color=status), size=2) +
  scale_colour_manual(values=c("black", "gray"))+
  geom_abline(slope=1, intercept=0)+
  labs(y = "Fall 2022", x="Spring 2022", title = "Outlet")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_bw()

growth_plot.OT


#scatter plot
growth_plot.IT <-
  ggplot(growth_data_wide.IT, aes(T1_spring, T2_fall)) +
  geom_point(aes(color=status), size=2) +
  scale_colour_manual(values=c("black", "gray"))+
  geom_abline(slope=1, intercept=0)+
  labs(y = "Fall 2022", x="Spring 2022", title = "Inlet")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_bw()

growth_plot.IT




####4) Analysis####

#size lm 
size_lm <- lm(size ~ site*timepoint, data = growth_data_long)
Anova(size_lm, type = "III")

#check assumptions
plot(fitted(size_lm),residuals(size_lm))
hist(residuals(size_lm))
qqnorm(residuals(size_lm))


#use sqrt transformation to better meet assumptions
size_lm_sqrt <- lm(sqrt(size) ~ site*timepoint, data = growth_data_long)
Anova(size_lm, type = "III")

#check assumptions
plot(fitted(size_lm_sqrt),residuals(size_lm_sqrt))
hist(residuals(size_lm_sqrt))
qqnorm(residuals(size_lm_sqrt))

#posthoc 
size.posthoc1 <- emmeans(size_lm_sqrt, pairwise ~ site|timepoint, type="response")
summary(size.posthoc1)

size.posthoc2 <- emmeans(size_lm_sqrt, pairwise ~ timepoint|site, type="response")
summary(size.posthoc2)



####____________________#### 

####<LAI BLEACH>####


####1) Libraries####

library(tidyverse)
library(emmeans)
library(car)
library(chisq.posthoc.test)
library(patchwork)
library(corrplot)


####2) Data####

lai_bleach <- read_csv("lai_bleach.csv")

summary(lai_bleach)

lai_bleach$site <- as.factor(lai_bleach$site)
lai_bleach$plot_id <- as.factor(lai_bleach$plot_id)
lai_bleach$bleach_cat <- as.factor(lai_bleach$bleach_cat)

summary(lai_bleach)
head(lai_bleach)

bleach_percent <-
  lai_bleach%>%
  group_by(site, plot_id, bleach_cat)%>%
  summarise(n=n())

bleach_percent <-
  bleach_percent%>%
  group_by(plot_id)%>%
  mutate(sum = sum(n))  

bleach_percent <-
  bleach_percent%>%
  group_by(plot_id)%>%
  mutate(percent = n/sum*100) 

HK_subset.bleach <-
  bleach_percent%>%
  filter(site=='HK')

OT_subset.bleach <-
  bleach_percent%>%
  filter(site=='OT')

IT_subset.bleach <-
  bleach_percent%>%
  filter(site=='IT')


####3) Plots####

#score order
bleach_percent$bleach_cat <- factor(bleach_percent$bleach_cat, levels=c("4", "3", "2", "1"))

#stacked bar plot; all sites
lai_bleach_plot <-
  ggplot(bleach_percent, aes(fill=bleach_cat, y=percent, x=plot_id)) + 
  geom_bar(position="stack", stat="identity", colour="black")+
  scale_fill_manual(values=c("white", "bisque1","bisque3", "bisque4"))+
  theme_classic()+
  labs(y = "Percentage of colonies", x="Plot")+
  theme(text=element_text(size=16,  family="sans"))+
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100))

lai_bleach_plot

lai_bleach_plot + theme(legend.position="top")

#subset score order
HK_subset.bleach$bleach_cat <- factor(HK_subset.bleach$bleach_cat, levels=c("4", "3", "2", "1"))
IT_subset.bleach$bleach_cat <- factor(IT_subset.bleach$bleach_cat, levels=c("4", "3", "2", "1"))
OT_subset.bleach$bleach_cat <- factor(OT_subset.bleach$bleach_cat, levels=c("4", "3", "2", "1"))

#HK
lai_bleach_plot.HK <-
  ggplot(HK_subset.bleach, aes(fill=bleach_cat, y=percent, x=plot_id)) + 
  geom_bar(position="stack", stat="identity", colour="black")+
  scale_fill_manual(values=c("white", "bisque1","bisque3", "bisque4"))+
  theme_classic()+
  labs(y = "Percentage of colonies", x="Plot", title = "Haikou reef")+
  theme(text=element_text(size=16,  family="sans"))+
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100))

lai_bleach_plot.HK

lai_bleach_plot.HK + theme(legend.position="none")


#OT
lai_bleach_plot.OT <-
  ggplot(OT_subset.bleach, aes(fill=bleach_cat, y=percent, x=plot_id)) + 
  geom_bar(position="stack", stat="identity", colour="black")+
  scale_fill_manual(values=c("white", "bisque1","bisque3", "bisque4"))+
  theme_classic()+
  labs(y = "Percentage of colonies", x="Plot", title = "Outlet reef")+
  theme(text=element_text(size=16,  family="sans"))+
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100))

lai_bleach_plot.OT

lai_bleach_plot.OT + theme(legend.position="none")

#IT
lai_bleach_plot.IT <-
  ggplot(IT_subset.bleach, aes(fill=bleach_cat, y=percent, x=plot_id)) + 
  geom_bar(position="stack", stat="identity", colour="black")+
  scale_fill_manual(values=c("white", "bisque1","bisque3", "bisque4"))+
  theme_classic()+
  labs(y = "Percentage of colonies", x="Plot", , title = "Inlet reef")+
  theme(text=element_text(size=16,  family="sans"))+
  scale_y_continuous(limits=c(0, 100), breaks=c(0, 25, 50, 75, 100))

lai_bleach_plot.IT

lai_bleach_plot.IT + theme(legend.position="none")

LAI_BLEACH_COMBO <-
  lai_bleach_plot.HK  + theme(legend.position="none")| lai_bleach_plot.IT + theme(legend.position="none") | lai_bleach_plot.OT  + theme(legend.position="right") 

LAI_BLEACH_COMBO



####4) Analysis####

LAI_bleach_check <- table(lai_bleach$site, lai_bleach$bleach_cat)

LAI_bleach_check

mosaicplot(LAI_bleach_check)

prop.table(LAI_bleach_check, margin =1)

addmargins(LAI_bleach_check)

LAI_chisq <- chisq.test(LAI_bleach_check)
#LAI_chisq2 <- chisq.test(LAI_bleach_check,simulate.p.value = TRUE)
#gives same overall results without an error message

LAI_chisq

corrplot(LAI_chisq$residuals, is.cor = FALSE)

chisq.posthoc.test(LAI_bleach_check)


####____________________#### 

####<RUGOSITY>####


####1) Libraries####

library(tidyverse)
library(emmeans)
library(car)

####2) Data####

rugosity_data <- read_csv("rugosity_data.csv")

summary(rugosity_data)

rugosity_data$site <-as.factor(rugosity_data$site)
rugosity_data$plot_id <-as.factor(rugosity_data$plot_id)
rugosity_data$season <-as.factor(rugosity_data$season)
rugosity_data$scale <-as.factor(rugosity_data$scale)

summary(rugosity_data)
head(rugosity_data)

####3) Plots####

#box plot
rugosity_boxplot <- 
  ggplot(rugosity_data, aes(x=plot_id, y=rugosity, fill=site)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("seagreen", "steelblue","red"))+
  geom_point(position = position_jitterdodge(0.4), alpha=0.4, size =1)+
  scale_y_continuous(limits=c(0, 2.5), breaks=c(0, .5, 1, 1.5, 2, 2.5))+
  labs(y = "Rugosity", x="Plot")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_classic()

rugosity_boxplot 
rugosity_boxplot  + theme(legend.position = "none")
rugosity_boxplot  +  facet_wrap(~season) + theme_bw()


#bar plot

#first calculare mean, sd, se
rug_sum.plot <-
  rugosity_data%>%
  group_by(site, season, plot_id)%>%
  summarise(mean = mean(rugosity), sd=sd(rugosity), n=n(), se=sd/sqrt(n))

rug_sum.plot

rug_sum <-
  rugosity_data%>%
  group_by(site, season)%>%
  summarise(mean = mean(rugosity), sd=sd(rugosity), n=n(), se=sd/sqrt(n))

rug_sum

#by plot
rug_sum.plot$season <- factor(rug_sum.plot$season, levels=c("Spring", "Fall"))

rugosity_barplot.plot <- 
  ggplot(rug_sum.plot, aes(x=plot_id, y=mean, fill=site)) + 
  geom_bar(stat="identity", colour="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.3,position=position_dodge(.9))+ 
  scale_fill_manual(values=c("seagreen", "steelblue","red"))+
  scale_y_continuous(limits=c(0, 2), breaks=c(0, .5, 1, 1.5, 2))+
  labs(y = "Rugosity", x="Plot")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_bw()

rugosity_barplot.plot + facet_wrap(~season)

rugosity_barplot.plot  + facet_wrap(~season) + theme(legend.position = "top") +  theme(text=element_text(size=16,  family="sans"))

rugosity_barplot.plot  + facet_wrap(~season) + theme(legend.position = "none") +  theme(text=element_text(size=16,  family="sans"))

#by site
rugosity_barplot <- 
  ggplot(rug_sum, aes(x=site, y=mean, fill=site)) + 
  geom_bar(stat="identity", colour= "black" , position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), size= 0.7, width=.1, position=position_dodge(.9))+ 
  scale_fill_manual(values=c("seagreen", "steelblue","red"))+
  scale_y_continuous(limits=c(0, 2), breaks=c(0, .5, 1, 1.5, 2))+
  labs(y = "Rugosity", x="Site")+
  theme(text=element_text(size=16,  family="sans"))+
  theme_bw()

rugosity_barplot + facet_wrap(~season) + theme(legend.position = "bottom")


####4) Analysis####

plot(rugosity_data$rugosity)
hist(rugosity_data$rugosity)

rugosity.lm <- lm(rugosity ~ site * season, data = rugosity_data)
Anova(rugosity.lm, type = "III")

#check assumptions
plot(fitted(rugosity.lm),residuals(rugosity.lm))
hist(residuals(rugosity.lm))
qqnorm(residuals(rugosity.lm))
vif(rugosity.lm)

#posthoc 
rugosity.posthoc1 <- emmeans(rugosity.lm, pairwise ~ site|season, type="response")
summary(rugosity.posthoc1)

rugosity.posthoc2 <- emmeans(rugosity.lm, pairwise ~ season|site, type="response")
summary(rugosity.posthoc2)


####____________________#### 

####<BENTHIC COMPOSTION>####


####1) Libraries####

library(tidyverse)
library(ggfortify)
library(vegan)
library(lme4)
library(performance)
library(patchwork)
library(car)
library(emmeans)

####2) Data####

benthic_comp <- read_csv("benthic_comp.csv")

summary(benthic_comp)

benthic_comp$site <- as.factor(benthic_comp$site)
benthic_comp$plot_id <- as.factor(benthic_comp$plot_id)
benthic_comp$timepoint <- as.factor(benthic_comp$timepoint)
benthic_comp$cover_type <- as.factor(benthic_comp$cover_type)

summary(benthic_comp)
head(benthic_comp)

#spring subset
benthic_comp_spring <-
  benthic_comp%>%
  filter(timepoint =="spring_2022")

benthic_comp_spring

#percent cover (mean and sd)
benthic_comp_spring_sum <-
  benthic_comp_spring%>%
  group_by(site, cover_type)%>%
  summarise(mean_spring = mean(percentage_cover), sd_sping = sd (percentage_cover))

benthic_comp_spring_sum



#fall subset
benthic_comp_fall <-
  benthic_comp%>%
  filter(timepoint =="fall_2022")

benthic_comp_fall

#percent cover (mean and sd)
benthic_comp_fall_sum <-
  benthic_comp_fall%>%
  group_by(site, cover_type)%>%
  summarise(mean_fall = mean(percentage_cover), sd_fall = sd (percentage_cover))

benthic_comp_fall_sum



####2) Plots####

#set orders
benthic_comp_spring$cover_type <- 
  factor(benthic_comp_spring$cover_type, 
         levels=c("other", "sand_sediment", "sponge", "turf", "fleshy_macro", "calcified_macro", "CCA",  "soft_coral", "hard_coral"))

benthic_comp_spring$site <- 
  factor(benthic_comp_spring$site, 
         levels=c("HK", "IT", "OT"))

benthic_comp_fall$cover_type <- 
  factor(benthic_comp_fall$cover_type, 
         levels=c("other", "sand_sediment", "sponge", "turf", "fleshy_macro", "calcified_macro", "CCA",  "soft_coral", "hard_coral"))

benthic_comp_fall$site <- 
  factor(benthic_comp_fall$site, 
         levels=c("HK", "IT", "OT"))



#spring lolipop plot
cover_spring_plot <-
  ggplot(benthic_comp_spring, aes(x=cover_type, y=percentage_cover, colour= site)) +
  geom_segment( aes(x=cover_type, xend=cover_type, y=0, yend=percentage_cover)) +
  geom_point(size=2.5, alpha=0.8) +
  coord_flip() +
  scale_colour_manual(values=c("seagreen", "steelblue", "red"))+
  labs(title="Spring 2022",y="Percent cover", x = "Benthic category")+
  theme_bw()+
  facet_wrap(~site)

cover_spring_plot 
COVER_FALL <-cover_spring_plot  + facet_wrap(~plot_id) + theme(text=element_text(size=16,  family="sans")) + theme(legend.position="none")

COVER_FALL

#fall lolipop plot
cover_fall_plot <-
  ggplot(benthic_comp_fall, aes(x=cover_type, y=percentage_cover, colour= site)) +
  geom_segment( aes(x=cover_type, xend=cover_type, y=0, yend=percentage_cover)) +
  geom_point(size=2.5, alpha=0.8) +
  coord_flip() +
  scale_colour_manual(values=c("seagreen", "steelblue", "red"))+
  labs(title="B. Fall 2022",y="Percent cover", x = "Benthic category")+
  theme_bw()+
  facet_wrap(~site)

cover_fall_plot 
COVER_SPRING <- cover_fall_plot + facet_wrap(~plot_id) + theme(text=element_text(size=16,  family="sans")) + theme(legend.position="none")

COVER_SPRING

COVER_COMBO <- COVER_FALL + COVER_SPRING

COVER_COMBO



####*PCA: spring####

#change format from long to wide
benthic_comp_spring_wide <- spread(benthic_comp_spring, cover_type, percentage_cover)
benthic_comp_spring_wide 

#remove NAs
benthic_comp_spring_wide <- benthic_comp_spring_wide[complete.cases(benthic_comp_spring_wide), ]

#check the correlation matrix of explanatory variables
corr_mat <- as.matrix(round(cor(benthic_comp_spring_wide[, c(5:13)]), 2))

#remove top half for ease of interpretation
corr_mat[upper.tri(corr_mat)] <- NA
corr_mat

#plot by site
pca_values <- 
  prcomp(benthic_comp_spring_wide[, c(5:13)],
         center = TRUE, 
         scale = TRUE)

summary(pca_values)
str(pca_values)

#simple PCA plot
autoplot(pca_values, 
         loadings = TRUE,
         loadings.label = TRUE) 

#add these PCA values to our original dataframe
pca_points <- 
  data.frame(benthic_comp_spring_wide, pca_values$x)

head(pca_points)

#pull out the eigenvectors (the arrows showing the loadings)
pca_load <- 
  data.frame(variables = rownames(pca_values$rotation), pca_values$rotation)

pca_load

#create a convex hull - the smallest polygon that includes all the points of a given level
pca_hull <- 
  pca_points %>% 
  group_by(site) %>% 
  slice(chull(PC1, PC2))

#spring PCA plot
pca_plot <- 
  ggplot(pca_points, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = site),
             size = 2, alpha = 0.7) +
  geom_polygon(data = pca_hull,
               aes(fill = site),
               alpha = 0.3,
               show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = (pca_load$PC1*5.7), y = (pca_load$PC2*5.2),
           label = pca_load$variables,
           size = 2.5) +
  theme_classic() 

pca_plot

col_list<-c("seagreen", "steelblue", "red") 

SPRING_PCA <-
  pca_plot +
  scale_color_manual(values = col_list) +
  scale_fill_manual(values = col_list) + theme(text=element_text(size=16,  family="sans")) +
  labs(title="Spring 2022",y="PC2 (21.33%)", x = "PC1 (42.68%)")

SPRING_PCA

SPRING_PCA <-
  pca_plot +
  scale_color_manual(values = col_list) +
  scale_fill_manual(values = col_list) + theme(text=element_text(size=16,  family="sans")) +
  labs(title="Spring 2022",y="PC2 (21.33%)", x = "PC1 (42.68%)")+ theme(legend.position = "none")

SPRING_PCA


####*PCA: fall####

#change format from long to wide
benthic_comp_fall_wide <- spread(benthic_comp_fall, cover_type, percentage_cover)
benthic_comp_fall_wide 

#remove NAs
benthic_comp_fall_wide <- benthic_comp_fall_wide[complete.cases(benthic_comp_fall_wide), ]

#check the correlation matrix of explanatory variables
corr_mat <- as.matrix(round(cor(benthic_comp_fall_wide[, c(5:13)]), 2))

#remove top half for ease of interpretation
corr_mat[upper.tri(corr_mat)] <- NA
corr_mat

#plot by site
pca_values <- 
  prcomp(benthic_comp_fall_wide[, c(5:13)],
         center = TRUE, 
         scale = TRUE)

summary(pca_values)
str(pca_values)

#quick PCA plot
autoplot(pca_values, 
         loadings = TRUE,
         loadings.label = TRUE) 

#add these PCA values to our original dataframe
pca_points <- 
  data.frame(benthic_comp_fall_wide, pca_values$x)

head(pca_points)


#pull out the eigenvectors (the arrows showing the loadings)
pca_load <- 
  data.frame(variables = rownames(pca_values$rotation), pca_values$rotation)

pca_load

#create a convex hull - the smallest polygon that includes all the points of a given level
pca_hull <- 
  pca_points %>% 
  group_by(site) %>% 
  slice(chull(PC1, PC2))

#fall PCA plot
pca_plot <- 
  ggplot(pca_points, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = site),
             size = 2, alpha = 0.7) +
  geom_polygon(data = pca_hull,
               aes(fill = site),
               alpha = 0.3,
               show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = (pca_load$PC1*5.7), y = (pca_load$PC2*5.2),
           label = pca_load$variables,
           size = 2.5) +
  theme_classic() 

pca_plot

col_list<-c("seagreen", "steelblue", "red") 

FALL_PCA <-
  pca_plot +
  scale_color_manual(values = col_list) +
  scale_fill_manual(values = col_list) + theme(text=element_text(size=16,  family="sans")) +
  labs(title="Fall 2022",y="PC2 (23.92%)", x = "PC1 (36.53%)")

FALL_PCA

FALL_PCA <-
  pca_plot +
  scale_color_manual(values = col_list) +
  scale_fill_manual(values = col_list) + theme(text=element_text(size=16,  family="sans")) +
  labs(title="Fall 2022",y="PC2 (23.92%)", x = "PC1 (36.53%)")+ theme(legend.position = "none")

FALL_PCA


####4) Analysis####

####*Spring#### 
#IMPORTANT: re-run the spring PCA plot code before running the analysis below 
#(this is due to similar naming in the PCA plot code for spring and fall time points)

pca_load

summary(pca_values)

metadata <- benthic_comp_spring_wide[, c(1)] #the explanatory variables we want to test
test <- envfit(pca_values, metadata, permu= 999)
test

pca_points

spring_comp_lm <- lm(PC1 ~ site, data = pca_points)
summary(spring_comp_lm)
Anova(spring_comp_lm, type = "III")

#check assumptions
plot(fitted(spring_comp_lm),residuals(spring_comp_lm))
hist(residuals(spring_comp_lm))
qqnorm(residuals(spring_comp_lm))
vif(spring_comp_lm)

#posthoc 
spring_comp.posthoc1 <- emmeans(spring_comp_lm, pairwise ~ site, type="response")
summary(spring_comp.posthoc1)


####*Fall####
#IMPORTANT: re-run the fall PCA plot code before running the analysis below 
#(this is due to similar naming in the PCA plot code for spring and fall time points)

pca_load

summary(pca_values)

metadata <- benthic_comp_fall_wide[, c(1)] #the explanatory variables we want to test
test <- envfit(pca_values, metadata, permu= 999)
test

fall_comp_lm <- lm(PC1 ~ site, data = pca_points)
summary(fall_comp_lm)
Anova(fall_comp_lm, type = "III")

#check assumptions
plot(fitted(fall_comp_lm),residuals(fall_comp_lm))
hist(residuals(fall_comp_lm))
qqnorm(residuals(fall_comp_lm))

#posthoc 
fall_comp.posthoc1 <- emmeans(fall_comp_lm, pairwise ~ site, type="response")
summary(fall_comp.posthoc1)


####__________END__________#### 


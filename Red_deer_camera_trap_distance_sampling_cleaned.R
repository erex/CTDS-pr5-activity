############################################ Camera trap distance  sampling ########################################################################

# Photos and deployment days bootstrapped per camera
# Delay between triggers (with SE)
# Detection probability overall (with SE)
# Activity per month (with SE)
# Group size per season (summer and winter), bootstrapped
# Angle constant
# Area constant

setwd("D:/Users/maik-/Documents/Maik/Rothirschprojekt Bayerischer Wald/Distance Sampling/Distance Sampling mit Kamerafallen")

library(plyr)
library(dplyr)
library(lubridate)
library(reshape)
library(boot)
library(activity)
library(Distance)
library(ggplot2)
library(car)
library(reshape)

########## Deployments, record table, and group size#################################################################################

# Deployment days
Deployment_table<-read.csv2("D:/Users/maik-/Documents/Maik/Rothirschprojekt Bayerischer Wald/Interreg Camera trapping/Deployments/Deployment_matrix_gesamt.csv")
Camera_operations<-data.frame(Station=c(Deployment_table$X, "Total"), 
                              Mai_2018=c(rowSums(Deployment_table[,2:26], na.rm = T),NA), 
                              Jun_2018=c(rowSums(Deployment_table[,27:56], na.rm=T),NA), 
                              Jul_2018=c(rowSums(Deployment_table[,57:87], na.rm=T),NA), 
                              Aug_2018=c(rowSums(Deployment_table[,88:118], na.rm=T),NA), 
                              Sep_2018=c(rowSums(Deployment_table[,119:148], na.rm=T),NA),
                              Okt_2018=c(rowSums(Deployment_table[,149:179], na.rm=T),NA),
                              Nov_2018=c(rowSums(Deployment_table[,180:209], na.rm=T),NA),
                              Dez_2018=c(rowSums(Deployment_table[,210:240], na.rm=T),NA),
                              Jan_2019=c(rowSums(Deployment_table[,241:271], na.rm=T),NA),
                              Feb_2019=c(rowSums(Deployment_table[,272:299], na.rm=T),NA),
                              Mrz_2019=c(rowSums(Deployment_table[,300:330], na.rm=T),NA),
                              Apr_2019=c(rowSums(Deployment_table[,331:360], na.rm=T),NA),
                              Mai_2019=c(rowSums(Deployment_table[,361:391], na.rm=T),NA),
                              Jun_2019=c(rowSums(Deployment_table[,392:421], na.rm=T),NA),
                              Jul_2019=c(rowSums(Deployment_table[,422:452], na.rm=T),NA)) 
Camera_operations[109,2:16]<-colSums(Camera_operations[1:108,2:16])
Camera_operations2<-melt(Camera_operations[1:108,], id="Station")
names(Camera_operations2)[1]<-"Station"
names(Camera_operations2)[2]<-"month_year"
names(Camera_operations2)[3]<-"Cam_days"
Camera_operations2$Cam_secs<-Camera_operations2$Cam_days*86400
Camera_operations2$Cam_days<-NULL

names(Deployment_table)[1]<-"Station"
Deployment_table2_BFNP<-melt(Deployment_table[,1:488], id="Station")
Deployment_table2_BFNP<-Deployment_table2_BFNP%>%arrange(Station,variable)
Deployment_table2_BFNP$Depl_IDs<-paste(Deployment_table2_BFNP$Station, substr(Deployment_table2_BFNP$variable,2,11))
Deployment_table2_BFNP_active<-subset(Deployment_table2_BFNP,Deployment_table2_BFNP$value==1)

BFNP<-read.csv2("D:/Users/maik-/Documents/Maik/Rothirschprojekt Bayerischer Wald/Interreg Camera trapping/Record tables/BFNP/Record_table_BFNP_29072020_Events_edited_new.csv")
BFNP$Grid_cell_number[is.na(BFNP$Grid_cell_number)==T]<-BFNP$Station[is.na(BFNP$Grid_cell_number)==T]
BFNP$Depl_IDs<-paste(BFNP$Grid_cell_number, BFNP$Date)
BFNP<-subset(BFNP,BFNP$Depl_IDs%in%Deployment_table2_BFNP_active$Depl_IDs)
BFNP$delta.time.mins<-as.numeric(as.character(BFNP$delta.time.mins))
BFNP$delta.time.hours<-as.numeric(as.character(BFNP$delta.time.hours))
BFNP$delta.time.days<-as.numeric(as.character(BFNP$delta.time.days))
BFNP$Date<-as.POSIXct(BFNP$Date, format="%d.%m.%Y")
BFNP$month<-as.factor(format.Date(BFNP$Date, format="%h"))
BFNP$DateTimeOriginal<-as.POSIXct(paste(BFNP$Date, BFNP$Time), format="%Y-%m-%d %H:%M:%S")
BFNP$month_year<-as.factor(format.Date(BFNP$Date, format="%h_%Y"))

BFNP_red_deer<-subset(BFNP,BFNP$Species=="Red_deer")
BFNP_photos_series_location_full<-aggregate(Species ~ Station +  month_year, data=BFNP_red_deer, function(x)length(x))
names(BFNP_photos_series_location_full)[3]<-"Nr_photo_series"

Red_deer_times_photos_BFNP<-left_join(Camera_operations2, BFNP_photos_series_location_full)
Red_deer_times_photos_BFNP$Nr_photo_series[is.na(Red_deer_times_photos_BFNP$Nr_photo_series)==T]<-0
rm(BFNP, BFNP_photos_series_location_full, Deployment_table, Camera_operations, Camera_operations2)

Red_deer_times_photos_BFNP$month<-factor(substr(Red_deer_times_photos_BFNP$month_year,1,3))
Red_deer_times_photos_BFNP$Season<-ifelse(Red_deer_times_photos_BFNP$month%in%c("Mai","Jun","Jul","Aug","Sep","Okt"),"Summer","Winter")



# Group size
BFNP_red_deer$Season<-ifelse(BFNP_red_deer$month%in%c("Mai","Jun","Jul","Aug","Sep","Okt"),"Summer","Winter")
BFNP_red_deer$Number_of_animals_species_1<-as.numeric(BFNP_red_deer$Number_of_animals_species_1)

set.seed(12345)
boot_grp_size_winter<-boot(data=BFNP_red_deer[BFNP_red_deer$Season=="Winter",], statistic = Grp_size_est, R=10000)
set.seed(12345)
boot_grp_size_summer<-boot(data=BFNP_red_deer[BFNP_red_deer$Season=="Summer",], statistic = Grp_size_est, R=10000)
Group_size<-data.frame(Season=rep(c("Winter","Summer"),each=10000),Group_size=c(boot_grp_size_winter$t,boot_grp_size_summer$t),Group_size_t0=rep(c(boot_grp_size_winter$t0, boot_grp_size_summer$t0),each=10000)) 




# Activity
BFNP_red_deer$radian_time<-gettime(BFNP_red_deer$DateTimeOriginal,  format="%d.%m.%Y %H:%M", scale ="radian")
BFNP_red_deer$hour<-gettime(BFNP_red_deer$DateTimeOriginal,  format="%d.%m.%Y %H:%M", scale ="hour")

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mai"])
Red_deer_BFNP_activity_Mai<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mai"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Mai)
Red_deer_BFNP_activity_Mai_df<-data.frame(month="Mai", activity=Red_deer_BFNP_activity_Mai@act[1], se=Red_deer_BFNP_activity_Mai@act[2], Lower_95_CI=Red_deer_BFNP_activity_Mai@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Mai@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jun"])
Red_deer_BFNP_activity_Jun<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jun"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Jun)
Red_deer_BFNP_activity_Jun_df<-data.frame(month="Jun", activity=Red_deer_BFNP_activity_Jun@act[1], se=Red_deer_BFNP_activity_Jun@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jun@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jun@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jul"])
Red_deer_BFNP_activity_Jul<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jul"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Jul)
Red_deer_BFNP_activity_Jul_df<-data.frame(month="Jul", activity=Red_deer_BFNP_activity_Jul@act[1], se=Red_deer_BFNP_activity_Jul@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jul@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jul@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Aug"])
Red_deer_BFNP_activity_Aug<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Aug"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Aug)
Red_deer_BFNP_activity_Aug_df<-data.frame(month="Aug", activity=Red_deer_BFNP_activity_Aug@act[1], se=Red_deer_BFNP_activity_Aug@act[2], Lower_95_CI=Red_deer_BFNP_activity_Aug@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Aug@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Sep"])
Red_deer_BFNP_activity_Sep<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Sep"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Sep)
Red_deer_BFNP_activity_Sep_df<-data.frame(month="Sep", activity=Red_deer_BFNP_activity_Sep@act[1], se=Red_deer_BFNP_activity_Sep@act[2], Lower_95_CI=Red_deer_BFNP_activity_Sep@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Sep@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Okt"])
Red_deer_BFNP_activity_Okt<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Okt"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Okt)
Red_deer_BFNP_activity_Okt_df<-data.frame(month="Okt", activity=Red_deer_BFNP_activity_Okt@act[1], se=Red_deer_BFNP_activity_Okt@act[2], Lower_95_CI=Red_deer_BFNP_activity_Okt@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Okt@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Nov"])
Red_deer_BFNP_activity_Nov<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Nov"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Nov)
Red_deer_BFNP_activity_Nov_df<-data.frame(month="Nov", activity=Red_deer_BFNP_activity_Nov@act[1], se=Red_deer_BFNP_activity_Nov@act[2], Lower_95_CI=Red_deer_BFNP_activity_Nov@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Nov@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Dez"])
Red_deer_BFNP_activity_Dez<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Dez"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Dez)
Red_deer_BFNP_activity_Dez_df<-data.frame(month="Dez", activity=Red_deer_BFNP_activity_Dez@act[1], se=Red_deer_BFNP_activity_Dez@act[2], Lower_95_CI=Red_deer_BFNP_activity_Dez@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Dez@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jan"])
Red_deer_BFNP_activity_Jan<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jan"], reps=10000, sample="data", show=TRUE)
plot(Red_deer_BFNP_activity_Jan)
Red_deer_BFNP_activity_Jan_df<-data.frame(month="Jan", activity=Red_deer_BFNP_activity_Jan@act[1], se=Red_deer_BFNP_activity_Jan@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jan@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jan@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Feb"])
Red_deer_BFNP_activity_Feb<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Feb"], reps=10000, sample="data", show=TRUE)
plot(Red_deer_BFNP_activity_Feb)
Red_deer_BFNP_activity_Feb_df<-data.frame(month="Feb", activity=Red_deer_BFNP_activity_Feb@act[1], se=Red_deer_BFNP_activity_Feb@act[2], Lower_95_CI=Red_deer_BFNP_activity_Feb@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Feb@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mrz"])
Red_deer_BFNP_activity_Mrz<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mrz"], reps=10000, sample="data", show=TRUE)
plot(Red_deer_BFNP_activity_Mrz)
Red_deer_BFNP_activity_Mrz_df<-data.frame(month="Mrz", activity=Red_deer_BFNP_activity_Mrz@act[1], se=Red_deer_BFNP_activity_Mrz@act[2], Lower_95_CI=Red_deer_BFNP_activity_Mrz@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Mrz@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Apr"])
Red_deer_BFNP_activity_Apr<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Apr"], reps=10000, sample="model", show=TRUE)
plot(Red_deer_BFNP_activity_Apr)
Red_deer_BFNP_activity_Apr_df<-data.frame(month="Apr", activity=Red_deer_BFNP_activity_Apr@act[1], se=Red_deer_BFNP_activity_Apr@act[2], Lower_95_CI=Red_deer_BFNP_activity_Apr@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Apr@act[4])

Red_deer_BFNP_activity_df<-rbind(Red_deer_BFNP_activity_Mai_df,Red_deer_BFNP_activity_Jun_df, Red_deer_BFNP_activity_Jul_df, Red_deer_BFNP_activity_Aug_df, Red_deer_BFNP_activity_Sep_df,  Red_deer_BFNP_activity_Okt_df, Red_deer_BFNP_activity_Nov_df, Red_deer_BFNP_activity_Dez_df, Red_deer_BFNP_activity_Jan_df, Red_deer_BFNP_activity_Feb_df, Red_deer_BFNP_activity_Mrz_df, Red_deer_BFNP_activity_Apr_df)

rm(Red_deer_BFNP_activity_Mai_df,Red_deer_BFNP_activity_Jun_df, Red_deer_BFNP_activity_Jul_df, Red_deer_BFNP_activity_Aug_df, Red_deer_BFNP_activity_Sep_df,  Red_deer_BFNP_activity_Okt_df, Red_deer_BFNP_activity_Nov_df, Red_deer_BFNP_activity_Dez_df, Red_deer_BFNP_activity_Jan_df, Red_deer_BFNP_activity_Feb_df, Red_deer_BFNP_activity_Mrz_df, Red_deer_BFNP_activity_Apr_df)

rm(BFNP_red_deer, boot_grp_size_summer, boot_grp_size_winter)

########## Delays between triggers ########################################################################################
load("Kameraverz?gerung.RData")
summary(rec_test_photos$delta.time.secs)
hist(rec_test_photos$delta.time.secs, xlim=c(1,20),breaks=110)
hist(rec_test_photos$delta.time.secs, xlim=c(1,50),breaks=110)

# Option 1: Draw from normal distribution
set.seed(12345)
Camera_lag<-rnorm(10000, mean=6, sd=(1.5/1.96))
hist(Camera_lag)

# Option 2: Draw from original distribution with replacement:
set.seed(12345)
Camera_lag_2<-sample(subset(rec_test_photos$delta.time.sec, rec_test_photos$delta.time.sec<20),10000, replace=TRUE)
hist(Camera_lag_2)
summary(Camera_lag_2)
Camera_lag_3<-sample(subset(rec_test_photos$delta.time.sec, rec_test_photos$delta.time.sec<=9),10000, replace=TRUE)
hist(rec_test_photos$delta.time.secs, xlim=c(3,9),breaks=220)


########################################## Detection probability ##########################################################

Red_deer_detection_area_table<-read.csv2("D:/Users/maik-/Documents/Maik/Rothirschprojekt Bayerischer Wald/Interreg Camera trapping/Detection zones/Measurement_table_red_deer_BFNP_complete_10_station_month_csv.csv") 
Red_deer_detection_area_table<-subset(Red_deer_detection_area_table, Red_deer_detection_area_table$Y_m!="")
Red_deer_detection_area_table<-subset(Red_deer_detection_area_table, Red_deer_detection_area_table$Species=="Red_deer")
Red_deer_detection_area_table$month_year<-paste(Red_deer_detection_area_table$month, format.Date(as.POSIXct(as.character(Red_deer_detection_area_table$DateTimeOriginal),format="%d.%m.%Y %H:%M"), "%Y"),sep='-')
Red_deer_detection_area_table$Season<-as.factor(ifelse(Red_deer_detection_area_table$month%in%c("Mai","Jun","Jul","Aug","Sep","Okt"),"Summer","Winter")) 

Red_deer_total_radius<-data.frame(siteID=Red_deer_detection_area_table$Station,distance=Red_deer_detection_area_table$radius_m, Species=Red_deer_detection_area_table$Species) 
Red_deer_total_radius$distance<-as.numeric(as.character(Red_deer_total_radius$distance)) 
Red_deer_total_radius$siteID<-as.factor(Red_deer_total_radius$siteID)

# Proportion of observations within the truncation distance of 14 m
Prop_within_trunc<-length(Red_deer_total_radius$dist[Red_deer_total_radius$dist<15])/length(Red_deer_total_radius$dist)

# Fitting the detection function and extracting the probability of detection
Detfunc1<-ds(Red_deer_total_radius, truncation=14, transect="point",formula=~1, key="hn")
Detfunc2<-ds(Red_deer_total_radius, truncation=14, transect="point",formula=~1, key="hr")

summarize_ds_models(Detfunc1, Detfunc2)
# half normal detection function with expansions

summary(Detfunc1) # Distance package (in contrast to the Rdistance-package gives the detection probability with a standard error, but not the effective detection radius)
gof_ds(Detfunc1, plot=TRUE) # model fits well
plot(Detfunc1)

Det_prop<-0.2214239 # copy from the summary of the detection function
Det_prop_SE<-0.01252955


# Angle
theta<-(55*pi)/180

# Area in km?
Area<- 230



############################################################# Function #############################################################################################

# use camera trap delay as snapshot moment interval!

Dist_red_deer_months<-Dist_conf_function_months_base(Camera_lag_3,Red_deer_times_photos_BFNP, Det_prop, Det_prop_SE, trunc_dist=14, Prop_within_trunc, angle=55,Red_deer_BFNP_activity_df, Group_size, 2000, 0.05 )

Dist_red_deer_months_plotting<-Dist_conf_function_months_plotting(Camera_lag_3, Red_deer_times_photos_BFNP, Det_prop, Det_prop_SE, trunc_dist=14, Prop_within_trunc, angle=55,Red_deer_BFNP_activity_df, Group_size, 2000)

Dist_red_deer_months_plotting$month<-factor(Dist_red_deer_months_plotting$month, levels=c("Mai_2018","Jun_2018","Jul_2018","Aug_2018","Sep_2018","Okt_2018", "Nov_2018", "Dez_2018","Jan_2019", "Feb_2019", "Mrz_2019", "Apr_2019", "Mai_2019","Jun_2019", "Jul_2019"),ordered = T)


Fig1<-ggplot(Dist_red_deer_months_plotting,aes(x=month, y=Range))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x=element_text(size=12, vjust=0, margin=margin(1,0,0,0,"cm")),axis.title.y=element_text(size=12, vjust=1,margin=margin(0,1,0,0,"cm")), 
        axis.text.x=element_text(size=12, colour="black",angle=90, vjust=0.5),axis.text.y=element_text(size=12, colour="black"),
        plot.title=element_text(size=14,vjust=0,hjust=0.5, margin=margin(0,0,1,0,"cm")),
        plot.margin = margin(1,2,1,2,"cm"))+ 
  ylim(0,8)+
  xlab("Month")+
  ylab("Population density [Animals/km?]")
Fig1



Dist_red_deer_months_test<-Dist_conf_function_months_base(rep(6,10000),Red_deer_times_photos_BFNP, Det_prop, Det_prop_SE, trunc_dist=14, Prop_within_trunc, angle=55,Red_deer_BFNP_activity_df, Group_size, 2000, 0.05 )

Dist_red_deer_months_plotting_test<-Dist_conf_function_months_plotting(rep(6,10000), Red_deer_times_photos_BFNP, Det_prop, Det_prop_SE, trunc_dist=14, Prop_within_trunc, angle=55,Red_deer_BFNP_activity_df, Group_size, 2000)
Dist_red_deer_months_plotting_test$month<-factor(Dist_red_deer_months_plotting$month, levels=c("Mai_2018","Jun_2018","Jul_2018","Aug_2018","Sep_2018","Okt_2018", "Nov_2018", "Dez_2018","Jan_2019", "Feb_2019", "Mrz_2019", "Apr_2019", "Mai_2019","Jun_2019", "Jul_2019"),ordered = T)

Fig2<-ggplot(Dist_red_deer_months_plotting_test,aes(x=month, y=Range))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x=element_text(size=12, vjust=0, margin=margin(1,0,0,0,"cm")),axis.title.y=element_text(size=12, vjust=1,margin=margin(0,1,0,0,"cm")), 
        axis.text.x=element_text(size=12, colour="black",angle=90, vjust=0.5),axis.text.y=element_text(size=12, colour="black"),
        plot.title=element_text(size=14,vjust=0,hjust=0.5, margin=margin(0,0,1,0,"cm")),
        plot.margin = margin(1,2,1,2,"cm"))+ 
  ylim(0,8)+
  xlab("Month")+
  ylab("Population density [Animals/km?]")
Fig2

write.csv(Dist_red_deer_months, file="Density_estimates_distances_red_deer_months_new.csv")
write.csv(Dist_red_deer_months_test, file="Density_estimates_distances_red_deer_months_new_wo_variance_in_trigger_speed.csv")


#################################### Add-on: Activity data per sex  #############################################################################

(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mai" & is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Mai_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mai"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Mai_females_df_females<-data.frame(month="Mai", activity=Red_deer_BFNP_activity_Mai_females@act[1], se=Red_deer_BFNP_activity_Mai_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Mai_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Mai_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jun"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Jun_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jun"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Jun_females_df_females<-data.frame(month="Jun", activity=Red_deer_BFNP_activity_Jun_females@act[1], se=Red_deer_BFNP_activity_Jun_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jun_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jun_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jul"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Jul_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jul"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Jul_females_df_females<-data.frame(month="Jul", activity=Red_deer_BFNP_activity_Jul_females@act[1], se=Red_deer_BFNP_activity_Jul_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jul_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jul_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Aug"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Aug_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Aug"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Aug_females_df_females<-data.frame(month="Aug", activity=Red_deer_BFNP_activity_Aug_females@act[1], se=Red_deer_BFNP_activity_Aug_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Aug_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Aug_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Sep"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Sep_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Sep"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Sep_females_df_females<-data.frame(month="Sep", activity=Red_deer_BFNP_activity_Sep_females@act[1], se=Red_deer_BFNP_activity_Sep_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Sep_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Sep_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Okt"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Okt_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Okt"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Okt_females_df_females<-data.frame(month="Okt", activity=Red_deer_BFNP_activity_Okt_females@act[1], se=Red_deer_BFNP_activity_Okt_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Okt_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Okt_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Nov"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Nov_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Nov"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Nov_females_df_females<-data.frame(month="Nov", activity=Red_deer_BFNP_activity_Nov_females@act[1], se=Red_deer_BFNP_activity_Nov_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Nov_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Nov_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Dez"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Dez_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Dez"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Dez_females_df_females<-data.frame(month="Dez", activity=Red_deer_BFNP_activity_Dez_females@act[1], se=Red_deer_BFNP_activity_Dez_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Dez_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Dez_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jan"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Jan_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jan"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Jan_females_df_females<-data.frame(month="Jan", activity=Red_deer_BFNP_activity_Jan_females@act[1], se=Red_deer_BFNP_activity_Jan_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jan_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jan_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Feb"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Feb_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Feb"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Feb_females_df_females<-data.frame(month="Feb", activity=Red_deer_BFNP_activity_Feb_females@act[1], se=Red_deer_BFNP_activity_Feb_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Feb_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Feb_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mrz"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Mrz_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mrz"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Mrz_females_df_females<-data.frame(month="Mrz", activity=Red_deer_BFNP_activity_Mrz_females@act[1], se=Red_deer_BFNP_activity_Mrz_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Mrz_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Mrz_females@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Apr"& is.na(BFNP_red_deer$Number_females)==F])
Red_deer_BFNP_activity_Apr_females<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Apr"& is.na(BFNP_red_deer$Number_females)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Apr_females_df_females<-data.frame(month="Apr", activity=Red_deer_BFNP_activity_Apr_females@act[1], se=Red_deer_BFNP_activity_Apr_females@act[2], Lower_95_CI=Red_deer_BFNP_activity_Apr_females@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Apr_females@act[4])


Red_deer_BFNP_activity_df_females<-rbind(Red_deer_BFNP_activity_Mai_females_df_females, Red_deer_BFNP_activity_Jun_females_df_females, Red_deer_BFNP_activity_Jul_females_df_females, Red_deer_BFNP_activity_Aug_females_df_females, Red_deer_BFNP_activity_Sep_females_df_females,  Red_deer_BFNP_activity_Okt_females_df_females, Red_deer_BFNP_activity_Nov_females_df_females, Red_deer_BFNP_activity_Dez_females_df_females, Red_deer_BFNP_activity_Jan_females_df_females, Red_deer_BFNP_activity_Feb_females_df_females, Red_deer_BFNP_activity_Mrz_females_df_females, Red_deer_BFNP_activity_Apr_females_df_females)

rm(Red_deer_BFNP_activity_Mai_females_df_females, Red_deer_BFNP_activity_Jun_females_df_females, Red_deer_BFNP_activity_Jul_females_df_females, Red_deer_BFNP_activity_Aug_females_df_females, Red_deer_BFNP_activity_Sep_females_df_females,  Red_deer_BFNP_activity_Okt_females_df_females, Red_deer_BFNP_activity_Nov_females_df_females, Red_deer_BFNP_activity_Dez_females_df_females, Red_deer_BFNP_activity_Jan_females_df_females, Red_deer_BFNP_activity_Feb_females_df_females, Red_deer_BFNP_activity_Mrz_females_df_females, Red_deer_BFNP_activity_Apr_females_df_females)






(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mai" & is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Mai_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mai"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Mai_males_df_males<-data.frame(month="Mai", activity=Red_deer_BFNP_activity_Mai_males@act[1], se=Red_deer_BFNP_activity_Mai_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Mai_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Mai_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jun"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Jun_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jun"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Jun_males_df_males<-data.frame(month="Jun", activity=Red_deer_BFNP_activity_Jun_males@act[1], se=Red_deer_BFNP_activity_Jun_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jun_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jun_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jul"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Jul_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jul"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Jul_males_df_males<-data.frame(month="Jul", activity=Red_deer_BFNP_activity_Jul_males@act[1], se=Red_deer_BFNP_activity_Jul_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jul_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jul_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Aug"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Aug_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Aug"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Aug_males_df_males<-data.frame(month="Aug", activity=Red_deer_BFNP_activity_Aug_males@act[1], se=Red_deer_BFNP_activity_Aug_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Aug_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Aug_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Sep"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Sep_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Sep"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Sep_males_df_males<-data.frame(month="Sep", activity=Red_deer_BFNP_activity_Sep_males@act[1], se=Red_deer_BFNP_activity_Sep_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Sep_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Sep_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Okt"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Okt_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Okt"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Okt_males_df_males<-data.frame(month="Okt", activity=Red_deer_BFNP_activity_Okt_males@act[1], se=Red_deer_BFNP_activity_Okt_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Okt_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Okt_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Nov"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Nov_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Nov"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Nov_males_df_males<-data.frame(month="Nov", activity=Red_deer_BFNP_activity_Nov_males@act[1], se=Red_deer_BFNP_activity_Nov_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Nov_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Nov_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Dez"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Dez_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Dez"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Dez_males_df_males<-data.frame(month="Dez", activity=Red_deer_BFNP_activity_Dez_males@act[1], se=Red_deer_BFNP_activity_Dez_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Dez_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Dez_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jan"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Jan_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Jan"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Jan_males_df_males<-data.frame(month="Jan", activity=Red_deer_BFNP_activity_Jan_males@act[1], se=Red_deer_BFNP_activity_Jan_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Jan_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Jan_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Feb"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Feb_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Feb"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Feb_males_df_males<-data.frame(month="Feb", activity=Red_deer_BFNP_activity_Feb_males@act[1], se=Red_deer_BFNP_activity_Feb_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Feb_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Feb_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mrz"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Mrz_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Mrz"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="data", show=TRUE)
Red_deer_BFNP_activity_Mrz_males_df_males<-data.frame(month="Mrz", activity=Red_deer_BFNP_activity_Mrz_males@act[1], se=Red_deer_BFNP_activity_Mrz_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Mrz_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Mrz_males@act[4])

length(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Apr"& is.na(BFNP_red_deer$Number_males)==F])
Red_deer_BFNP_activity_Apr_males<-fitact(BFNP_red_deer$radian_time[BFNP_red_deer$month=="Apr"& is.na(BFNP_red_deer$Number_males)==F], reps=10000, sample="model", show=TRUE)
Red_deer_BFNP_activity_Apr_males_df_males<-data.frame(month="Apr", activity=Red_deer_BFNP_activity_Apr_males@act[1], se=Red_deer_BFNP_activity_Apr_males@act[2], Lower_95_CI=Red_deer_BFNP_activity_Apr_males@act[3] , Upper_95_CI=Red_deer_BFNP_activity_Apr_males@act[4])
plot(Red_deer_BFNP_activity_Apr_males)

Red_deer_BFNP_activity_df_males<-rbind(Red_deer_BFNP_activity_Mai_males_df_males, Red_deer_BFNP_activity_Jun_males_df_males, Red_deer_BFNP_activity_Jul_males_df_males, Red_deer_BFNP_activity_Aug_males_df_males, Red_deer_BFNP_activity_Sep_males_df_males,  Red_deer_BFNP_activity_Okt_males_df_males, Red_deer_BFNP_activity_Nov_males_df_males, Red_deer_BFNP_activity_Dez_males_df_males, Red_deer_BFNP_activity_Jan_males_df_males, Red_deer_BFNP_activity_Feb_males_df_males, Red_deer_BFNP_activity_Mrz_males_df_males, Red_deer_BFNP_activity_Apr_males_df_males)

rm(Red_deer_BFNP_activity_Mai_males_df_males, Red_deer_BFNP_activity_Jun_males_df_males, Red_deer_BFNP_activity_Jul_males_df_males, Red_deer_BFNP_activity_Aug_males_df_males, Red_deer_BFNP_activity_Sep_males_df_males,  Red_deer_BFNP_activity_Okt_males_df_males, Red_deer_BFNP_activity_Nov_males_df_males, Red_deer_BFNP_activity_Dez_males_df_males, Red_deer_BFNP_activity_Jan_males_df_males, Red_deer_BFNP_activity_Feb_males_df_males, Red_deer_BFNP_activity_Mrz_males_df_males, Red_deer_BFNP_activity_Apr_males_df_males)

Red_deer_BFNP_activity_df$sex<-"all"
Red_deer_BFNP_activity_df_females$sex<-"f"
Red_deer_BFNP_activity_df_males$sex<-"m"

Red_deer_BFNP_activity_df_all<-rbind(Red_deer_BFNP_activity_df, Red_deer_BFNP_activity_df_females, Red_deer_BFNP_activity_df_males)
Red_deer_BFNP_activity_df_all$month<-factor(Red_deer_BFNP_activity_df_all$month, levels=c("Jan","Feb","Mrz","Apr","Mai","Jun", "Jul", "Aug","Sep", "Okt", "Nov", "Dez"),ordered = T)
Red_deer_BFNP_activity_df_all<-Red_deer_BFNP_activity_df_all%>%group_by(sex)%>%arrange(month, .by_group=TRUE)



Trapping_rates<-read.csv("D:/Users/maik-/Documents/Maik/Rothirschprojekt Bayerischer Wald/Interreg Camera trapping/Result tables/BFNP/Trapping_rates_roe_deer_red_deer_per_sex_and_month.csv")


par(mar = c(5, 5, 3, 5))
plot(Red_deer_BFNP_activity_df_all$activity[13:24], type="o", pch=19,  ylim=c(0,1), xlab="Month", ylab="Activity", xaxt="n", bty="l")
axis(side=1, at=1:12, labels=levels(Red_deer_BFNP_activity_df_all$month))
lines(Red_deer_BFNP_activity_df_all$activity[25:36], type="o", pch=15, lty="dotted")
lines(Red_deer_BFNP_activity_df_all$Lower_95_CI[13:24], type="l", pch=15, lty="solid")
lines(Red_deer_BFNP_activity_df_all$Upper_95_CI[13:24], type="l", pch=15, lty="solid")
lines(Red_deer_BFNP_activity_df_all$Lower_95_CI[25:36], type="l", pch=15, lty="dotted")
lines(Red_deer_BFNP_activity_df_all$Upper_95_CI[25:36], type="l", pch=15, lty="dotted")
par(new = TRUE) 
plot(Trapping_rates$Trapping_rate_red_deer_females, type="l", pch=15, col="red",lty="solid", axes=FALSE, xlab="",ylab="", ylim=c(0,14))
lines(Trapping_rates$Trapping_rate_red_deer_males, type="l", pch=15, col="red",lty="dotted")
axis(side = 4, at = seq(0,14,by=1))      # Add second axis
box(col = 'black')
mtext("Events/100 days", side = 4, line = 3)  
legend("topright", legend=c("Activity females", "Activity males","CT rate females", "CT rate males"), lty=c("solid", "dotted","solid", "dotted"), col=c("black","black","red","red"),pch=c(19,15,NA,NA), inset=0.01, cex=1)
title("a)", adj = 0, line = 0.5)







plot(Red_deer_BFNP_activity_df_all$activity[1:12], type="o", pch=19,  ylim=c(0,1), xlab="Month", ylab="Activity", xaxt="n", bty="l")
axis(side=1, at=1:12, labels=levels(Red_deer_BFNP_activity_df_all$month))
lines(Red_deer_BFNP_activity_df_all$Lower_95_CI[1:12], type="l", pch=19, lty="dotted")
lines(Red_deer_BFNP_activity_df_all$Upper_95_CI[1:12], type="l", pch=19, lty="dotted")

legend("topright", legend=c("Activity females", "CT rate females"), lty=c("solid","solid"), col=c("black","red"),pch=c(19,NA), inset=0.01, cex=1)
title("a)", adj = 0, line = 0.5)







plot(Red_deer_BFNP_activity_df_all$activity[13:24], type="o", pch=19,  ylim=c(0,1), xlab="Month", ylab="Activity", xaxt="n", bty="l")
axis(side=1, at=1:12, labels=levels(Red_deer_BFNP_activity_df_all$month))
lines(Red_deer_BFNP_activity_df_all$Lower_95_CI[13:24], type="l", pch=19, lty="dotted")
lines(Red_deer_BFNP_activity_df_all$Upper_95_CI[13:24], type="l", pch=19, lty="dotted")
par(new = TRUE) 
plot(Trapping_rates$Trapping_rate_red_deer_females, type="l", pch=19, col="red",lty="solid", axes=FALSE, xlab="",ylab="", ylim=c(0,14))
axis(side = 4, at = seq(0,14,by=1))      # Add second axis
box(col = 'black')
mtext("Events/100 days", side = 4, line = 3)  
legend("topright", legend=c("Activity females", "CT rate females"), lty=c("solid","solid"), col=c("black","red"),pch=c(19,NA), inset=0.01, cex=1)
title("a)", adj = 0, line = 0.5)


plot(Red_deer_BFNP_activity_df_all$activity[25:36], type="o", pch=15,  ylim=c(0,1), xlab="Month", ylab="Activity", xaxt="n", bty="l")
axis(side=1, at=1:12, labels=levels(Red_deer_BFNP_activity_df_all$month))
lines(Red_deer_BFNP_activity_df_all$Lower_95_CI[25:36], type="l", pch=15, lty="dotted")
lines(Red_deer_BFNP_activity_df_all$Upper_95_CI[25:36], type="l", pch=15, lty="dotted")
par(new = TRUE) 
plot(Trapping_rates$Trapping_rate_red_deer_males, type="l", pch=15, col="red",lty="solid", axes=FALSE, xlab="",ylab="", ylim=c(0,14))
axis(side = 4, at = seq(0,14,by=1))      # Add second axis
box(col = 'black')
mtext("Events/100 days", side = 4, line = 3)  
legend("topright", legend=c("Activity males", "CT rate males"), lty=c("solid","solid"), col=c("black","red"),pch=c(19,NA), inset=0.01, cex=1)
title("b)", adj = 0, line = 0.5)






##################################### Sensitivity analysis #############################################################################
set.seed(45701)
Dist_red_deer_sensitivity_May_2018<-Dist_conf_function_months_sensitivity(Camera_lag_3,Red_deer_times_photos_BFNP[Red_deer_times_photos_BFNP$month_year=="Mai_2018",], Det_prop, Det_prop_SE, trunc_dist=14, Prop_within_trunc, angle=55,Red_deer_BFNP_activity_df, Group_size, 2000)


Dist_red_deer_sensitivity_May_2018$Trap_rate_scaled<-as.numeric(scale(Dist_red_deer_sensitivity_May_2018$Trap_rate))
Dist_red_deer_sensitivity_May_2018$Det_prob_scaled<-as.numeric(scale(Dist_red_deer_sensitivity_May_2018$Det_prob))
Dist_red_deer_sensitivity_May_2018$Activity_scaled<-as.numeric(scale(Dist_red_deer_sensitivity_May_2018$Activity))
Dist_red_deer_sensitivity_May_2018$Delay_scaled<-as.numeric(scale(Dist_red_deer_sensitivity_May_2018$Delay))
Dist_red_deer_sensitivity_May_2018$Group_size_scaled<-as.numeric(scale(Dist_red_deer_sensitivity_May_2018$Group_size)) 
Reg_mod_red_deer_May_2019_scaled<-lm(Density ~ Trap_rate_scaled + Det_prob_scaled + Activity_scaled + Delay_scaled + Group_size_scaled + (Trap_rate_scaled + Det_prob_scaled + Activity_scaled + Delay_scaled + Group_size_scaled)^2, data=Dist_red_deer_sensitivity_May_2018)
summary(Reg_mod_red_deer_May_2019_scaled)
anova(Reg_mod_red_deer_May_2019_scaled)
Anova(Reg_mod_red_deer_May_2019_scaled,type=c("III"))
Prop_ss_red_deer_May_2019<-round(anova(Reg_mod_red_deer_May_2019_scaled)[1:15,2]*100/sum(anova(Reg_mod_red_deer_May_2019_scaled)[1:15,2]),2)












######################## Functions #################################################################################################################################

Grp_size_est<-function(Grp.dataframe ,indices){
  d<-Grp.dataframe[indices,]
  Grp_size<-mean(d$Number_of_animals_species_1, na.rm=TRUE)
  return(Grp_size)
}




Dist_conf_function_months_base=function(Camera_lag, Cam_tab_all, Det_prop, Det_prop_SE, trunc_dist, Prop_within_trunc, angle, Activity_tab_all, Group_tab_all,  nboot, alpha)
{
  id=levels(droplevels(Cam_tab_all$month_year)) 
  n=length(id)
  Distlist=list()
  
  for(i in 1:n){
    Cam_tab=Cam_tab_all[Cam_tab_all$month_year==id[i],]
    Activity_tab=Activity_tab_all[Activity_tab_all$month==unique(Cam_tab$month),]
    Group_tab=Group_tab_all[Group_tab_all$Season==unique(Cam_tab$Season),]
    whi=trunc_dist/1000
    
    n_cam=length(Cam_tab$Cam_secs)
    
    Snap<-median(Camera_lag)
    Nr_photos_obs=sum(Cam_tab$Nr_photo_series)                                                 
    Secs_total_obs=sum(Cam_tab$Cam_secs)
    theta<-(angle*pi)/180
    Activity<-Activity_tab$activity
    g_obs<-unique(Group_tab$Group_size_t0) 
    
    Dens_obs<-(2*Snap * Nr_photos_obs*Prop_within_trunc)/(theta*whi^2*Secs_total_obs*Det_prop)*(1/Activity)*g_obs

    Dens_est<-numeric(nboot)
    
    for (nsim in 1: nboot) {
      Snap_boot<-Camera_lag[nsim]
      obs_boot = sample_n(Cam_tab, n_cam, replace=TRUE) 
      Nr_photos_boot<-sum(obs_boot$Nr_photo_series)                                                                                   
      Secs_total_boot<-sum(obs_boot$Cam_secs)
      Det_prop_boot<-rnorm(1,Det_prop, 1.96 *Det_prop_SE)
      Activity_boot<-rnorm(1,Activity_tab$activity, 1.96 *Activity_tab$se)
      g_boot<-Group_tab$Group_size[nsim]                
      
      Dens_est[nsim]<-(2*Snap_boot * Nr_photos_boot *Prop_within_trunc)/(theta*whi^2*Secs_total_boot*Det_prop_boot)*(1/Activity_boot)*g_boot
      
    }
    b               = qnorm((sum(Dens_est>Dens_obs)+sum(Dens_est==Dens_obs)/2)/length(Dens_est))
    z               = qnorm(c(alpha/2,1-alpha/2))      
    p               = pnorm(z-2*b)                      
    Conf_biasCorrected = quantile(Dens_est ,p=p)  
    Dist_est<-data.frame(month=id[i],Est=Dens_obs, Lower95conf=Conf_biasCorrected[1], Upper95conf=Conf_biasCorrected[2],Photos=sum(Cam_tab$Nr_photo_series), Deployment_secs=sum(Cam_tab$Cam_secs))
    Distlist[[i]]<-Dist_est
  }
  Dist_dat<-Distlist[[1]]
  for (i in 2:n){
    Dist_dat=rbind(Dist_dat,Distlist[[i]])}
  return(Dist_dat)
}





Dist_conf_function_months_plotting=function(Camera_lag, Cam_tab_all, Det_prop, Det_prop_SE, trunc_dist, Prop_within_trunc, angle, Activity_tab_all, Group_tab_all,  nboot)
{
  id=levels(droplevels(Cam_tab_all$month_year)) 
  n=length(id)
  Distlist=list()
  
  for(i in 1:n){
    Cam_tab=Cam_tab_all[ Cam_tab_all$month_year==id[i],]
    Activity_tab=Activity_tab_all[Activity_tab_all$month==unique(Cam_tab$month),]
    Group_tab=Group_tab_all[Group_tab_all$Season==unique(Cam_tab$Season),]
    whi=trunc_dist/1000
    
    n_cam=length(Cam_tab$Cam_secs)
    
    Snap<-median(Camera_lag)
    Nr_photos_obs=sum(Cam_tab$Nr_photo_series)                                                 
    Secs_total_obs=sum(Cam_tab$Cam_secs)
    theta<-(angle*pi)/180
    Activity<-Activity_tab$activity
    g_obs<-unique(Group_tab$Group_size_t0) 
    
    Dens_obs<-(2*Snap * Nr_photos_obs*Prop_within_trunc)/(theta*whi^2*Secs_total_obs*Det_prop)*(1/Activity)*g_obs
    
    Dens_est<-numeric(nboot)
    
    for (nsim in 1: nboot) {
      Snap_boot<-Camera_lag[nsim]
      obs_boot = sample_n(Cam_tab, n_cam, replace=TRUE) 
      Nr_photos_boot<-sum(obs_boot$Nr_photo_series)                                                                                   
      Secs_total_boot<-sum(obs_boot$Cam_secs)
      Det_prop_boot<-rnorm(1,Det_prop, 1.96 *Det_prop_SE)
      Activity_boot<-rnorm(1,Activity_tab$activity, 1.96 *Activity_tab$se)
      g_boot<-Group_tab$Group_size[nsim]                
      
      Dens_est[nsim]<-(2*Snap_boot * Nr_photos_boot *Prop_within_trunc)/(theta*whi^2*Secs_total_boot*Det_prop_boot)*(1/Activity_boot)*g_boot
      
    }
    Dist_data<-data.frame(month=id[i], Range=Dens_est)
    Distlist[[i]]<-Dist_data
  }
  Dist_dat<-Distlist[[1]]
  for (i in 2:n){
    Dist_dat=rbind(Dist_dat,Distlist[[i]])}
  return(Dist_dat)
}




Dist_conf_function_months_sensitivity=function(Camera_lag, Cam_tab, Det_prop, Det_prop_SE, trunc_dist, Prop_within_trunc, angle, Activity_tab_all, Group_tab_all,  nboot)
  
{
  Sim_list=list()
  Activity_tab=Activity_tab_all[Activity_tab_all$month==unique(Cam_tab$month),]
  Group_tab=Group_tab_all[Group_tab_all$Season==unique(Cam_tab$Season),]
  whi=trunc_dist/1000
  theta<-(angle*pi)/180
  n_cam=length(Cam_tab$Cam_secs)
  
  for (nsim in 1: nboot) {
    Snap_boot<-Camera_lag[nsim]
    obs_boot = sample_n(Cam_tab, n_cam, replace=TRUE) 
    Nr_photos_boot<-sum(obs_boot$Nr_photo_series)                                                                                   
    Secs_total_boot<-sum(obs_boot$Cam_secs)
    Det_prop_boot<-rnorm(1,Det_prop, 1.96 *Det_prop_SE)
    Activity_boot<-rnorm(1,Activity_tab$activity, 1.96 *Activity_tab$se)
    g_boot<-Group_tab$Group_size[nsim]                
    
    Dens_est<-(2*Snap_boot * Nr_photos_boot *Prop_within_trunc)/(theta*whi^2*Secs_total_boot*Det_prop_boot)*(1/Activity_boot)*g_boot
    
Sim_list[[nsim]]<-data.frame(month=unique(obs_boot$month_year), Trap_rate=(Nr_photos_boot/Secs_total_boot)*86400, Det_prob=Det_prop_boot, Activity= Activity_boot, Group_size=g_boot, Delay=Snap_boot ,Density=Dens_est)
  }
  
  Sim_dat<-bind_rows(Sim_list)
  return(Sim_dat)
}

















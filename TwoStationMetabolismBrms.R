
## TWO-STATION METABOLISM ESTIMATION WITH BRMS
#############################################################################################################################
#############################################################################################################################

#NOTES:
###########################################################################################################################
  #This R script is greatly influenced by the one reported in Hall et al.(2016,Ecosystems)
  #Before fitting the models it is important to check quality data. I usually do that with "xts" and "dygraphs" packages
  #Sometimes it is necessary to remove outliers, and can be imputated using na.approx() or na.spline()
  #In my experience, light can strongly affect fitting output. I think it is worth doing some tests with different lights 
  #I usually compare light from sensors with open-sky modeled light calculated with StreamMetabolizer::calc_light()
  #I will import a dataset with DOsat already in it, but it is also possible to import barometric pressure (bp) and calculate DOsat in R
  #Check Chapter34- Stream Metabolism- Hall & Hotchkiss (2017)in Methods in Stream Ecology:
      #osat<- function(temp, bp) {
      
      #tstd<-log((298.15-temp) / (273.15 + temp))
      #a0<-2.00907
      #a1<-3.22014
      #a2<-4.0501
      #a3<-4.94457
      #a4<- -0.256847
      #a5<- 3.88767
      
      #u<-10^(8.10765-(1750.286/(235+temp)))
      
      #sato<-(exp(a0 + a1*tstd + a2*tstd^2 + a3*tstd^3 + a4*tstd^4+a5*tstd^5))*((bp-u)/(760-u))*1.42905
      #sato
      #}

#Packages and functions
###########################################################
library(tidyverse)  #includes dplyr, ggplot2, readr
library(tibbletime)
library(clipr)
library(magrittr)
library(ggpubr)

#library(xts)
#library(dygraphs)

library(brms)


#The following function is to get to round 2.5 to 3, not 2 (and similar cases), as the normal round() function would do
  round2 = function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5 + sqrt(.Machine$double.eps)
    z = trunc(z)
    z = z/10^n
    z*posneg
  }


##1. Import data & subset observations
#############################################################################################################################
#############################################################################################################################

#1a) Import
#############################################################

  dat_all<-read_clip_tbl()
  dat_all$DateTime<-as.POSIXct(dat_all$DateTime)
  dat_all<- as_tbl_time(dat_all, index = DateTime)
  
  dat_all[,c("Stream","Reach","Period")] <- lapply(dat_all[,c("Stream","Reach","Period")],as.factor)


#1b) Subset observations
############################################################

  dat<-dat_all%>% 
          filter( Stream=="Aarhus",
                  Reach=="Control",
                  Period=="Before")
  
  #or filter by time:
      #dat<-dat%>%
      #       filter_time('2018-06-14' ~ '2018-06-16')
  
  #example of initial dataset:
  
    # A time tibble: 6 x 13
    # Index: DateTime
    # DateTime            Stream Reach   Period temp_up DOsat_up DO_up temp_dw DOsat_dw DO_dw light     z    tt
    # <dttm>              <fct>  <fct>   <fct>    <dbl>    <dbl> <dbl>   <dbl>    <dbl> <dbl> <dbl> <dbl> <int>
    # 1 2018-06-14 14:00:00 Aarhus Control Before    18.9     9.21  11.0    18.9     9.21  10.7 1160. 0.074    15
    # 2 2018-06-14 14:05:00 Aarhus Control Before    18.9     9.20  10.9    18.9     9.21  10.6 1026. 0.074    15
    # 3 2018-06-14 14:10:00 Aarhus Control Before    18.9     9.20  10.9    18.8     9.22  10.6 1323. 0.074    15
    # 4 2018-06-14 14:15:00 Aarhus Control Before    18.9     9.20  10.9    18.8     9.22  10.6  776. 0.074    15
    # 5 2018-06-14 14:20:00 Aarhus Control Before    19.0     9.19  10.8    18.9     9.21  10.6  979  0.074    15
    # 6 2018-06-14 14:25:00 Aarhus Control Before    19.0     9.18  10.9    18.9     9.21  10.5  955. 0.074    15
  
  summary(dat)


##2. Data wrangling: align upstream and downstream data, and convertion factors
###################################################################################################################
####################################################################################################################

  meas_interval<-5 #minutes between DO measures
  dat$tt_day<-dat$tt/(60*24)  #travel time from minutes to days
  
  lag_tt<-round2(dat$tt[1]/meas_interval,0) #lag between up and down signals in number of rows
  
  
  
  #DateTime can be changed as the metabolism estimates refer to what occurred to the water between when it was upstream at time t
  #and when it was downstream at time t+lag_tt
  
  dat%<>%mutate(DateTime=DateTime+(dat$tt[1]*60/2),  
                   temp_up=lag(temp_up,n=lag_tt),
                   DOsat_up=lag(DOsat_up,n=lag_tt),
                   DO_up=lag(DO_up,n=lag_tt))
  
  
  #The light that is producing the photosynthesis that will increase DO downstream should be integrated between t and t+lag_tt;
  #However, after several tests, I usually get better models when light is just integrated for half of lag time
  #I interpret this as a greater part of downstream DO signal is produced closer to dw station and time
  #I apply the same formula for temperature:
  
  int_t<-round2((lag_tt+1)/2,0)
  rolling_mean_int_t <- rollify(mean, window = int_t)
  
  #
  dat%<>%mutate(light=rolling_mean_int_t(light),
                   temp=rolling_mean_int_t(temp_dw))     # or temp=rowMeans(select(., c("temp_up","temp_dw")))

  #Optional: plot DO signals and light (packages xts and dygraphs)
  ###############################################################################
    #dat.xts<-xts(dat[-1],order.by=dat$DateTime)
  
    #dygraph(dat.xts[,c("DO_up","DO_dw","light")])%>% 
      #dyAxis("y2",drawGrid = F) %>%
      #dySeries("light", axis = 'y2')%>%dyRangeSelector()


#Reaeration flux convertion factors
#################################################################################

  dat$O2def<-with(dat, (DOsat_up-DO_up+DOsat_dw)/2)
  
  dat$conv<- with(dat,((1568-(temp*86.04)+(2.142*temp^2)-(0.0216*temp^3))/600)^-0.5)


##3. Model fitting
###################################################################################################################
####################################################################################################################

#Subset definitive range period and predictors
########################################################################################
  
  datdef<-dat%>%filter_time('2018-06-14 15:00' ~ '2018-06-16 03:00')
  #datdef<-datdef%>%filter(light==0) #sometimes is useful to compare your K600 and ER20 estimates with a night model,
                                     # where the photosynthetic part (Pmax, alpha) is removed
  
  datdef<-mutate(datdef, date = collapse_index(DateTime, "daily"))
  datdef$date<-as.Date(datdef$date)
  datdef$date<-as.factor(datdef$date)
  
  datdef$time<-seq(0, 1, length.out = length(datdef$DateTime)) 
  
#Priors
##########################################################################################
  
  priors<-c( prior(normal(20,20),nlpar="K600",lb=0)+   
             prior(normal(7,7),nlpar="Pmax",lb=0)+
             prior(normal(20,20),nlpar="alpha",lb=0)+
             prior(normal(7,7),nlpar="ER20",lb=0)
  )
  
#Formula
##########################################################################################
  

  formula<- bf(DO_dw ~ (DO_up+
                         ((tt_day/z)*Pmax*tanh(alpha*light/(1000*Pmax)))+  #1000 factor, from umols to mmols;  Pmax*tanh(alpha*light/Pmax))
                         (O2def*tt_day*K600*conv)-
                         ((tt_day/z)*ER20)*1.045^(temp-20))/            #*1.045^(temp-20)
                (1+(K600*conv*tt_day/2)),  
              K600~1,  #+s(time,k=24)  +data
              ER20~1,
              alpha~1,
              Pmax~1,
              nl=TRUE)
  
  
#Fitting
#################################################################################################################
  
  fit1<-brm(formula=formula,
            data = datdef,
            family=gaussian(),
            prior = priors ,
            control = list(adapt_delta = 0.99,max_treedepth=10),iter=2000,cores=4)


##4. Model evaluation
###################################################################################################################
####################################################################################################################

  print(summary(fit1),digits=6)
  plot(fit1,ask=FALSE)
  
  pp_check(fit1)
  loo(fit1)
  
  pp_fit1<-cbind(datdef,
                  fitted(fit1),
                  fitted(fit1,nlpar="K600"),
                  fitted(fit1,nlpar="alpha"),
                  fitted(fit1,nlpar="Pmax"),
                  fitted(fit1,nlpar="ER20"))
  
  colnames(pp_fit1)[(ncol(datdef)+1):(ncol(pp_fit1))]<-c(
                       "oxy_dw_Estimate","oxy_dw_Est.Error","oxy_dw_Q2.5","oxy_dw_Q97.5",
                       "K600_Estimate","K600_Est.Error","K600_Q2.5","K600_Q97.5",
                       "alpha_Estimate","alpha_Est.Error","alpha_Q2.5","alpha_Q97.5",
                       "Pmax_Estimate","Pmax_Est.Error","Pmax_Q2.5","Pmax_Q97.5",
                       "ER20_Estimate","ER20_Est.Error","ER20_Q2.5","ER20_Q97.5")
  
  pp_fit1<- pp_fit1 %>% as_tibble %>%as_tbl_time(DateTime)
  
  
#GPP and ER instantaneous
#############################################################################################################################
  
  pp_fit1%<>%
    mutate(
      GPP=Pmax_Estimate*tanh(alpha_Estimate*datdef$light/(Pmax_Estimate*1000)),
      GPP_low=Pmax_Q2.5*tanh(alpha_Q2.5*datdef$light/(Pmax_Estimate*1000)),
      GPP_high=Pmax_Q97.5*tanh(alpha_Q97.5*datdef$light/(Pmax_Estimate*1000)),
      ER= ER20_Estimate*1.045^(datdef$temp-20),
      ER_low= ER20_Q2.5*1.045^(datdef$temp-20),
      ER_high= ER20_Q97.5*1.045^(datdef$temp-20)
    )
  
#Errors ans Summary
###########################################################################################################################
  
  pp_fit1%<>% mutate(AE=abs(DO_dw-oxy_dw_Estimate))
  
  pp_fit1 %>%
    #group_by(date)%>%
    filter_time('2018-06-15' ~ '2018-06-15')%>%
    select(GPP,GPP_low,GPP_high,ER,ER_low,ER_high,AE)%>%
    summarise_all(mean)
  
  #Rocher-Ros et al. (2019,Global Change Biol.) set a threshold of 0.2 for MAE
  
  
#Export dataset
###################################################################################
  write.csv(pp_fit1,file="fit1.csv")
#################################################################################
  
  
#Viusal check: estimated and observed DO
#############################################################################################################################
  
  pp_fit1%>%
    select(DateTime,DO_up,DO_dw,oxy_dw_Estimate)%>%
    gather(key=DO_type,value=DO,c("DO_up","DO_dw","oxy_dw_Estimate"))%>%
    ggplot(aes(x = DateTime,y=DO,group=DO_type)) + 
    geom_vline(xintercept=as.POSIXct("2018-06-25 11:00"),alpha=0.6,size=1.5)+
    geom_point(aes(color=DO_type,size=DO_type,alpha=DO_type))+
    geom_line(aes(linetype=DO_type,color=DO_type),size=1.3)+
    
    scale_alpha_manual(values=c(0.25,0.25,0),labels=c("DOdw","DOup","DOdw estimated"))+
    scale_size_manual(values=c(3.5,3.5,2.5),labels=c("DOdw","DOup","DOdw estimated"))+
    scale_color_manual(values=c("#f16913","#cccccc","#bd0026"),labels=c("DOdw","DOup","DOdw estimated"))+
    scale_linetype_manual(values=c("blank","blank","solid"),labels=c("DOdw","DOup","DOdw estimated"))+
    scale_x_datetime(date_breaks = "1 day",date_labels = "%d/%m",expand = c(0.08,0.08))+
    xlab(bquote(''))+
    ylab(bquote('DO (mg'*~~O[2]~~ L^-1*')'))+
    theme_classic()+
    theme(
      legend.position = c(0.9, 0.93),
      legend.title=element_blank())
  
  #ggsave("fit1_DO.tiff",width=12,height=6,dpi=300)
  
  
#Plot instantaneous metabolism
#############################################################################################################################
  
  pp_fit1%>%
    select(DateTime,GPP,GPP_low,GPP_high,ER,ER_low,ER_high)%>%
    ggplot(aes(x=DateTime)) + 
    geom_ribbon(aes(ymin=GPP_low,ymax=GPP_high), fill="#31a354",alpha = 0.8)+
    geom_line(aes(y = GPP), color="black") +
    geom_ribbon(aes(ymin=abs(ER_low),ymax=abs(ER_high)),  fill="#980043",alpha = 0.8)+
    geom_line(aes(y = abs(ER)), color="black") +
    scale_y_continuous(expand=c(0,0),breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))+
    scale_x_datetime(date_breaks = "1 day",date_labels = "%d/%m",expand = c(0.09,0.09))+
    ylab(bquote('Metabolism  (g'*~~O[2]~~ m^-2~~day^-1*')'))+
    xlab(bquote(''))+
    theme_classic()
  
  
  #ggsave("fit1_metinst.tiff",width=12,height=6,dpi=300)
  
  
#Estimated parameters  (K600,ER20,alpha,Pmax)
#############################################################################################################################
  
  pp_fit1%>%
    ggplot(aes(x=DateTime))+
    geom_ribbon(aes(ymin = K600_Q2.5, ymax = K600_Q97.5), fill = "grey70") +
    geom_line(aes(y = K600_Estimate),size=1.3)+
    scale_x_datetime(date_breaks = "2 days",date_labels = "%d/%m",expand = c(0.1,0.1))+
    xlab(bquote(''))+
    ylab(bquote('K600 ('*d^-1*')'))+
    theme_classic()
  
  #ggsave("fit1_K600.tiff",width=4,height=6,dpi=300)
  
  
  
  pp_fit1%>%
    ggplot(aes(x=DateTime))+
    geom_ribbon(aes(ymin = alpha_Q2.5, ymax = alpha_Q97.5), fill = "grey70") +
    geom_line(aes(y = alpha_Estimate),size=1.3)+
    scale_x_datetime(date_breaks = "2 days",date_labels = "%d/%m",expand = c(0.1,0.1))+
    xlab(bquote(''))+
    ylab(bquote(alpha*'  (g'*~O[2]~~s~~mmolPAR^-1~~day^-1*')'))+
    theme_classic()
  
  
  #ggsave("fit1_alpha.tiff",width=12,height=6,dpi=300)
  
  
  pp_fit1%>%
    ggplot(aes(x=DateTime))+
    geom_ribbon(aes(ymin = Pmax_Q2.5, ymax = Pmax_Q97.5), fill = "grey70") +
    geom_line(aes(y = Pmax_Estimate),size=1.3)+
    scale_x_datetime(date_breaks = "2 days",date_labels = "%d/%m",expand = c(0.1,0.1))+
    xlab(bquote(''))+
    ylab(bquote('Pmax (g'*~O[2]~~m^-2~~day^-1*')'))+
    theme_classic()
  
  #ggsave("fit1_Pmax.tiff",width=12,height=6,dpi=300)
  
  
  pp_fit1%>%
    ggplot(aes(x=DateTime))+
    geom_ribbon(aes(ymin = ER20_Q2.5, ymax = ER20_Q97.5), fill = "grey70") +
    geom_line(aes(y = ER20_Estimate),size=1.3)+
    scale_x_datetime(date_breaks = "2 days",date_labels = "%d/%m",expand = c(0.1,0.1))+
    xlab(bquote(''))+
    ylab(bquote('ER20 (g'*~O[2]~~m^-2~~day^-1*')'))+
    theme_classic()
  
  #ggsave("fit1_ER20.tiff",width=12,height=6,dpi=300)
  
  
#############################################################################################################################
#############################################################################################################################
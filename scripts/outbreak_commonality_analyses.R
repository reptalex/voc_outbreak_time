library(tidyverse)
library(COVID19)
library(padr)
library(lubridate)
library(progress)
library(RcppRoll)
library(parallel)
library(data.table)
library(ggpubr)
library(usmap)
library(tsoutliers)
library(mgcv)
source("scripts/delta_omicron_utils.R")



# Settings & params -------------------------------------------------------

# county_population_threshold <- 1e6 # Can be useful when conducting county-level analyses whose 
                                    ## oultier detection + growth rate estimation are time-consuming
theme_set(theme_bw(base_size=12))
ncores=6


# Initialize cluster for parallelized outlier detection --------------------------------------------

cl <- makeCluster(ncores)
clusterEvalQ(cl,expr = {library(tsoutliers)
  library(magrittr)
  library(tidyverse)})
clusterExport(cl,'outlier_detection')


# Countries ---------------------------------------------------------------

World <- covid19() %>% as.data.table

countries <- c('United Kingdom','South Africa','India')

World[,country:=administrative_area_level_1]
World <- World[country %in% countries]
World[,new_confirmed:=dfs(confirmed),by=country]
World[,new_deaths:=dfs(deaths),by=country]
World <- World[date>=as.Date('2020-10-01')]

######### Outlier removal
World[,new_confirmed:=outlier_detection_par(new_confirmed,country,cl)]
World <- covid19_nbss(World,mc.cores = ncores) %>% as.data.table
## Note of acknowledgement: Justin Silverman wrote the covid19_nbss function for our earlier work
## analyzing COVID outbreaks on a timescale of burden. He's an exceptionally talented statistician
## and his functions have proven useful throughout the pandemic, allowing growth rate estimation
## to be standardized, with filtering/smoothign made easy. Give him all your grants! :-D

# US states ---------------------------------------------------------------

US_states <- covid19('United States',level=2) %>% as.data.table
US_states[,country:=administrative_area_level_1]
US_states[,state:=administrative_area_level_2]
US_states <- US_states[!state %in% c('Guam','Northern Mariana Islands','Virgin Islands','American Samoa')]
US_states[,new_confirmed:=dfs(confirmed),by=state]
US_states[,new_deaths:=dfs(deaths),by=state]
US_states <- US_states[date>as.Date('2021-04-01')]

### some states changed their reporting practices after the B 1.1.7 wave to not report at all on e.g. weekends.
### Below isn't a perfect fix, but hopefully helps the algorithm converge to reliable growth rate estimates for recent weeks.
### We make an additional fix later by separately estimating growth rates for Delta & Omicron waves.
US_states[date>as.Date('2021-06-01') & new_confirmed==0,new_confirmed:=NA]
US_states[date>as.Date('2021-06-01') & new_deaths==0,new_deaths:=NA]

US_states[state=='Nevada' & date==as.Date('2021-10-06'),new_confirmed:=NA] 
## Nevada outlier not picked up well by outlier_detection
thanksgiving <- c(seq(as.Date('2020-11-26'),as.Date('2020-11-27'),by='day'),
                  seq(as.Date('2021-11-25'),as.Date('2021-11-29'),by='day'))
christmas <- as.Date(c('2020-12-25','2020-12-26'))
ny<- c(seq(as.Date('2021-01-01'),as.Date('2021-01-02'),by='day'),
       seq(as.Date('2022-01-01'),as.Date('2022-01-02'),by='day'))
lbr <- as.Date(c('2021-09-06','2021-09-07'))
ipd <- as.Date(c('2021-10-11','2021-10-12'))
hlwn <- as.Date(c('2021-10-30','2021-10-31','2021-11-01'))
ind_day <- seq(as.Date('2021-07-03'),as.Date('2021-07-06'),by='day')
US_states[state=='Missouri' & new_confirmed>4e4 & date<as.Date('2021-05-01'),new_confirmed:=NA] ## Missouri outlier

#### Outlier removal
US_states[,new_confirmed:=outlier_detection_par(new_confirmed,state,cl)]


### fourth of July holiday - monkey wrench of data delays+dumps in critical period of Delta growth
US_states[state=='California' & date %in% as.Date(c('2021-07-01','2021-07-08')),new_confirmed:=NA]


US_states[date %in% c(thanksgiving,christmas,ny,lbr,ipd,hlwn,ind_day),new_confirmed:=NA]
US_states[date %in% c(thanksgiving,christmas,ny,lbr,ipd,hlwn,ind_day),new_deaths:=NA]
US_states <- covid19_nbss(US_states,mc.cores = ncores) %>% as.data.table

US_states[state=='North Dakota' & date==as.Date('2021-08-31'),hosp:=NA]
US_states[,hosp_gr:=nbs(hosp,dispersion=1),by=state]

# Populous US counties ----------------------------------------------------

### This takes a little bit longer. 

# US_counties <- covid19("United States",level=3) %>% as.data.table
# US_counties[,country:=administrative_area_level_1]
# US_counties[,state:=administrative_area_level_2]
# US_counties[,county:=administrative_area_level_3]
# US_counties <- US_counties[!state %in% c('Guam','Northern Mariana Islands','Virgin Islands','American Samoa')]
# US_counties[,new_confirmed:=dfs(confirmed),by=c('state','county')]
# US_counties[,new_deaths:=dfs(deaths),by=c('state','county')]
# US_counties[,lbl:=paste(state,county,sep=',')]
# US_counties <- US_counties[population>=county_population_threshold]
# 
# US_counties[,new_confirmed:=outlier_detection_par(new_confirmed,lbl,cl)]
# 
# US_counties <- covid19_nbss(US_counties,mc.cores = ncores)


# South African Provinces -------------------------------------------------

ZA <- covid19('South Africa',level=2) %>% as.data.table
ZA[,province:=administrative_area_level_2]
ZA[,new_confirmed:=dfs(confirmed),by=province]
ZA[,new_deaths:=dfs(deaths),by=province]

##### Oytlier removal
ZA[,new_confirmed:=outlier_detection_par(new_confirmed,province,cl)]

ZA[date==as.Date('2021-11-23'),new_confirmed:=NA]
ZA[date==as.Date('2021-11-23'),new_deaths:=NA]

ZA <- covid19_nbss(ZA,mc.cores = ncores) %>% as.data.table

stopCluster(cl)
rm('cl')
gc()

save(list=ls(),file='data/outbreak_duration_workspace.Rd')

# Delta wave classification -----------------------------------------------------

### name "surge/decline/peak/valley" - used here to ID valleys for outbreak start-dates.
trend_namer <- function(r){
  y <- sign(r)
  dy <- shift(y)
  xx <- 1:length(r)
  
  valleys <- y>0 & dy<0
  peaks <- y<0 & dy>0
  
  valleys[is.na(valleys)] <- FALSE
  peaks[is.na(peaks)] <- FALSE
  declines <- which(cumsum(peaks)-cumsum(valleys)==1)
  waves <- which(cumsum(peaks)-cumsum(valleys)==0)
  
  wave_category <- rep('decline',length(r))
  wave_category[waves] <- 'surge'
  wave_category[peaks] <- 'peak'
  wave_category[valleys] <- 'valley'
  return(wave_category)
}

### Michigan's case reporting is very abnormal. This is what we've used for Michigan, even though Michigan
### doesn't play a prominent role in our current analyses + correspondence.
michigan_namer <- function(x){
  nas <- min(which(!is.na(x)))-1
  
  xx=1:length(x)
  
  fit <- mgcv::gam(x~s(xx),family=nb)
  r <- rep(NA,length(x))
  r[!is.na(x)] <- c(NA,diff(fit$fitted.values))
  
  y <- sign(r)
  dy <- shift(y)
  xx <- 1:length(r)
  
  valleys <- y>0 & dy<0
  peaks <- y<0 & dy>0
  
  valleys[is.na(valleys)] <- FALSE
  peaks[is.na(peaks)] <- FALSE
  declines <- which(cumsum(peaks)-cumsum(valleys)==1)
  waves <- which(cumsum(peaks)-cumsum(valleys)==0)
  
  wave_category <- rep('decline',length(r))
  wave_category[waves] <- 'surge'
  wave_category[peaks] <- 'peak'
  wave_category[valleys] <- 'valley'
  return(wave_category)
}

### Delta wave finder - valley with the lowest mean_position within the date range [min_date,max_date]
delta_wave_finder <- function(date,mean_position,wave,min_date=as.Date('2021-04-01'),max_date=Inf){
  valleys <- which(wave=='valley' & date>=min_date & date<=max_date)
  if (length(valleys)==1){
    return(date>=date[valleys])
  } else {
    lower_valley <- valleys[which.min(mean_position[valleys])]
    return(date>=date[lower_valley])
  }
}

World[,wave:=trend_namer(growth_rate),by=country]
US_states[,wave:=trend_namer(growth_rate),by=state]
US_states[state=='Michigan',wave:=michigan_namer(new_confirmed)]


World[,delta_wave:=delta_wave_finder(date,mean_position,wave,max_date = as.Date('2021-09-01')),by=country]

World[country %in% c('India','South Africa'),
      delta_wave:=delta_wave_finder(date,mean_position,wave,min_date=as.Date('2021-02-01'),max_date=as.Date('2021-11-01')),
      by=country]
World[country=='South Africa' & date<as.Date('2021-04-24'),delta_wave:=FALSE]

##### A note on Delta wave classification in South Africa #####
### Visual inspection of South Africa's growth rates, estimated with filtering=TRUE, show a period of 
### rising growth rates starting early April followed by decelerations ending near Freedom day prior to the bulk of their Delta outbreak.
### This early wobble in growth rates produces uncertainty in the start-date of the SA Delta outbreak.
### In South Africa, the lowest valley of cases occured on 4/9/21, accelerated, then decelerated, then almost reached another saddle
### in late Paril. By 5-01-2021, cases growth rates were sustained & accelerating up to the peak of cases in the SA Delta wave.
### This wobbling growth rate is clearly impacted by Freedom Day reporting blips around 4-27-2021,
### yet there was also a clear deceleration of growth rate prior to Freedom day starting 4-14-21.

### This unfortunate timing of the South African Delta outbreak's start-date adds uncertainty about
### the exact timing of the SA Delta outbreak. Our cutoff of 4-24 was inspired by the hockey-stick start of Delta-like growth in cases
### visible in the scatter plot below.

##### Plot of South African cases & growth rates, highlighting uncertainty of Delta outbreak start-date
ggarrange(
  
  ggplot(World[country=='South Africa' & date<as.Date('2021-08-01') & date>as.Date('2021-03-01')],
         aes(date,new_confirmed,color=delta_wave,fill=delta_wave))+
    geom_point(cex=4)+scale_y_continuous(trans='log')+
    geom_vline(xintercept = as.Date(c('2021-04-09','2021-4-27')))+
    geom_vline(xintercept = as.Date(c('2021-04-24')),lty=2),
  
  ggplot(World[country=='South Africa' & date<as.Date('2021-08-01') & date>as.Date('2021-03-01')],aes(date,growth_rate,color=delta_wave))+
    geom_line(lwd=2)+
    geom_vline(xintercept = as.Date(c('2021-04-09','2021-4-27')))+
    geom_vline(xintercept = as.Date(c('2021-04-24')),lty=2)+
    geom_hline(yintercept = 0),
  
nrow=2,align='v')

### These dates are chosen to cut off the US Alpha wave from impacting our estimation.
US_states[,delta_wave:=delta_wave_finder(date,mean_position,wave,
                                         min_date=as.Date('2021-04-18'),
                                         max_date=as.Date('2021-09-01')),by=state]

### Visual sanity-checking of US Delta outbreak classifications
ggplot(US_states[date<as.Date('2021-11-23')],aes(date,new_confirmed,color=delta_wave,fill=delta_wave))+
  geom_bar(stat='identity')+
  facet_wrap(.~state,scales = 'free_y')

### the same, but for early Delta outbreaks in India, UK and South Africa.
ggplot(World[date<as.Date('2021-11-01')],aes(date,new_confirmed,color=delta_wave,fill=delta_wave))+
  geom_bar(stat='identity')+
  facet_wrap(.~country,scales = 'free_y')


### Delta outbreak time, t_Delta
US_states[delta_wave==TRUE,t_Delta:=1:.N,by=state]
World[delta_wave==TRUE,t_Delta:=1:.N,by=country]


# Omicron wave classification -----------------------------------------------------------------

### First date with at least mn_growth_time successive days of positive case growth rates.
za_omicron_wave <- function(r,date,min_date=as.Date('2021-11-10'),max_date=as.Date('2021-12-20'),mn_growth_time=10){
  y <- rep(FALSE,length(r))
  
  ix=which(date>min_date & date<max_date & rollapply(r,FUN=function(x) all(x>0),width=mn_growth_time,align='left',fill=NA))
  y[min(ix):length(y)] <- TRUE
  return(y)
}

### US outbreak rule, tricker due to ongoing non-Omicron case rises in some states prior to arrival of Omicron.
### We use Omicron's characteristically faster-than-delta growth to ID Omicron outbreaks, and then start at the 
### most recent date when r<=r_start with r_start=0.02. This rule applies only to outbreaks after mn_date.
us_omicron <- function(r,date,new_confirmed,r_threshold=0.125,r_start=0.02,mn_date=as.Date('2021-11-15')){
  if (!any(r[!is.na(r)]>r_threshold & date[!is.na(r)]>mn_date)){
    y <- rep(FALSE,length(r))
  } else {
    ix <- min(which(r>r_threshold & date>mn_date))
    ix <- max(which(r[1:ix]<r_start))
    y <- (1:length(r))>=ix
  }
  return(y)
}

#### South African omicron wave classification.
ZA[,omicron:=za_omicron_wave(growth_rate,date),by=province]
ZA[omicron==TRUE,t_omicron:=1:.N,by=province]
ZA[,state:=province] ### this makes ggplotting easier later.

#### US state omicron wave classification
US_states[,omicron:=us_omicron(growth_rate,date,exp(mean_position)),by=state]
US_states[omicron==TRUE,t_omicron:=1:.N,by=state]
US_states[,any_omicron:=any(omicron==TRUE),by=state]
US_states[,tt:=date-max(date,na.rm=T)+max(t_omicron,na.rm=T),by=state]

US_states[any_omicron==TRUE & tt>= -21,gr:=nbs(new_confirmed,filtering=TRUE,remove_outliers=FALSE),by=state]
US_states[any_omicron==TRUE & tt>= -21,hgr:=nbs(hosp,filtering=TRUE,remove_outliers=FALSE),by=state]

### Plotting key early US Omicron outbreaks
ggarrange(
ggplot(US_states[state %in% c('District of Columbia','Puerto Rico','Hawaii') & date>as.Date('2021-11-15')],aes(date,new_confirmed))+
  geom_bar(stat='identity',aes(color=omicron,fill=omicron))+
  facet_wrap(.~state,scales='free_y'),
ggplot(US_states[state %in% c('District of Columbia','Puerto Rico','Hawaii') & date>as.Date('2021-11-15')],aes(date,gr))+
  geom_line(aes(color=omicron))+
  facet_wrap(.~state,scales='free_y'),
nrow=2,align='v')

# Delta outbreak comparisons  -------------------------------------------------------------------

### The function below, first_peak, is used to label the dates of first peak for plotting + labelling.
### "peak" for Delta outbreaks was defined by first point at which case growth rates were negative for at least 7 days.
first_peak <- function(t,r,n_min=7){
  ids <- data.table('ix'=1:length(r),'sgn'=sign(r),'id'=rleid(sign(r)))
  ids[,n:=.N,by=id]
  t[ids[sgn<0 & n>=n_min,min(ix)]] %>% return
}

World[delta_wave==TRUE,first_delta_peak:=first_peak(t_Delta,growth_rate),by=country]
US_states[delta_wave==TRUE,first_delta_peak:=first_peak(t_Delta,growth_rate),by=state]


delta_cls <- viridis::plasma(4)
g_delta_early <- ggplot(World[country %in% c("United Kingdom",'India','South Africa') & t_Delta<(first_delta_peak+5)],aes(t_Delta,growth_rate))+
  geom_line(aes(color=country),lwd=2)+
  scale_y_continuous('Growth Rate',limits=c(-0.05,0.4))+
  scale_x_continuous('Outbreak Time',limits=c(0,90))+
  geom_hline(yintercept=0)+
  annotate(geom='segment',x = 80,y=0.12,
           xend = World[country=='India' & t_Delta==first_delta_peak,t_Delta],
           yend=World[country=='India' & t_Delta==first_delta_peak,growth_rate],col='black')+
  annotate(geom='label',x=80,y=0.12,label=paste('India peak: \n',World[country=='India' & t_Delta==first_delta_peak,date]),
           col='white',fill=delta_cls[1],size=5)+
  
  annotate(geom='segment',x = 55,y=0.25,
           xend = World[country=='United Kingdom' & t_Delta==first_delta_peak,t_Delta],
           yend=World[country=='United Kingdom' & t_Delta==first_delta_peak,growth_rate],col='black')+
  annotate(geom='label',x=55,y=0.25,label=paste('UK peak: \n',World[country=='United Kingdom' & t_Delta==first_delta_peak,date]),
           col='white',fill=delta_cls[3],size=5)+
  
  annotate(geom='segment',x = 30,y=0.15,
           xend = World[country=='South Africa' & t_Delta==first_delta_peak,t_Delta],
           yend=World[country=='South Africa' & t_Delta==first_delta_peak,growth_rate],col='black')+
  annotate(geom='label',x=30,y=0.15,label=paste('South Africa peak: \n',World[country=='South Africa' & t_Delta==first_delta_peak,date]),
           col='white',fill=delta_cls[2],size=5)+
  scale_color_manual(values=delta_cls)+
  theme(legend.position=c(0.15,0.8))+
  ggtitle('Delta Outbreaks: Early Archetypes')
  
us_cls <- c('darkred','darkgreen')
g_delta_us <- ggplot(US_states[delta_wave==TRUE & t_Delta <90 & !state %in% c('Puerto Rico','District of Columbia')],aes(t_Delta,growth_rate))+
  geom_line(data=World[country %in% c("United Kingdom",'India','South Africa') & t_Delta<(first_delta_peak+4)],aes(color=country),lwd=2,alpha=0.6)+
  scale_color_manual(values=c(us_cls[1],delta_cls[1],us_cls[2],delta_cls[2:3]))+
  scale_y_continuous('Growth Rate',limits=c(-0.05,0.4))+
  scale_x_continuous('Outbreak Time',limits=c(0,90))+
  geom_line(aes(group=state),lwd=2,col=rgb(0,0,0,0.1))+
  theme(legend.position='none')+
  geom_line(data=US_states[delta_wave==TRUE & t_Delta<90 & state %in% c('Arkansas','Missouri')],aes(color=state),lwd=2.5)+
  geom_hline(yintercept=0)+
  
  annotate(geom='segment',x = 80,y=0.2,
           xend = US_states[state=='Arkansas' & t_Delta==first_delta_peak,t_Delta],
           yend=US_states[state=='Arkansas' & t_Delta==first_delta_peak,growth_rate],col='black')+
  annotate(geom='label',x=80,y=0.2,label=paste('Arkansas peak: \n',US_states[state=='Arkansas' & t_Delta==first_delta_peak,date]),
           col='white',fill=us_cls[1],size=5)+
  
  annotate(geom='segment',x = 30,y=.3,
           xend = US_states[state=='Missouri' & t_Delta==first_delta_peak,t_Delta],
           yend=US_states[state=='Missouri' & t_Delta==first_delta_peak,growth_rate],col='black')+
  annotate(geom='label',x=30,y=.3,label=paste('Missouri peak: \n',US_states[state=='Missouri' & t_Delta==first_delta_peak,date]),
           col='white',fill=us_cls[2],size=5)+
  ggtitle('Delta Outbreaks: US States')


# Omicron outbreak comparisons --------------------------------------------

g_omicron_za <- ggplot(ZA[omicron==TRUE & date<as.Date('2021-12-17')],aes(t_omicron,growth_rate))+
    geom_line(lwd=2,col='orange',alpha=0.3,aes(group=province))+
    geom_line(data=ZA[omicron==TRUE & province=='Gauteng' & date<as.Date('2021-12-17')],col='brown1',lwd=2)+
    geom_smooth(col='orange',lwd=2)+
    ggtitle('Omicron Outreaks: RSA Provinces Prior to 12/21/25')+
    scale_y_continuous('Growth Rate',limits=c(-0.1,0.4))+
    geom_hline(yintercept = 0)+
    annotate(geom='segment',x = 62,y=0.25,
             xend = ZA[province=='Gauteng' & date==as.Date('2021-12-16'),t_omicron],
             yend=ZA[province=='Gauteng' & date==as.Date('2021-12-16'),growth_rate],col='brown1')+
    annotate(geom='label',x=62,y=0.25,label='Gauteng peak: \n 12/16/21',col='white',fill='brown1',size=5)+
  scale_x_continuous('Outbreak Time',limits=c(0,90))
state_cls <- c('purple','darkgreen')
g_omicron_us <- ggplot(US_states[omicron==TRUE],aes(t_omicron,gr))+
    geom_smooth(data=ZA[omicron==TRUE & date<as.Date('2021-12-25')],aes(y=growth_rate),col='orange',lwd=2)+
    geom_line(col='darkgrey',lwd=2,alpha=0.4,aes(group=state))+
    geom_smooth(col='darkgrey',lwd=2)+
    # geom_line(aes(y=hgr,group=state),col='steelblue',lwd=2,alpha=0.4)+
    # geom_smooth(data=US_states[omicron==TRUE & t_omicron>10],aes(y=hgr),col='steelblue',lwd=2)+
    ggtitle('Omicron Outbreaks: US States')+
    scale_y_continuous('Growth Rate',limits=c(-0.1,0.4))+
    geom_point(data=US_states[state %in% c('District of Columbia','Puerto Rico') & t_omicron==max(t_omicron,na.rm=T)])+
    geom_hline(yintercept=0)+
    geom_line(data=US_states[state %in% c('District of Columbia','Puerto Rico')],lwd=2,aes(group=state,color=state))+
  
    annotate(geom='segment',x = 37,y=0.3,
             xend = US_states[state=='District of Columbia' & date==as.Date('2022-01-05'),t_omicron],
             yend=US_states[state=='District of Columbia' & date==as.Date('2022-01-05'),gr],col='black')+
    annotate(geom='label',x=37,y=0.3,label='DC peak: \n 1/5/22',col='white',fill=state_cls[1],size=5)+
  
  annotate(geom='segment',x = 50,y=0.15,
           xend = US_states[state=='Puerto Rico' & date==as.Date('2022-01-05'),t_omicron],
           yend=US_states[state=='Puerto Rico' & date==as.Date('2022-01-05'),gr],col='black')+
  annotate(geom='label',x=50,y=0.15,label='PR peak: \n 1/5/22',col='white',fill=state_cls[2],size=5)+
  scale_color_manual(values=state_cls)+
  theme(legend.position='none')+
  scale_x_continuous('Outbreak Time',limits=c(0,90))
ggarrange(
ggarrange(g_delta_early,g_delta_us,ncol=2,align='h',labels=c('A','B')),
ggarrange(g_omicron_za,g_omicron_us,ncol=2,align='h',labels=c('C','D')),nrow=2)
ggsave('figures/voc_outbreak_durations.png',height=10,width=14)


# Comparing Delta & Omicron -----------------------------------------------

#### In the script below, we combine UK, India, South Africa and US states for Delta outbreaks
#### and we combine South African provinces + US states for Omicron outbreaks.
#### We then compute mean +/- 2sd growth rates as a function of outbreak time & plot to compare Delta vs. Omicron.
ZA[,gr:=growth_rate]
Omics <- rbind(ZA[omicron==TRUE & date<as.Date('2022-12-20'),c('t_omicron','gr','id')][t_omicron<40],
               US_states[omicron==TRUE,c('t_omicron','gr','id')])
Omics[,growth_rate:=gr]
Omics[,outbreak:='Omicron']
Omics[,time:=t_omicron]
Deltas <- rbind(World[t_Delta<(first_delta_peak+5),c('t_Delta','growth_rate','id')],
                US_states[t_Delta<=90,c('t_Delta','growth_rate','id')])
Deltas[,outbreak:='Delta']
Deltas[,time:=t_Delta]

DD <- rbind(Omics[,c('time','growth_rate','outbreak','id')],Deltas[,c('time','growth_rate','outbreak','id')])

d <- DD[,list('growth_rate'=mean(growth_rate,na.rm=T),
              'sd'=sd(growth_rate,na.rm=T)),by=c('outbreak','time')]

trnd=ggplot(data=d[!is.na(sd)],aes(time,growth_rate,color=outbreak))+
  geom_hline(yintercept = 0)+
  geom_line(lwd=2)+
  geom_ribbon(aes(fill=outbreak,ymin=growth_rate-2*sd,ymax=growth_rate+2*sd),alpha=0.3)+
  theme(legend.position=c(0.5,.8))+
  ggtitle('Omicron vs Delta Outbreak Comparison')+
  scale_x_continuous(limits=c(0,90))



# Combining all plots ----------------------------------------------------
ggarrange(
  ggarrange(g_delta_early,g_delta_us,ncol=2,align='h',labels=c('A','B')),
  ggarrange(g_omicron_za,g_omicron_us,ncol=2,align='h',labels=c('C','D')),
  trnd,nrow=3,labels=c(NA,NA,'E'))
ggsave('figures/voc_outbreak_durations_2.png',height=16,width=14)

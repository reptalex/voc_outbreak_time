library(gganimate)
library(ggplot2)
library(data.table)
library(magrittr)
library(plyr)
library(COVID19)
library(parallel)
library(ggpubr)
source('scripts/delta_omicron_utils.R')
theme_set(theme_bw(base_size=15))

ncores <- 6

DD <- readRDS('data/outbreak_data_table.Rds')


p=ggplot(DD[outbreak=='Omicron' & ! state %in% c('California','Rhode Island')],
         aes(time,growth_rate,group=state,frame=date,color=country))+
  transition_reveal(date)+
  geom_line(aes(size=country,alpha=country))+
  enter_fade()+
  scale_x_continuous('Outbreak Time',limits=c(0,40))+
  exit_fade()+
  scale_y_continuous('Growth Rate',limits=c(-0.15,0.4))+
  ggtitle('Omicron Outbreak Tracker | Date: {frame_along}')+
  scale_size_manual(values=c(2,1))+
  scale_alpha_manual(values=c(0.5,0.8))+
  scale_color_manual(values=c('forestgreen','black'))+
  theme(legend.position=c(0.8,0.8))+
  geom_hline(yintercept=0,color='darkgrey')

animate(p,nframes = length(unique(DD$date)),fps = 16) %>%
  anim_save(filename='figures/Omicron_outbreak.gif',animation=.)


ggplot(DD[outbreak=='Omicron' & ! state %in% c('California','Rhode Island')],
       aes(time,growth_rate,color=country))+
  geom_line(aes(size=country,alpha=country,group=state))+
  geom_smooth()+
  scale_x_continuous('Outbreak Time',limits=c(0,40))+
  scale_y_continuous('Growth Rate',limits=c(-0.15,0.4))+
  ggtitle('Omicron Outbreaks')+
  scale_size_manual(values=c(2,1))+
  scale_alpha_manual(values=c(0.1,0.2))+
  scale_color_manual(values=c('forestgreen','black'))+
  theme(legend.position=c(0.8,0.8))+
  geom_hline(yintercept=0,color='darkgrey')
ggsave('figures/Omicron_outbreak.png')



# smoothed trajectories ---------------------------------------------------


smooth <- function(y,x,...) gam(y~s(x,bs='cr'),...)$fitted.values
DD[,r_smooth:=smooth(growth_rate,time),by=id]

# DD[,r_smooth:=nbs(new_confirmed,filtering=FALSE,dispersion=2),by=id]


ggplot(DD[outbreak=='Omicron' & ! state %in% c('California','Rhode Island')],
       aes(time,r_smooth,color=country))+
  geom_line(aes(size=country,alpha=country,group=state))+
  geom_smooth()+
  scale_x_continuous('Outbreak Time',limits=c(0,40))+
  scale_y_continuous('Growth Rate',limits=c(-0.15,0.4))+
  ggtitle('Omicron Outbreaks')+
  scale_size_manual(values=c(2,1))+
  scale_alpha_manual(values=c(0.1,0.2))+
  scale_color_manual(values=c('forestgreen','black'))+
  theme(legend.position=c(0.8,0.8))+
  geom_hline(yintercept=0,color='darkgrey')

ggsave('figures/gifs/Omicron_outbreak_smooth.png')



p_smooth=ggplot(DD[outbreak=='Omicron' & ! state %in% c('California','Rhode Island')],
         aes(time,r_smooth,group=state,frame=date,color=country))+
  transition_reveal(date)+
  geom_line(aes(size=country,alpha=country))+
  enter_fade()+
  scale_x_continuous('Outbreak Time',limits=c(0,40))+
  exit_fade()+
  scale_y_continuous('Case Exponential Growth Rate',limits=c(-0.15,0.4))+
  ggtitle('Omicron Outbreak Tracker | Date: {frame_along}')+
  scale_size_manual(values=c(2,1))+
  scale_alpha_manual(values=c(0.5,0.8))+
  scale_color_manual(values=c('forestgreen','black'))+
  theme(legend.position=c(0.8,0.8))+
  geom_hline(yintercept=0,color='darkgrey')

animate(p_smooth,nframes = length(unique(DD$date)),fps = 16) %>%
  anim_save(filename='figures/gifs/Omicron_outbreak_smooth.gif',animation=.)


# State gif-making function ------------------------------------------------------------

gif_my_states <- function(sts,filename=NULL,DD.=DD){
  p_st <- ggplot(DD[outbreak=='Omicron' & 
                      (country=='South Africa'|state %in% sts)],
         aes(time,growth_rate,group=state,frame=date,color=country))+
    transition_reveal(date)+
    geom_line(aes(size=country,alpha=country))+
    geom_text(aes(label=state),subset= .(state %in% sts))+
    enter_fade()+
    scale_x_continuous('Outbreak Time',limits=c(0,40))+
    exit_fade()+
    scale_y_continuous('Case Exponential Growth Rate',limits=c(-0.15,0.4))+
    ggtitle('Omicron Outbreak Tracker | Date: {frame_along}')+
    scale_size_manual(values=c(2,2))+
    scale_alpha_manual(values=c(0.5,0.8))+
    scale_color_manual(values=c('forestgreen','black'))+
    theme(legend.position=c(0.8,0.8))+
    geom_hline(yintercept=0,color='darkgrey')
  
  if (is.null(filename)){
    filename=paste('figures/gifs/',paste(gsub(' ','-',sts),collapse='_'),'_Omicron_outbreak_smooth.gif')
  }
  animate(p_st,nframes = length(unique(DD$date)),fps = 16) %>%
    anim_save(filename=filename,animation=.)
}



# US state gifs -----------------------------------------------------------


gif_my_states('District of Columbia')
gif_my_states(c('Florida','Puerto Rico','Hawaii')) ### early US outbreaks
gif_my_states(c('Florida','Puerto Rico','Hawaii','District of Columbia'))


# US_counties -------------------------------------------------------------
county_population_threshold <- 1e6

cl <- makeCluster(ncores)
clusterEvalQ(cl,expr = {library(tsoutliers)
  library(magrittr)
  library(tidyverse)})
clusterExport(cl,'outlier_detection')

thanksgiving <- c(seq(as.Date('2020-11-26'),as.Date('2020-11-27'),by='day'),
                  seq(as.Date('2021-11-25'),as.Date('2021-11-29'),by='day'))
christmas <- as.Date(c('2020-12-25','2020-12-26'))
ny<- c(seq(as.Date('2021-01-01'),as.Date('2021-01-02'),by='day'),
       seq(as.Date('2022-01-01'),as.Date('2022-01-02'),by='day'))
lbr <- as.Date(c('2021-09-06','2021-09-07'))
ipd <- as.Date(c('2021-10-11','2021-10-12'))
hlwn <- as.Date(c('2021-10-30','2021-10-31','2021-11-01'))
ind_day <- seq(as.Date('2021-07-03'),as.Date('2021-07-06'),by='day')

US_counties <- covid19("United States",level=3) %>% as.data.table
US_counties <- US_counties[population>county_population_threshold]
US_counties[,country:=administrative_area_level_1]
US_counties[,state:=administrative_area_level_2]
US_counties[,county:=administrative_area_level_3]
US_counties <- US_counties[!state %in% c('Guam','Northern Mariana Islands','Virgin Islands','American Samoa')]
US_counties[,new_confirmed:=dfs(confirmed),by=c('state','county')]
US_counties[,new_deaths:=dfs(deaths),by=c('state','county')]
US_counties[,lbl:=paste(state,county,sep=',')]
US_counties <- US_counties[population>=county_population_threshold]

US_counties[date %in% c(thanksgiving,christmas,ny,lbr,ipd,hlwn,ind_day),new_confirmed:=NA]
US_counties <- US_counties[date>as.Date('2021-10-01')]
US_counties[,new_confirmed:=outlier_detection_par(new_confirmed,lbl,cl)]


stopCluster(cl)
rm('cl')
gc()

US_counties[new_confirmed==0,new_confirmed:=NA]

US_counties[,growth_rate:=nbs(new_confirmed,filtering=TRUE,dispersion=2),by=c('state','county')]
US_counties <- US_counties[,c('date','new_confirmed','new_deaths','country','state','county',
                              'growth_rate','id','lbl')]


us_omicron <- function(r,date,new_confirmed,r_threshold=0.125,r_start=0.027,mn_date=as.Date('2021-11-23')){
  if (!any(r[!is.na(r)]>r_threshold & date[!is.na(r)]>mn_date)){
    y <- rep(FALSE,length(r))
  } else {
    ix <- min(which(r>r_threshold & date>mn_date))
    ix <- max(which(r[1:ix]<r_start))
    y <- (1:length(r))>=ix
  }
  return(y)
}

US_counties[,omicron:=us_omicron(growth_rate,date,new_confimred),by=id]
US_counties[omicron==TRUE,time:=1:.N,by=id]

saveRDS(US_counties,file = 'data/populous_US_counties_data_table.Rds')
write.csv(US_counties,file = 'data/populous_US_counties_data_table.csv')



# County visualization functions ----------------------------------------------------
plot_state <- function(x,v=NA,US_states.=US_states){
  ggarrange(
    ggplot(US_states[state==x],aes(date,new_confirmed))+
      geom_bar(stat='identity',aes(color=omicron,fill=omicron))+
      geom_vline(xintercept = v)+
      ggtitle(paste(x,'new cases')),
    ggplot(US_states[state==x],aes(date,growth_rate))+
      geom_line(lwd=2,aes(color=omicron))+
      geom_vline(xintercept = v)+
      ggtitle(paste(x,'r(t)')),
    nrow=2,align='v') %>%
    return()
}
plot_county <- function(x,v=NA,US_counties.=US_counties){
  ggarrange(
    ggplot(US_counties[county==x],aes(date,new_confirmed))+
      geom_bar(stat='identity',aes(color=omicron,fill=omicron))+
      geom_vline(xintercept = v)+
      ggtitle(paste(x,'county new cases')),
    ggplot(US_counties[county==x],aes(date,growth_rate))+
      geom_line(lwd=2,aes(color=omicron))+
      geom_vline(xintercept = v)+
      ggtitle(paste(x,'county r(t)')),
      nrow=2,align='v') %>%
    return()
}


gif_my_county <- function(counties,county_labels=NULL,states=NULL,
                          filename=NULL,omicron_start_date=NULL,
                          DD.=DD,US_counties.=US_counties,fps=16){
  if (!is.null(county_labels)){
    if(length(county_labels)!=length(counties)){
      stop('county_labels must be same length as counties')
    }
  } else {
    county_labels <- counties
  }
  ucs <- US_counties[(county %in% counties)]
  ucs[,label:=county_labels[match(county,counties)],by=county]
  ucs[,country:='USA']
  
  if (!is.null(omicron_start_date)){
    ucs$omicron <- NULL
    ucs[,omicron:=date>=as.Date(omicron_start_date)]
    ucs[omicron==TRUE,time:=1:.N]
  }
  ucs <- ucs[omicron==TRUE]
  if (!is.null(states)){
    ucs <- ucs[state %in% states]
  }
  xx <- DD[outbreak=='Omicron' & country=='South Africa']
  xx[,label:=state]
  
  xx <- rbind(xx[,c('time','growth_rate','label','date','country')],
              ucs[,c('time','growth_rate','label','date','country')])
  
  p_county <- ggplot(xx,aes(time,growth_rate,group=label,frame=date,color=country))+
    transition_reveal(date)+
    geom_line(aes(size=country,alpha=country))+
    geom_text(aes(label=label))+
    enter_fade()+
    scale_x_continuous('Outbreak Time',limits=c(0,45))+
    exit_fade()+
    scale_y_continuous('Case Exponential Growth Rate',limits=c(-0.15,0.4))+
    ggtitle('Omicron Outbreak Tracker | Date: {frame_along}')+
    scale_size_manual(values=c(2,2))+
    scale_alpha_manual(values=c(0.5,0.8))+
    scale_color_manual(values=c('forestgreen','black'))+
    theme(legend.position=c(0.8,0.8))+
    geom_hline(yintercept=0,color='darkgrey')
  
  if (is.null(filename)){
    filename=paste('figures/gifs/counties_',paste(gsub(' ','-',counties),collapse='_'),'_Omicron_outbreak.gif')
  }
  animate(p_county,nframes = length(unique(xx$date)),fps = fps) %>%
    anim_save(filename=filename,animation=.)
}



# County gifs -------------------------------------------------------------

gif_my_county('New York City',fps=6,US_counties=US_counties)
gif_my_county(c('New York City','Montgomery','Orange'),
              county_labels = c('NYC','Montgomery, MD','Orange, FL'),
              states=c('New York','Maryland','Florida'))

### the following counties don't have Omicron waves classified by our
### fast-growth classifier above. We'll chose start-dates
### as the approximate date where their recent growth rate rises would 
### project backwards to zero.

plot_county('Cook',v=as.Date('2021-12-11'))
plot_county('Middlesex',v=as.Date('2021-12-15'))

gif_my_county('Cook',omicron_start_date = '2021-12-11',fps=6)
gif_my_county('New York City',fps=6)


gif_my_county(c('Cook','Middlesex'),omicron_start_date = '2021-12-11',fps=6)


US_counties[county=='Cook',omicron:=date>=as.Date('2021-12-11')]
US_counties[county=='Cook' & omicron==TRUE,time:=1:.N]
US_counties[county=='Middlesex',omicron:=date>=as.Date('2021-12-15')]
US_counties[county=='Middlesex' & omicron==TRUE,time:=1:.N]

gif_my_county(c('New York City','Cook','Middlesex'),
              county_labels=c('NYC','Boston','Chicago'))



# Delta gif ---------------------------------------------------------------

DD[country=='USA',label:=state]
DD[country!='USA' & outbreak=='Delta',label:=country]
DD[country=='South Africa' & outbreak=='Omicron',label:=state]

delta_gif_my_states <- function(sts,filename=NULL,DD.=DD){
  p_st <- ggplot(DD[outbreak=='Delta' & 
                      (country %in% c('United Kingdom','India','South Africa') |
                         state %in% sts)],
                 aes(time,growth_rate,group=label,frame=date,color=country))+
    transition_reveal(date)+
    geom_line(size=2)+
    geom_text(aes(label=label))+
    enter_fade()+
    scale_x_continuous('Outbreak Time',limits=c(0,90))+
    exit_fade()+
    scale_y_continuous('Case Exponential Growth Rate',limits=c(-0.15,0.4))+
    ggtitle('Delta Outbreak Tracker | Date: {frame_along}')+
    scale_size_manual(values=c(2,2))+
    scale_color_manual(values=c('#FF7722','#007a4d','#00247d','black'))+
    theme(legend.position=c(0.8,0.8))+
    geom_hline(yintercept=0,color='darkgrey')
  
  ### colors: India, South Africa, UK, USA
  
  if (is.null(filename)){
    filename=paste('figures/gifs/Delta_',paste(gsub(' ','-',sts),collapse='_'),'outbreak.gif')
  }
  animate(p_st,nframes = length(unique(DD$date)),fps = 16) %>%
    anim_save(filename=filename,animation=.)
}
#devtools::install_github("m-clark/visibly")
library(visibly)
library(mgcv)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotly)
library(here)

##source log-normal CI function
source(here('R','log-normCIs.R'))

##get a version of the data:
dat=read.csv(here('data','TT_align_allData_13102025.csv'))

##processing
#times
dat$timeStamp=as.POSIXct(dat$dt,format='%d/%m/%Y %H:%M', tz='UTC')
dat$hour=as.numeric(format(dat$timeStamp, "%H"))
dat$day_of_year=as.numeric(format(dat$timeStamp, "%j"))
summary(dat$day_of_year)

#krill presence-absence
dat$krill_PA=0
dat$krill_PA[dat$rho_v>0]=1
table(dat$krill_PA)
#minimum detectable krill density:
min_det=min(dat$rho_v[dat$krill_PA==1])
dat$rho_min=dat$rho_v+min_det/2

#deal with centre of mass NaN and assign a 'background value':
dat$Center_of_mass_min=200
dat$Center_of_mass_min[which(!is.nan(dat$Center_of_mass))]=
  dat$Center_of_mass[which(!is.nan(dat$Center_of_mass))]
plot(dat$Center_of_mass_min)

#deal with Corrected length NaN and negative values:
summary(dat$Corrected_length)
dat$Corrected_length[which(dat$Corrected_length<0)]=NA
dat$Corrected_length[which(is.nan(dat$Corrected_length))]=NA
min_Corrected_length=min(dat$Corrected_length,na.rm=T)
#assign a background value for Corrected length:
dat$Corrected_length_min=dat$Corrected_length+min_Corrected_length/2

#deal with Proportion occupied NaN and negative values:
dat$Proportion_occupied[which(is.nan(dat$Proportion_occupied))]=NA
#assign a background value for Proportion occupied:
Proportion_occupied_min=min(dat$Proportion_occupied,na.rm=T)
dat$Proportion_occupied[is.na(dat$Proportion_occupied)]=Proportion_occupied_min/2

##Modelling
#gam for "Bm_Ant_ABZ"
summary(dat$Bm_Ant_ABZ)
m1 <- gam(Bm_Ant_ABZ ~  
            s(hour, k=15,bs='cc') + 
            s(day_of_year, k=15,bs='cc')+
            s(solar_lux,k=6,bs='cr')+
            s(lunar_lux,k=6,bs='cr')+
            s(ice,k=6,bs='cr'),
          data=dat,
          method='ML',
          family=gaussian(link='log'),
          gamma=1.4,
          select=TRUE,
          knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(m1)
plot(m1)
gam.check(m1)
concurvity(m1)

m1a <- gam(Bm_Ant_ABZ ~  
            s(hour, k=15,bs='cc') + 
            s(day_of_year, k=20,bs='cc')+
            #s(solar_lux,k=6,bs='cr')+
            s(lunar_lux,k=6,bs='cr')+
            s(ice,k=6,bs='cr'),
          data=dat,
          method='ML',
          family=gaussian(link='log'),
          gamma=1.4,
          select=TRUE,
          knots=list(hour=c(0,24),day_of_year=c(0,365)))
concurvity(m1a)
gam.check(m1a)
visibly::plot_gam_check(m1a)
#gam for "Bm_D"

m2 <- gam(Bm_D ~  
            s(hour, k=15,bs='cc') + 
            s(day_of_year, k=15,bs='cc')+
            s(solar_lux,k=6,bs='cr')+
            s(lunar_lux,k=6,bs='cr')+
            s(ice,k=6,bs='cr'),
          data=dat,
          method='ML',
          family=gaussian(link='log'),
          gamma=1.4,
          select=TRUE,
          knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(m2)
concurvity(m2)

m2a <- gam(Bm_D ~  
            s(hour, k=15,bs='cc') + 
            s(day_of_year, k=15,bs='cc')+
            #s(solar_lux,k=6,bs='cr')+
            s(lunar_lux,k=6,bs='cr')+
            s(ice,k=6,bs='cr'),
          data=dat,
          method='ML',
          family=gaussian(link='log'),
          gamma=1.4,
          select=TRUE,
          knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(m2a)
concurvity(m2a)
gam.check(m2a)
#gam for "BD"
m3 <- gam(BD ~  
            s(hour, k=15,bs='cc') + 
            s(day_of_year, k=15,bs='cc')+
            s(solar_lux,k=6,bs='cr')+
            s(lunar_lux,k=6,bs='cr')+
            s(ice,k=6,bs='cr'),
          data=dat,
          method='ML',
          family=gaussian(link='log'),
          gamma=1.4,
          select=TRUE,
          knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(m3)
concurvity(m3)
gam.check(m3)
m3a <- gam(BD ~  
            s(hour, k=15,bs='cc') + 
            s(day_of_year, k=15,bs='cc')+
            #s(solar_lux,k=6,bs='cr')+
            s(lunar_lux,k=6,bs='cr')+
            s(ice,k=6,bs='cr'),
          data=dat,
          method='ML',
          family=gaussian(link='log'),
          gamma=1.4,
          select=TRUE,
          knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(m3a)
concurvity(m3a)
gam.check(m3a)

#predictions for the call models (m1,m2 and  m3)
preds=data.frame(hour=mean(dat$hour),
                 day_of_year=0:365,
                 solar_lux=mean(dat$solar_lux,na.rm=T),
                 lunar_lux=mean(dat$lunar_lux,na.rm=T),
                 ice=mean(dat$ice,na.rm=T))
preds_m1=as.data.frame(predict(m1a,newdata=preds,type='response',se.fit=TRUE))
preds_m1$model='Bm_Ant_ABZ'
preds_m2=as.data.frame(predict(m2a,newdata=preds,type='response',se.fit=TRUE))
preds_m2$model='Bm_D'
preds_m3=as.data.frame(predict(m3a,newdata=preds,type='response',se.fit=TRUE))
preds_m3$model='BD'
preds_all=rbind(preds_m1,preds_m2,preds_m3)
preds_all$day_of_year=preds$day_of_year
preds_all$lower=preds_all$fit-2*preds_all$se.fit
preds_all$upper=preds_all$fit+2*preds_all$se.fit

#plot predictions of the call models:
ggplot(data=preds_all,aes(x=day_of_year,y=fit,colour=model,group=model))+
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=model),alpha=0.4,colour=NA)+
  theme_minimal(base_size = 15)+
  labs(y='Predicted call rate, nbr/hr',x='Day of year')+
  theme(legend.position = 'top')
ggsave(here('plots','callsTimeSeries.png'))

##Krill model:
mkrill= gam(rho_min ~  
            s(hour, k=15,bs='cc') + 
            s(day_of_year, k=15,bs='cc')+
            s(solar_lux,k=6,bs='cr')+
            s(lunar_lux,k=6,bs='cr')+
            s(ice,k=6,bs='cr')+
            s(Center_of_mass_min,k=6,bs='cr')+
            s(Corrected_length_min,k=6,bs='cr')+
            s(Proportion_occupied,k=3,bs='cr'),
          data=dat,
          method='ML',
          family=Gamma(link='log'),
          gamma=1.4,
          select=TRUE,
          knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(mkrill)
gam.check(mkrill)
concurvity(mkrill)

mkrilla= gam(rho_min ~  
              s(hour, k=15,bs='cc') + 
              s(day_of_year, k=15,bs='cc')+
              #s(solar_lux,k=6,bs='cr')+
              s(lunar_lux,k=6,bs='cr')+
              s(ice,k=6,bs='cr')+
              s(Center_of_mass_min,k=6,bs='cr')+
              s(Corrected_length_min,k=6,bs='cr')+
              s(Proportion_occupied,k=3,bs='cr'),
            data=dat,
            method='ML',
            family=Gamma(link='log'),
            gamma=1.4,
            select=TRUE,
            knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(mkrilla)
gam.check(mkrilla)
concurvity(mkrilla)

mkrillb= gam(rho_min ~  
               s(day_of_year,hour, k=15) + 
               #s(day_of_year, k=15,bs='cc')+
               #s(solar_lux,k=6,bs='cr')+
               s(lunar_lux,k=6,bs='cr')+
               s(ice,k=6,bs='cr')+
               s(Center_of_mass_min,k=6,bs='cr')+
               s(Corrected_length_min,k=6,bs='cr')+
               s(Proportion_occupied,k=3,bs='cr'),
             data=dat,
             method='ML',
             family=Gamma(link='log'),
             gamma=1.4,
             select=TRUE,
             knots=list(hour=c(0,24),day_of_year=c(0,365)))
summary(mkrillb)

preds_k=data.frame(hour=mean(dat$hour),
                   day_of_year=0:365,
                   solar_lux=mean(dat$solar_lux,na.rm=T),
                   lunar_lux=mean(dat$lunar_lux,na.rm=T),
                   ice=mean(dat$ice,na.rm=T),
                   Center_of_mass_min=mean(dat$Center_of_mass_min,na.rm=T),
                   Corrected_length_min=mean(dat$Corrected_length_min,na.rm=T),
                   Proportion_occupied=mean(dat$Proportion_occupied,na.rm=T))
preds_k_all=as.data.frame(predict(mkrillb,newdata=preds_k,type='response',se.fit=TRUE))
preds_k_all$model='krill_density'
preds_k_all$day_of_year=preds_k$day_of_year


preds_k_all <- cbind(preds_k_all,
                     lognormal_ci_orig(mean=preds_k_all$fit, se=preds_k_all$se.fit))

head(preds_all)
head(preds_k_all)
predictions_df=rbind(preds_all,preds_k_all)

head(predictions_df)

ggplot(data=predictions_df,aes(x=day_of_year,y=fit,colour=model,group=model))+
  geom_line(size=1) + 
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=model),alpha=0.4,colour=NA)+ 
  theme_minimal(base_size = 15)+
  labs(y='Predicted cue rate (nbr/hr) and krill density (g/m**3)',x='Day of year')+
  theme(legend.position = 'top')
ggsave(here('plots','calls_krillTimeSeries.png'))
########################
##PCA
df_pca <- predictions_df %>%  dplyr::select(day_of_year, model, fit) %>%
  tidyr::pivot_wider(
    names_from = model,       # column whose values become new column names
    values_from = fit  # column whose values fill those new columns
  )

summary(df_pca)
pca <- prcomp(df_pca[, -1], center = TRUE, scale. = TRUE)

# Extract PC scores (observations)
pc_df <- as.data.frame(pca$x[, 1:3])
pc_df$Sample <- rownames(df_pca)

# Extract loadings (eigenvectors)
loadings <- as.data.frame(pca$rotation[, 1:3])

# Variance explained
var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
pc_labels <- paste0("PC", 1:3, " (", round(var_explained[1:3], 1), "%)")


# Arrow scale factor for visibility
arrow_scale <- 5

# Compute arrow start/end points
arrow_data <- data.frame(
  var = rownames(loadings),
  x0 = 0, y0 = 0, z0 = 0,
  x1 = loadings[, 1] * arrow_scale,
  y1 = loadings[, 2] * arrow_scale,
  z1 = loadings[, 3] * arrow_scale
)

# Build interactive 3D plot
# Base plot
p <- plot_ly() %>%
  # PCA sample points
  add_markers(
    data = pc_df,
    x = ~PC1, y = ~PC2, z = ~PC3,
    marker = list(size = 5, color = 'steelblue'),
    text = ~Sample, hoverinfo = "text", name = "Samples"
  )
# Add 3D arrows (one per variable)
for (i in seq_len(nrow(arrow_data))) {
  p <- p %>%
    add_trace(
      x = c(arrow_data$x0[i], arrow_data$x1[i]),
      y = c(arrow_data$y0[i], arrow_data$y1[i]),
      z = c(arrow_data$z0[i], arrow_data$z1[i]),
      type = "scatter3d", mode = "lines",
      line = list(color = "red", width = 6),
      name = arrow_data$var[i],
      hoverinfo = "text",
      text = arrow_data$hover[i],
      showlegend = FALSE
    ) %>%
    add_text(
      x = arrow_data$x1[i],
      y = arrow_data$y1[i],
      z = arrow_data$z1[i],
      text = arrow_data$var[i],
      textposition = "top center",
      showlegend = FALSE,
      textfont = list(color = "red", size = 14)
    )
}

# Axis labels and layout
p <- p %>%
  layout(scene = list(
    xaxis = list(title = pc_labels[1]),
    yaxis = list(title = pc_labels[2]),
    zaxis = list(title = pc_labels[3])
  ))

p

#htmlwidgets::saveWidget(p, here('plots','PCA_3D.html'))

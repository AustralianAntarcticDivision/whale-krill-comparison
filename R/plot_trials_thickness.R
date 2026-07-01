##################################################
# Thickness
plot(mthickness,pages=1) #
concurvity(mthickness) #ice should go
mthickness2<-update(mthickness,select=TRUE) #shrinks curves where not important, contribute 0
plot(mthickness2,pages=1) # suggest ice U/S

mthickness2_drop_ice<-update(mthickness2,.~.-s(ice, k = 6, bs = "cr"))
plot(mthickness2_drop_ice,pages=1) #all good.

vars <- c("hour", "day_of_year", "solar_lux", "lunar_lux")
pairs <- combn(vars, 2, simplify = FALSE)


biv_df <- lapply(pairs, function(p) make_bivariate(p[1], p[2],model=mthickness2_drop_ice,dat=dat)) |> bind_rows()

cat("Pairs computed:", length(pairs), "\n")
cat("Total rows:", nrow(biv_df), "\n")

z_rng<-range(biv_df$fit)
library(patchwork)

n_vars <- length(vars)

# Build 5×5 grid: lower triangle = heatmap, diagonal = label, upper = spacer
grid_plots <- vector("list", n_vars ^ 2)
range(biv_df$fit)
for (col in seq_len(n_vars)) {
  for (row in seq_len(n_vars)) {
    idx <- (row - 1) * n_vars + col
    
    if (row == col) {
      grid_plots[[idx]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = short_labels[vars[row]], size = 3.5, fontface = "bold") +
        theme_void()
      
    } else if (row > col) {
      v1 <- vars[col]   # x-axis variable (column)
      v2 <- vars[row]   # y-axis variable (row)
      df <- biv_df |> filter(var1 == v1, var2 == v2)
      
      p <- ggplot(df, aes(x = .data[[v1]], y = .data[[v2]], fill = fit)) +
        geom_raster(interpolate = TRUE) +
        geom_contour(aes(z = fit), colour = "white", alpha = 0.5, linewidth = 0.3) +
        #scale_fill_viridis_c(option = "magma", trans = "log", 
        #                     limits=z_rng ,guide = "none") +
        scale_fill_viridis_c(
          option = "magma",
          trans = "log",
          limits = z_rng,
          name = "Krill thickness",#centre\nof mass",
          labels = scales::label_number(accuracy = 0.1)
        )
      theme_bw(base_size = 7) +
        theme(axis.title  = element_blank(),
              plot.margin = margin(2, 2, 2, 2))
      
      # Suppress axis text except on outer edges
      if (row != n_vars) p <- p + theme(axis.text.x  = element_blank(),
                                        axis.ticks.x = element_blank())
      if (col != 1)      p <- p + theme(axis.text.y  = element_blank(),
                                        axis.ticks.y = element_blank())
      
      grid_plots[[idx]] <- p
      
    } else {
      grid_plots[[idx]] <- plot_spacer()
    }
  }
}

p_final <- wrap_plots(grid_plots, ncol = n_vars) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Pairwise marginal predictions from krill thickness",
    subtitle = "Lower triangle: Krill thickness | others held at their means"
  ) &
  theme(
    legend.position = "right"
  )

p_final
#### add in the whale calls:
AIC_base<-AIC(mthickness2_drop_ice)

#add in Bm_Ant_ABZ
mthickness2_drop_ice_Bm_Ant_ABZ<-update (mthickness2_drop_ice, .~.+s(Bm_Ant_ABZ,bs='cr',k=5))
plot(mthickness2_drop_ice_Bm_Ant_ABZ,pages=1)
summary(mthickness2_drop_ice_Bm_Ant_ABZ)
concurvity(mthickness2_drop_ice_Bm_Ant_ABZ) #0.916 for Bm_Ant_ABZ is really high.
#can't use Bm_Ant_ABZ

#Bm_D: USE ME!

mthickness2_drop_ice_Bm_D<- update (mthickness2_drop_ice, .~.+s(Bm_D,bs='cr',k=5))
plot(mthickness2_drop_ice_Bm_D,pages=1)
summary(mthickness2_drop_ice_Bm_D)
concurvity(mthickness2_drop_ice_Bm_D) #all of these are fine.
AIC(mthickness2_drop_ice_Bm_D) - AIC_base#lets keep Bm_D in AIC -6.6

#BD
mthickness2_drop_ice_BD<- update (mthickness2_drop_ice, .~.+s(BD,bs='cr',k=5))
plot(mthickness2_drop_ice_BD,pages=1)
summary(mthickness2_drop_ice_BD)
concurvity(mthickness2_drop_ice_BD)
AIC(mthickness2_drop_ice_BD) - AIC_base #not much improvement in AIC (-7.7)

#BD + Bm_D
mthickness2_drop_ice_BDplusBm_D<- update (mthickness2_drop_ice_BD, .~.+s(Bm_D,bs='cr',k=5))
plot(mthickness2_drop_ice_BDplusBm_D,pages=1)
summary(mthickness2_drop_ice_BDplusBm_D)
concurvity(mthickness2_drop_ice_BDplusBm_D)
AIC(mthickness2_drop_ice_BDplusBm_D) - AIC_base #its an improvement of  -14.3

#add in BD + BmD into mcom
mthick_plot <- update(mthickness,.~.,select=TRUE)
mthick_plot <- update(mthick_plot, .~.+s(BD,bs='cr',k=5)++s(Bm_D,bs='cr',k=5))
plot(mthick_plot,pages=1)
##############################################################################

vars <- c("hour", "day_of_year",'solar_lux','ice','lunar_lux','Bm_D','BD') #c("hour", "day_of_year", "ice","Bm_D")
pairs <- combn(vars, 2, simplify = FALSE)


biv_df <- lapply(pairs, function(p) make_bivariate(p[1], p[2],model=mthick_plot,dat=dat)) |> bind_rows()

cat("Pairs computed:", length(pairs), "\n")
cat("Total rows:", nrow(biv_df), "\n")

z_rng<-range(biv_df$fit)
library(patchwork)

n_vars <- length(vars)

# Build 4×4 grid: lower triangle = heatmap, diagonal = label, upper = spacer
grid_plots <- vector("list", n_vars ^ 2)
range(biv_df$fit)
for (col in seq_len(n_vars)) {
  for (row in seq_len(n_vars)) {
    idx <- (row - 1) * n_vars + col
    
    if (row == col) {
      grid_plots[[idx]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = short_labels[vars[row]], size = 3.5, fontface = "bold") +
        theme_void()
      
    } else if (row > col) {
      v1 <- vars[col]   # x-axis variable (column)
      v2 <- vars[row]   # y-axis variable (row)
      df <- biv_df |> filter(var1 == v1, var2 == v2)
      
      p <- ggplot(df, aes(x = .data[[v1]], y = .data[[v2]], fill = fit)) +
        geom_raster(interpolate = TRUE) +
        geom_contour(aes(z = fit), colour = "white", alpha = 0.5, linewidth = 0.3) +
        #scale_fill_viridis_c(option = "magma", trans = "log", 
        #                     limits=z_rng ,guide = "none") +
        scale_fill_viridis_c(
          option = "magma",
          trans = "log",
          limits = z_rng,
          name = "Krill thickness",
          labels = scales::label_number(accuracy = 0.1)
        )
      theme_bw(base_size = 7) +
        theme(axis.title  = element_blank(),
              plot.margin = margin(2, 2, 2, 2))
      
      # Suppress axis text except on outer edges
      if (row != n_vars) p <- p + theme(axis.text.x  = element_blank(),
                                        axis.ticks.x = element_blank())
      if (col != 1)      p <- p + theme(axis.text.y  = element_blank(),
                                        axis.ticks.y = element_blank())
      
      grid_plots[[idx]] <- p
      
    } else {
      grid_plots[[idx]] <- plot_spacer()
    }
  }
}

p_final <- wrap_plots(grid_plots, ncol = n_vars) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Pairwise marginal predictions from krill thickness",
    subtitle = "Lower triangle: Krill thickness | others held at their means"
  ) &
  theme(
    legend.position = "right"
  )

p_final

#Compilation for MS figure:
thick_fig_list=list(solar_ice=p_final[[24]]+  theme_bw(),
                  solar_BD=p_final[[45]]+  theme_bw(),
                  solar_BmD=p_final[[38]]+   theme_bw())

thick_collect_plots<-wrap_plots(thick_fig_list, ncol = 3) + plot_layout(guides='collect',
                                                                    axis_titles = 'collect') 


den_collect_plots / thick_collect_plots / com_collect_plots +
  plot_annotation(tag_levels = 'a') &
  theme(text = element_text(size = 12))

ggsave(here::here('plots','density_thickness_and_com_comparison.png'))

##
#look at relative importance (krill centre of mass model)
hp<-NULL
hp <- gam.hp(
  mthick_plot,
  type = "dev"
)

hp
plot(hp, plot.perc = TRUE)
#calculate the proportion explained as a contribution to total deviance 
hp_tab <- as.data.frame(hp$hierarchical.partitioning)

total_dev_expl <- summary(
  mthick_plot
)$dev.expl

hp_tab$contribution_to_total_deviance <-
  total_dev_expl * hp_tab$I.perc / 100

hp_tab$percentage_total_dev=hp_tab$contribution_to_total_deviance *100
rn<-row.names(hp_tab)
f1=function(x){
  if(substr(x,1,1)=='s'){
    out<-strsplit(substr(x,3,nchar(x)),',')[[1]][1]
  } else {
    out=x}
  return(out)}
f2=Vectorize(f1)
hp_tab$var_name<-f2(rn)
#bar plot

hp_tab_thick<-hp_tab[order(hp_tab$percentage,decreasing=T),]
hp_tab_thick$var_nameF=factor(hp_tab_thick$var_name,levels=hp_tab_thick$var_name)
pthick <- ggplot(hp_tab_thick,aes(var_nameF,percentage_total_dev)) + geom_col(col='black',fill='lightgrey') + theme_bw() +
  ylab('Total deviance explained (%)') + xlab('Explanatory variable') + ggtitle("Swarm thickness")+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16)
  )

pdensity / pcom / pthick

ggsave(here::here('plots','krill all model relative contribution.png'))

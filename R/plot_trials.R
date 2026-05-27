[File: 10.1KB, 296 lines, showing lines 1-296]
# Marginal prediction plots for m1
# Each predictor varied over its observed range; all others held at their means.

library(tidyverse)

# Predictor means (used to hold non-focal variables constant)
means <- dat |>
  summarise(across(c(hour, day_of_year, solar_lux, lunar_lux, ice), mean, na.rm = TRUE))

# Build a marginal prediction data frame for a single variable
make_marginal <- function(var, n = 200) {
  seq_vals <- seq(min(dat[[var]], na.rm = TRUE),
                  max(dat[[var]], na.rm = TRUE),
                  length.out = n)
  
  grid <- means[rep(1, n), ]
  grid[[var]] <- seq_vals
  
  pr <- predict(m1, newdata = grid, type = "link", se.fit = TRUE)
  
  # CIs computed on link scale then back-transformed to avoid negative bounds
  data.frame(
    x     = seq_vals,
    fit   = exp(pr$fit),
    lower = exp(pr$fit - 1.96 * pr$se.fit),
    upper = exp(pr$fit + 1.96 * pr$se.fit),
    var   = var
  )
}

vars <- c("hour", "day_of_year", "solar_lux", "lunar_lux", "ice")

var_labels <- c(
  hour        = "Hour of day",
  day_of_year = "Day of year",
  solar_lux   = "Solar illuminance (lux)",
  lunar_lux   = "Lunar illuminance (lux)",
  ice         = "Ice concentration"
)

marg_df <- lapply(vars, make_marginal) |> bind_rows()

# Rug data: observed values of each predictor
rug_df <- dat |>
  pivot_longer(all_of(vars), names_to = "var", values_to = "x") |>
  mutate(var = factor(var, levels = vars))

# Austral season boundaries and midpoint labels for the day_of_year panel
season_vlines <- data.frame(
  var   = "day_of_year",
  x     = season_starts[2:5]  # transition days: 60, 152, 244, 335
)

season_text <- data.frame(
  var   = "day_of_year",
  x     = c(30, 106, 198, 290, 350),
  label = c("Sum.", "Aut.", "Win.", "Spr.", "Sum.")
)

marg_df |>
  mutate(var = factor(var, levels = vars)) |>
  ggplot(aes(x, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line() +
  geom_rug(data = rug_df, aes(x = x), sides = "b", alpha = 0.05,
           inherit.aes = FALSE) +
  geom_vline(data = season_vlines, aes(xintercept = x),
             linetype = "dashed", colour = "steelblue", linewidth = 0.5) +
  geom_text(data = season_text, aes(x = x, label = label),
            y = Inf, vjust = 1.4, size = 2.8, colour = "steelblue",
            inherit.aes = FALSE) +
  facet_wrap(~var, scales = "free",
             labeller = labeller(var = var_labels)) +
  labs(
    x        = NULL,
    y        = "Predicted Bm_Ant_ABZ",
    title    = "Marginal predictions from m1",
    subtitle = "All other predictors held at their mean"
  )

# --- Day-of-year panel with calendar date x-axis and observed monthly means --
# Reference year 2021 (non-leap); day 1 = Jan 1 2021.
doy_to_date <- function(doy) as.Date(doy - 1, origin = "2021-01-01")

doy_dates <- marg_df |>
  filter(var == "day_of_year") |>
  mutate(date = doy_to_date(x))

season_vlines_date <- data.frame(date = doy_to_date(season_starts[2:5]))

season_text_date <- data.frame(
  date  = doy_to_date(c(30, 106, 198, 290, 350)),
  label = c("Summer", "Autumn", "Winter", "Spring", "Summer")
)

rug_doy_date <- dat |>
  transmute(date = doy_to_date(day_of_year))

# Observed monthly means ± 95% CI, placed at mid-month
monthly_obs <- dat |>
  mutate(month = lubridate::month(timeStamp)) |>
  group_by(month) |>
  summarise(
    mean_obs = mean(Bm_Ant_ABZ, na.rm = TRUE),
    se_obs   = sd(Bm_Ant_ABZ,   na.rm = TRUE) / sqrt(sum(!is.na(Bm_Ant_ABZ))),
    .groups  = "drop"
  ) |>
  mutate(date = as.Date(sprintf("2021-%02d-15", month)))

ggplot(doy_dates, aes(date, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_pointrange(
    data    = monthly_obs,
    mapping = aes(x = date, y = mean_obs,
                  ymin = mean_obs - 1.96 * se_obs,
                  ymax = mean_obs + 1.96 * se_obs),
    colour  = "firebrick", size = 0.4,
    inherit.aes = FALSE
  ) +
  geom_rug(data = rug_doy_date, aes(x = date), sides = "b", alpha = 0.05,
           inherit.aes = FALSE) +
  geom_vline(data = season_vlines_date, aes(xintercept = date),
             linetype = "dashed", colour = "steelblue", linewidth = 0.5) +
  geom_text(data = season_text_date, aes(x = date, label = label),
            y = Inf, vjust = 1.4, size = 3, colour = "steelblue",
            inherit.aes = FALSE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    x        = NULL,
    y        = "Bm_Ant_ABZ",
    title    = "Seasonal marginal prediction from m1 with observed monthly means",
    subtitle = "Line/ribbon: model prediction ± 95% CI  |  Points: observed monthly mean ± 95% CI"
  )

# --- Pairwise marginal prediction heatmaps -----------------------------------
# For each pair of predictors, a 50×50 grid is evaluated; all other predictors
# are held at their means. Fill is sqrt-scaled to reduce dominance of the
# seasonal peak.

make_bivariate <- function(var1, var2, n = 50) {
  v1_seq <- seq(min(dat[[var1]], na.rm = TRUE), max(dat[[var1]], na.rm = TRUE), length.out = n)
  v2_seq <- seq(min(dat[[var2]], na.rm = TRUE), max(dat[[var2]], na.rm = TRUE), length.out = n)
  
  grid <- expand.grid(v1 = v1_seq, v2 = v2_seq)
  names(grid) <- c(var1, var2)
  for (v in setdiff(vars, c(var1, var2))) grid[[v]] <- means[[v]]
  
  pr <- predict(m1, newdata = grid, type = "link", se.fit = FALSE)
  grid$fit  <- exp(pr)
  grid$var1 <- var1
  grid$var2 <- var2
  grid
}

pairs   <- combn(vars, 2, simplify = FALSE)
biv_df  <- lapply(pairs, function(p) make_bivariate(p[1], p[2])) |> bind_rows()

short_labels <- c(
  hour        = "Hour",
  day_of_year = "Day of year",
  solar_lux   = "Solar lux",
  lunar_lux   = "Lunar lux",
  ice         = "Ice"
)

# Arrange into a lower-triangle scatterplot-matrix layout (5×5 grid)
# Lower triangle = heatmap, diagonal = variable label, upper triangle = empty
n_vars <- length(vars)

grid_plots <- vector("list", n_vars ^ 2)

for (col in seq_len(n_vars)) {
  for (row in seq_len(n_vars)) {
    idx <- (row - 1) * n_vars + col
    
    if (row == col) {
      grid_plots[[idx]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = short_labels[vars[row]], size = 3.5, fontface = "bold") +
        theme_void()
      
    } else if (row > col) {
      v1 <- vars[col]   # x-axis variable (column position)
      v2 <- vars[row]   # y-axis variable (row position)
      df <- biv_df |> filter(var1 == v1, var2 == v2)
      
      p <- ggplot(df, aes(x = .data[[v1]], y = .data[[v2]], fill = fit)) +
        geom_raster(interpolate = TRUE) +
        geom_contour(aes(z = fit), colour = "white", alpha = 0.5, linewidth = 0.3) +
        scale_fill_viridis_c(option = "magma", trans = "sqrt", guide = "none") +
        theme_bw(base_size = 7) +
        theme(axis.title  = element_blank(),
              plot.margin = margin(2, 2, 2, 2))
      
      # Suppress tick labels except on the outer edges
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

wrap_plots(grid_plots, ncol = n_vars) +
  plot_annotation(
    title    = "Pairwise marginal predictions from m1",
    subtitle = "Lower triangle: predicted Bm_Ant_ABZ (sqrt-scaled fill) | others held at their means"
  )

# --- Seasonal mean predictions tile plot -------------------------------------
# Each cell: focal predictor at its seasonal mean, day_of_year at the season
# midpoint, all other predictors at their overall mean.
# Rows = seasons, columns = predictors.

seasonal_means <- dat |>
  mutate(season = case_when(
    day_of_year >= 335 | day_of_year < 60  ~ "Summer",
    day_of_year >= 60  & day_of_year < 152 ~ "Autumn",
    day_of_year >= 152 & day_of_year < 244 ~ "Winter",
    day_of_year >= 244 & day_of_year < 335 ~ "Spring"
  )) |>
  group_by(season) |>
  summarise(across(c(hour, day_of_year, solar_lux, lunar_lux, ice),
                   mean, na.rm = TRUE),
            .groups = "drop") |>
  mutate(day_of_year = case_when(
    season == "Summer" ~ 15,
    season == "Autumn" ~ 106,
    season == "Winter" ~ 198,
    season == "Spring" ~ 289
  ))

focal_vars    <- c("hour", "solar_lux", "lunar_lux", "ice")
focal_labels  <- c(
  hour      = "Hour of day",
  solar_lux = "Solar illuminance (lux)",
  lunar_lux = "Lunar illuminance (lux)",
  ice       = "Ice concentration"
)
season_levels <- c("Summer", "Autumn", "Winter", "Spring")
overall_means <- dat |>
  summarise(across(c(hour, solar_lux, lunar_lux, ice), mean, na.rm = TRUE))

tile_df <- crossing(season = season_levels, focal_var = focal_vars) |>
  pmap_dfr(function(season, focal_var) {
    sm             <- seasonal_means |> filter(season == !!season)
    grid           <- overall_means
    grid$day_of_year  <- sm$day_of_year
    grid[[focal_var]] <- sm[[focal_var]]
    
    pr <- predict(mdensity_int, newdata = grid, type = "link", se.fit = FALSE)
    tibble(season = season, focal_var = focal_var, fit = exp(as.numeric(pr)))
  }) |>
  mutate(
    season    = factor(season, levels = season_levels),
    focal_var = factor(focal_var, levels = focal_vars,
                       labels = unname(focal_labels))
  )

ggplot(tile_df, aes(x = focal_var, y = season, fill = fit)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(fit, 1)), size = 3.5) +
  scale_fill_viridis_c( trans = "log",
                       name = "Krill density")#Predicted\nBm_Ant_ABZ") +
  labs(
    x        = NULL,
    y        = NULL,
    title    = "Krill density",#Predicted Bm_Ant_ABZ at seasonal means",
    subtitle = "Each cell: focal predictor at its seasonal mean; day_of_year at season midpoint;\nall other predictors at their overall mean"
  ) +
  theme_bw() +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

##interaction term
m1_ti <- gam(
  Bm_Ant_ABZ ~ s(hour, k = 15, bs = "cc") +
    s(day_of_year, k = 15, bs = "cc") +
    s(solar_lux, k = 6, bs = "cr") +
    s(lunar_lux, k = 6, bs = "cr") +
    s(ice, k = 6, bs = "cr") +
    ti(hour, solar_lux, bs = c("cc", "cr"), k = c(8, 6)),
  family = gaussian(link = "log"),
  data   = dat,
  method = "ML"
)

summary(m1_ti)

##################################################
#KRILL DENSITY
mdensity_int=update(mdensity, . ~ . +ti(hour, solar_lux, bs = c("cc", "cr"), k = c(8, 6)))
summary(mdensity_int)

#KRILL SWARMS(?) CENTRE OF MASS
mcom_int<-update(mcom, . ~ . +ti(hour, solar_lux, bs = c("cc", "cr"), k = c(8, 6)))
summary(mcom_int)

################################################
library(tidyverse)

# 2D marginal prediction grid for a pair of variables
make_bivariate <- function(var1, var2, n = 50,modObj) {
  v1_seq <- seq(min(dat[[var1]], na.rm = TRUE), max(dat[[var1]], na.rm = TRUE), length.out = n)
  v2_seq <- seq(min(dat[[var2]], na.rm = TRUE), max(dat[[var2]], na.rm = TRUE), length.out = n)

  grid <- expand.grid(v1 = v1_seq, v2 = v2_seq)
  names(grid) <- c(var1, var2)

  # Fill remaining predictors at their means
  for (v in setdiff(vars, c(var1, var2))) grid[[v]] <- means[[v]]

  pr <- predict(modObj, newdata = grid, type = "link", se.fit = FALSE)

  grid$fit <- exp(pr)
  grid$var1 <- var1
  grid$var2 <- var2
  grid
}

vars <- c("hour", "day_of_year", "solar_lux", "lunar_lux", "ice")
pairs <- combn(vars, 2, simplify = FALSE)

biv_df <- lapply(pairs, function(p) make_bivariate(p[1], p[2],modObj=mcom_int)) |> bind_rows()

cat("Pairs computed:", length(pairs), "\n")
cat("Total rows:", nrow(biv_df), "\n")

library(patchwork)

n_vars <- length(vars)

# Build 5×5 grid: lower triangle = heatmap, diagonal = label, upper = spacer
grid_plots <- vector("list", n_vars ^ 2)

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
        scale_fill_viridis_c(option = "magma", trans = "sqrt", guide = "none") +
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

wrap_plots(grid_plots, ncol = n_vars) +
  plot_annotation(
    title    = "Pairwise marginal predictions from krill centre of mass",
    subtitle = "Lower triangle: Krill centre of mass (sqrt-scaled fill) | others held at their means"
  )

##################################################
#MVN GAM:
dat2<-dat
dat2$rho_min[is.nan(dat2$rho_min)]=0
summary(dat2$rho_min)
summary(dat2$day_of_year)
summary(dat2$Bm_D)
dat2$Bm_D[is.na(dat2$Bm_D)]=0
m_mvn<- gam(list(Bm_D ~ s(day_of_year, k = 15, 
                          bs = "cc")+s(solar_lux, k = 6, bs = "cr")+ 
                   s(ice, k = 6, bs = "cr") + s(lunar_lux, 
                                                k = 6, bs = "cr"),
                 rho_min ~ s(day_of_year, k = 15, 
                             bs = "cc")+s(solar_lux, k = 6, bs = "cr")+
                    s(ice, k = 6, bs = "cr") + s(lunar_lux, 
                                                 k = 6, bs = "cr")),
                 data=dat2,
                 method='ML',
            family = mvn(d = 2),
            select=TRUE)
summary(m_mvn)

#check if bivariate GAMs are 'better' than two univariate GAMs:
m_Bm_D <- gam(Bm_D ~ s(day_of_year) + s(solar_lux),
              data = dat2, method = "REML")

m_rhomin <- gam(rho_min ~ s(day_of_year) + s(solar_lux),
                data = dat2, method = "REML")

# AIC of bivariate model vs sum of two univariate models
# (valid comparison: same data, same n, comparable log-likelihoods)
aic_mvn  <- AIC(m_mvn)
aic_sep  <- AIC(m_Bm_D) + AIC(m_rhomin)

# Degrees of freedom: mvn has 3 extra params (Cholesky R.1, R.2, R.3)
df_mvn  <- attr(logLik(m_mvn),  "df")
df_sep  <- attr(logLik(m_Bm_D), "df") + attr(logLik(m_rhomin), "df")

data.frame(
  model  = c("bivariate mvn", "two univariate GAMs"),
  AIC    = round(c(aic_mvn, aic_sep), 2),
  df     = c(df_mvn, df_sep),
  delta  = round(c(aic_mvn, aic_sep) - min(c(aic_mvn, aic_sep)), 2)
)

m_mvn_ml <- gam(
  list(Bm_D ~ s(day_of_year) + s(solar_lux),
       rho_min ~ s(day_of_year) + s(solar_lux)),
  family = mvn(d = 2), data = dat2, method = "ML"
)

m_Bm_D_ml   <- gam(Bm_D    ~ s(day_of_year) + s(solar_lux), data = dat2, method = "ML")
m_rhomin_ml <- gam(rho_min ~ s(day_of_year) + s(solar_lux), data = dat2, method = "ML")

ll_mvn <- as.numeric(logLik(m_mvn_ml))
ll_sep <- as.numeric(logLik(m_Bm_D_ml)) + as.numeric(logLik(m_rhomin_ml))

df_mvn <- attr(logLik(m_mvn_ml),  "df")
df_sep <- attr(logLik(m_Bm_D_ml), "df") + attr(logLik(m_rhomin_ml), "df")

aic_mvn <- -2 * ll_mvn + 2 * df_mvn
aic_sep <- -2 * ll_sep + 2 * df_sep

data.frame(
  model      = c("bivariate mvn", "two univariate GAMs"),
  logLik     = round(c(ll_mvn, ll_sep), 2),
  df         = round(c(df_mvn, df_sep), 2),
  AIC        = round(c(aic_mvn, aic_sep), 2),
  delta_AIC  = round(c(aic_mvn, aic_sep) - min(c(aic_mvn, aic_sep)), 2)
)


            
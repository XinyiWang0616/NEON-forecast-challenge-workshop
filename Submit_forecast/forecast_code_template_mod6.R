## install.packages('remotes')
## install.packages('tidyverse') # collection of R packages for data manipulation, analysis, and visualisation
## install.packages('lubridate') # working with dates and times
## remotes::install_github('eco4cast/neon4cast') # package from NEON4cast challenge organisers to assist with forecast building and submission
    
# ------ Load packages -----
library(tidyverse)
library(lubridate)
#--------------------------#
    
# Change this for your model ID
# Include the word "example" in my_model_id for a test submission
# Don't include the word "example" in my_model_id for a forecast that you have registered (see neon4cast.org for the registration form)
my_model_id <- 'example_ID'
    
# --Model description--- #
# Drivers include: air temperature, eastern wind, and solar radiation
# Target: Water temperature/oxygen
    
# -- Uncertainty representation -- #
# Uncertainty description: Process uncertainty, Drivers uncertainty, and Parameter uncertainty are included in the forecasting.
#------- Read data --------
# read in the targets data
targets <- read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz")
    
# read in the sites data
aquatic_sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ci/refs/heads/main/neon4cast_field_site_metadata.csv") |>
dplyr::filter(aquatics == 1)
    
focal_sites <- aquatic_sites |> 
  filter(field_site_subtype == 'Lake') |> 
  pull(field_site_id)
    
# Filter the targets
targets <- targets %>%
  filter(site_id %in% focal_sites,
          variable %in% c('temperature', 'oxygen'))
#--------------------------#
    
    
    
# ------ Weather data ------
met_variables <- c("air_temperature", 
                    "eastward_wind",
                    "surface_downwelling_shortwave_flux_in_air")
    
# Past stacked weather -----
weather_past_s3 <- neon4cast::noaa_stage3()
    
weather_past <- weather_past_s3  |> 
  dplyr::filter(site_id %in% focal_sites,
                datetime >= ymd('2017-01-01'),
                variable %in% met_variables) |> 
  dplyr::collect()
    
# aggregate the past to mean values
weather_past_daily <- weather_past |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction)
    
# Future weather forecast --------
# New forecast only available at 5am UTC the next day
forecast_date <- Sys.Date() 
noaa_date <- forecast_date - days(1)
    
weather_future_s3 <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))
    
weather_future <- weather_future_s3 |> 
  dplyr::filter(datetime >= forecast_date,
                site_id %in% focal_sites,
                variable %in% met_variables) |> 
  collect()
    
weather_future_daily <- weather_future |> 
  mutate(datetime = as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, site_id, parameter, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) |> 
  select(any_of(c('datetime', 'site_id', met_variables, 'parameter')))
    

# ----- Fit model & generate forecast----
    
# Generate a dataframe to fit the model to 
targets_lm <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, 
            by = c("datetime","site_id"))
    
# Loop through each site to fit the model
forecast_df <- NULL
n_drivers <- n_distinct(weather_future_daily$parameter) # 31
n_rep <- 10
n_members <- n_drivers * n_rep

for(i in 1:length(focal_sites)) {  
      
  curr_site <- focal_sites[i]
      
  # historical observations: training
  site_target <- targets_lm |>
    filter(site_id == curr_site) |>
    arrange(datetime) |>
    mutate(
      temp_lag = lag(temperature, 1),
      air_temperature_lag = lag(air_temperature, 1),
      shortwave_lag = lag(surface_downwelling_shortwave_flux_in_air, 1),
      air_temp_7days = as.numeric(stats::filter(air_temperature, rep(1/7, 7), sides = 1))) |>
    drop_na(temperature, temp_lag, air_temperature, eastward_wind, surface_downwelling_shortwave_flux_in_air, air_temperature_lag, shortwave_lag, air_temp_7days)
      
  # future weather forecasting: predicting
  noaa_future_site <- weather_future_daily |> 
    filter(site_id == curr_site) |>
    arrange(parameter, datetime) |>
    group_by(parameter) |>
    mutate(
      air_temperature_lag = lag(air_temperature, 1),
      shortwave_lag = lag(surface_downwelling_shortwave_flux_in_air, 1)) |>
    ungroup()
  last_hist <- site_target |>
    arrange(datetime) |>
    slice_tail(n = 1)
  noaa_future_site <- noaa_future_site |>
    group_by(parameter) |>
    mutate(air_temperature_lag = if_else(row_number() == 1,
                                         last_hist$air_temperature,
                                         air_temperature_lag),
           shortwave_lag = if_else(row_number() == 1, last_hist$surface_downwelling_shortwave_flux_in_air,
                                   shortwave_lag)) |>
    ungroup()
      
  #Fit linear model based on past data: water temperature = m * air temperature + b
  #you will need to change the variable on the left side of the ~ if you are forecasting oxygen or chla
  fit <- lm(temperature ~ temp_lag + air_temperature +
              eastward_wind +
              surface_downwelling_shortwave_flux_in_air +
              I(air_temperature^2) + 
              air_temperature_lag + 
              shortwave_lag + air_temp_7days, data=site_target)
      
  fit_summary <- summary(fit)
  coeffs <- round(fit$coefficients, 2)
  params_se <- fit_summary$coefficients[,2]
  
  param_df <- tibble(member_id = 1:n_members,
                     driver_parameter = rep(sort(unique(noaa_future_site$parameter)), each = n_rep),
                     rep_id = rep(1:n_rep, times = n_drivers),
                      beta1 = rnorm(n_members, coeffs[1], params_se[1]),
                          beta2 = rnorm(n_members, coeffs[2], params_se[2]),
                          beta3 = rnorm(n_members, coeffs[3], params_se[3]),
                          beta4 = rnorm(n_members, coeffs[4], params_se[4]),
                          beta5 = rnorm(n_members, coeffs[5], params_se[5]),
                          beta6 = rnorm(n_members, coeffs[6], params_se[6]),
                          beta7 = rnorm(n_members, coeffs[7], params_se[7]),
                          beta8 = rnorm(n_members, coeffs[8], params_se[8]),
                          beta9 = rnorm(n_members, coeffs[9], params_se[9]))
      
  mod <- predict(fit, newdata = site_target)
  # mod <- c(NA, mod)
  residuals <- mod - site_target$temperature
  sigma <- sd(residuals, na.rm = TRUE)
      
  #------------------ initial condition uncertainty---------------
  last_temp <- site_target |>
    arrange(datetime) |>
    slice_tail(n = 1) |>
    pull(temperature)
  sigma_ic <- sigma
  
  init_df <- tibble(member_id = 1:n_members,
                    init_temp = rnorm(n_members, mean = last_temp, sd = sigma_ic))
  forecast_list <- list()
  
  # noaa_future_site <- noaa_future_site |>
    # arrange(parameter, datetime)
  
  for (j in 1:n_members){
    current_driver <- param_df$driver_parameter[j]
    
    future_j <- noaa_future_site|>
      filter(parameter == current_driver) |>
      arrange(datetime)
  pred_j <- numeric(nrow(future_j))
  temp_prev <- init_df$init_temp[j]
  
  air_temp_hist <- site_target |>
    arrange(datetime) |>
    slice_tail(n = 6) |>
    pull(air_temperature)
  
  for (t in 1:nrow(future_j)){
    current_air_temp <- future_j$air_temperature[t]
    air_temp_window <- c(air_temp_hist, current_air_temp)
    air_temp_7days <- mean(air_temp_window, na.rm = TRUE) # rolling window
    
    pred <- param_df$beta1[j] + 
      param_df$beta2[j] * temp_prev +
      param_df$beta3[j] * future_j$air_temperature[t] +
      param_df$beta4[j] * future_j$eastward_wind[t] +
      param_df$beta5[j] * future_j$surface_downwelling_shortwave_flux_in_air[t] +
      param_df$beta6[j] * (future_j$air_temperature[t])^2 +
      param_df$beta7[j] * future_j$air_temperature_lag[t] +
      param_df$beta8[j] * future_j$shortwave_lag[t] + 
      param_df$beta9[j] * air_temp_7days
    
    pred_j[t] <- pred + rnorm(1, mean=0, sd=sigma)
    temp_prev <- pred_j[t]
    air_temp_hist <- c(tail(air_temp_window, 5), current_air_temp)
  }
  forecast_list[[j]] <- tibble(datetime = future_j$datetime,
                               site_id = curr_site,
                               parameter = as.character(j),
                               prediction = pred_j,
                               variable = "temperature")
  }
      
  # use linear regression to forecast water temperature for each ensemble member
  # You will need to modify this line of code if you add additional weather variables or change the form of the model
  # The model here needs to match the model used in the lm function above (or what model you used in the fit)
  # put all the relevant information into a tibble that we can bind together
  # forecasting result for single site
  curr_site_df <- bind_rows(forecast_list) #Change this if you are forecasting a different variable
      
  forecast_df <- dplyr::bind_rows(forecast_df, curr_site_df)
  message(curr_site, 'forecast run')
      }
    
    
#--------------------------#
    
    
#---- Covert to EFI standard ----
    
# Make forecast fit the EFI standards
forecast_df_EFI <- forecast_df %>%
  filter(datetime > forecast_date) %>%
  mutate(model_id = my_model_id,
          reference_datetime = forecast_date,
          family = 'ensemble',
          duration = 'P1D',
          parameter = as.character(parameter),
          project_id = 'neon4cast') %>%
  select(datetime, reference_datetime, duration, site_id, family, parameter, variable, prediction, model_id, project_id)
#---------------------------#
    
    
    
    
# ----- Submit forecast -----
# Write the forecast to file
theme <- 'aquatics'
date <- forecast_df_EFI$reference_datetime[1]
forecast_name <- paste0(forecast_df_EFI$model_id[1], ".csv")
forecast_file <- paste(theme, date, forecast_name, sep = '-')
    
write_csv(forecast_df_EFI, forecast_file)
    
neon4cast::forecast_output_validator(forecast_file)
    
    
neon4cast::submit(forecast_file =  forecast_file, ask = FALSE) # if ask = T (default), it will produce a pop-up box asking if you want to submit
    
#--------------------------#
Sys.setlocale("LC_TIME", "English")
forecast_df_EFI |> 
  ggplot(aes(x=datetime, y=prediction, group = parameter)) +
  geom_line() +
  facet_wrap(~site_id) +
  labs(title = paste0('Forecast generated for ', forecast_df_EFI$variable[1], ' on ', forecast_df_EFI$reference_datetime[1]))
    
plot_file_name <- paste0("Submit_forecast/figures/", forecast_df_EFI$variable[1], '-', forecast_df_EFI$reference_datetime[1], ".png")
ggsave(plot_file_name, width = 12, height = 8)
    
# plot needs to be submitted to Canvas
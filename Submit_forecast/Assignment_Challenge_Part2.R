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
         variable == 'temperature')
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

#--------------------------#



# ----- Fit model & generate forecast----

# Generate a dataframe to fit the model to 
targets_lm <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, by = c("datetime","site_id"))

# Loop through each site to fit the model
forecast_df <- NULL

for(i in 1:length(focal_sites)) {  
  
  curr_site <- focal_sites[i]
  
  # historical observations: training
  site_target <- targets_lm |>
    filter(site_id == curr_site) |>
    arrange(datetime) |>
    mutate(temp_lag = lag(temperature, 1)) |>
    drop_na(temperature, air_temperature, eastward_wind, surface_downwelling_shortwave_flux_in_air)
  
  # future weather forecasting: predicting
  noaa_future_site <- weather_future_daily |> 
    filter(site_id == curr_site) |>
    arrange(parameter, datetime)
  
  #Fit linear model based on past data: water temperature = m * air temperature + b
  #you will need to change the variable on the left side of the ~ if you are forecasting oxygen or chla
  # fit <- lm(site_target$temperature ~ site_target$air_temperature)
  fit <- lm(temperature ~ temp_lag + air_temperature + I(air_temperature^2) +
              eastward_wind + surface_downwelling_shortwave_flux_in_air, data=site_target)
  
  fit_summary <- summary(fit)
  coeffs <- coef(fit)
  params_se <- fit_summary$coefficients[,"Std. Error"]
  unique_params <- sort(unique(noaa_future_site$parameter))
  n_members <- length(unique_params)
  param_df <- tibble(parameter = unique_params,
                     beta1 = rnorm(n_members, coeffs[1], params_se[1]),
                     beta2 = rnorm(n_members, coeffs[2], params_se[2]),
                     beta3 = rnorm(n_members, coeffs[3], params_se[3]),
                     beta4 = rnorm(n_members, coeffs[4], params_se[4]),
                     beta5 = rnorm(n_members, coeffs[5], params_se[5]),
                     beta6 = rnorm(n_members, coeffs[6], params_se[6]))
  
  mod <- predict(fit, newdata = site_target)
  residuals <- site_target$temperature - mod
  sigma <- sd(residuals, na.rm = TRUE)
  
  last_temp <- site_target |>
    arrange(datetime) |>
    slice_tail(n=1) |>
    pull(temperature)
  
  sigma_ic <- sigma
  init_df <- tibble(parameter = unique_params, init_temp = rnorm(n_members, mean=last_temp, sd = sigma_ic))
  
  forecast_list <- list()
  for (j in seq_along(unique_params)){
    current_param <- unique_params[j]
      
    # future weather series
    future_j <- noaa_future_site |>
      filter(parameter==current_param)|>
      arrange(datetime)
    
    # if (nrow(future_j)==0){next}
    
    beta_j <- param_df |>
      filter(parameter == current_param)
    
    temp_prev <- init_df |>
      filter(parameter == current_param) |>
      pull(init_temp)
    
    pred_j <- numeric(nrow(future_j))
    
    for (t in 1:nrow(future_j)){
      predictions <- beta_j$beta1 + beta_j$beta2 * temp_prev + 
        beta_j$beta3 * future_j$air_temperature[t] +
        beta_j$beta4 * (future_j$air_temperature[t]^2) +
        beta_j$beta5 * future_j$eastward_wind[t] +
        beta_j$beta6 * future_j$surface_downwelling_shortwave_flux_in_air[t]
      
      pred_j[t] <- predictions + rnorm(1, mean=0, sd=sigma)
      temp_prev <- pred_j[t]
    }
    forecast_list[[j]] <- tibble(datetime=future_j$datetime, site_id = curr_site, parameter = as.character(current_param),
                                 prediction=pred_j, variable = "temperature")
  }
  
  # use linear regression to forecast water temperature for each ensemble member
  # You will need to modify this line of code if you add additional weather variables or change the form of the model
  # The model here needs to match the model used in the lm function above (or what model you used in the fit)
  
  curr_site_df <- dplyr::bind_rows(forecast_list)
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

plot_file_name <- paste0("Submit_forecast/", forecast_df_EFI$variable[1], '-', forecast_df_EFI$reference_datetime[1], ".png")
ggsave(plot_file_name, width = 12, height = 8)

# plot needs to be submitted to Canvas
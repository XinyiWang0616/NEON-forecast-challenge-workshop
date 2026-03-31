library(tidyverse)
library(scoringRules)
library(duckdbfs)

#Get this link from https://radiantearth.github.io/stac-browser/#/external/raw.githubusercontent.com/eco4cast/neon4cast-catalog/main/forecasts/Aquatics/Daily_Water_temperature/collection.json
# You will need to modify the link to remove the anonymous@ and the ?endpoint_override... part
all_results <- open_dataset("s3://bio230014-bucket01/challenges/forecasts/bundled-parquet/project_id=neon4cast/duration=P1D/variable=temperature", 
                            s3_endpoint = "sdsc.osn.xsede.org", 
                            anonymous = TRUE)

# Get today's forecast from two models for a particular site and reference_datetime
my_forecasts <- all_results |> 
  filter(reference_datetime == as_datetime("2026-03-24 00:00:00"),
         site_id == "BARC",
         model_id %in% c("climatology", "persistenceRW")) |> 
  collect()   # flareGLM

# Plot a parametric forecast
my_forecasts |> 
  filter(model_id == "climatology") |> 
  pivot_wider(names_from = parameter, values_from = prediction) |> 
  mutate(lower_2.5 = mu - 1.96*sigma,
         upper_97.5 = mu + 1.96*sigma) |> 
  ggplot(aes(x = datetime)) +
  geom_ribbon(aes(ymin = lower_2.5, ymax = upper_97.5), fill = "lightblue") +
  geom_line(aes(y = mu))

# Plot an ensemble forecast
my_forecasts |> 
  filter(model_id == "persistenceRW") |> 
  ggplot(aes(x = datetime, y = prediction, group = parameter)) +
  geom_line()


# Get a past forecast from two models for a particular site and reference_datetime
# We are getting a past forecast so we can evaluate it with data
my_forecasts <- all_results |> 
  filter(reference_datetime == as_datetime("2025-03-24 00:00:00"),
         site_id == "BARC",
         model_id %in% c("climatology", "persistenceRW")) |> 
  collect()

#Lets look at a single forecast
#This one site_id, one model_id, one variable, one reference_datetime, one datetime
# which is the unit of a single forecast
# my_forecast is already filtered for a single site and reference_datetime
# so we need to filter for a model_id and a datetime

single_forecast_ensemble <- my_forecasts |> 
  filter(model_id == "persistenceRW",
         datetime == as_datetime("2025-03-27 00:00:00"))

# We want to match the single forecast to an observation
# Get the target data and read in
# Target data is here https://projects.ecoforecast.org/neon4cast-ci/targets.html
# Aquatics tab
url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz"
targets <- read_csv(url)

# Filter to only the datetime, variable, and site_id that match our single forecast
targets_single_forecast <- targets |> 
  filter(datetime == as_datetime("2025-03-27 00:00:00"),
         variable == "temperature",
         site_id == "BARC")

#Calculate CRPS
crps_sample(y = targets_single_forecast$observation,
            dat = single_forecast_ensemble$prediction)

#Calculate Log Score (probability of observation given the forecast distribution)
logs_sample(y = targets_single_forecast$observation,
            dat = single_forecast_ensemble$prediction)

#Get a single parameteric forecast
single_forecast_parametric <- my_forecasts |> 
  filter(model_id == "climatology",
         datetime == as_datetime("2025-03-27 00:00:00")) |> 
  pivot_wider(names_from = parameter, values_from = prediction) 

#Calculate CRPS for the parametric distribution
crps_norm(y = targets_single_forecast$observation, 
          mean = single_forecast_parametric$mu,
          sd = single_forecast_parametric$sigma)

# Lets score multiple forecasts!

multi_forecast_ensemble <- my_forecasts |> 
  filter(model_id == "persistenceRW") |> 
  left_join(targets, by = c("datetime", "site_id", "variable")) |> 
  select(model_id, datetime, prediction, observation, parameter)

multi_forecast_parametric <- my_forecasts |> 
  filter(model_id == "climatology") |> 
  left_join(targets, by = c("datetime", "site_id", "variable")) |> 
  select(model_id, datetime, prediction, observation, parameter) |> 
  pivot_wider(names_from = parameter, values_from = prediction)

normal_crps <- multi_forecast_parametric |> 
  group_by(datetime, model_id) |> 
  summarize(crps = crps_norm(observation, mean = mu, sd = sigma))


# Scoring multiple forecaasts using crps_sample in a group_by summarize required
# making a small function to get a single observation
crps_tidy_ensemble <- function(prediction, observation) {
  obs <- observation[1] #all the observations are the same for each row, but we only need one
  crps_sample(obs, prediction)
}

#Use the function
ensemble_crps <- multi_forecast_ensemble |> 
  group_by(datetime, model_id) |> 
  summarize(crps = crps_tidy_ensemble(prediction, observation)) 

#Plot the histogram of the scores

bind_rows(normal_crps, ensemble_crps) |> 
  ggplot(aes(y = crps, fill = model_id)) +
  geom_histogram()


# load gaia
library(gaia)

# Path to the output folder. Change it to your desired location
output_dir <- 'impact'

# MRI-ESM2-0: Food-95----

pr_projection_file <- "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-95_pr_global_monthly_2015_2100.nc"
tas_projection_file <- "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-95_tas_global_monthly_2015_2100.nc"


gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = pr_projection_file,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = tas_projection_file,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/MRI-ESM2-0/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'MRI-ESM2-0',
  climate_scenario = 'Food-95',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = FALSE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

## step-by-step ----
gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = NULL,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = NULL,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/MRI-ESM2-0/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'MRI-ESM2-0',
  climate_scenario = 'Food-95',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = TRUE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

out_yield_shock <- yield_shock_projection(use_default_coeff = FALSE,
                                          climate_model = 'MRI-ESM2-0',
                                          climate_scenario = 'Food-95',
                                          base_year = 2015,
                                          start_year = 2015,
                                          end_year = 2100,
                                          smooth_window = 20,
                                          diagnostics = FALSE,
                                          output_dir = output_dir)


gcam_apg <- gcam_agprodchange(data = out_yield_shock,
                              climate_model = 'MRI-ESM2-0',
                              climate_scenario = 'Food-95',
                              member = 'r1i1p1f1',
                              bias_adj = 'w5e5',
                              cfe = 'no-cfe',
                              gcam_version = 'gcam7',
                              diagnostics = FALSE,
                              output_dir = output_dir)


write.csv(gcam_apg, file = "agyield_impact_mri-esm2-0_r1i1p1f1_w5e5_food-95.csv")





# CanESM5: Food-95----

pr_projection_file <- "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-95_pr_global_monthly_2015_2100.nc"
tas_projection_file <- "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-95_tas_global_monthly_2015_2100.nc"


gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = pr_projection_file,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = tas_projection_file,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/CanESM5/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'CanESM5',
  climate_scenario = 'Food-95',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = FALSE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

## step-by-step ----
gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = NULL,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = NULL,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/CanESM5/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'CanESM5',
  climate_scenario = 'Food-95',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = TRUE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

out_yield_shock <- yield_shock_projection(use_default_coeff = FALSE,
                                          climate_model = 'CanESM5',
                                          climate_scenario = 'Food-95',
                                          base_year = 2015,
                                          start_year = 2015,
                                          end_year = 2100,
                                          smooth_window = 20,
                                          diagnostics = FALSE,
                                          output_dir = output_dir)


gcam_apg <- gcam_agprodchange(data = out_yield_shock,
                              climate_model = 'CanESM5',
                              climate_scenario = 'Food-95',
                              member = 'r1i1p1f1',
                              bias_adj = 'w5e5',
                              cfe = 'no-cfe',
                              gcam_version = 'gcam7',
                              diagnostics = FALSE,
                              output_dir = output_dir)


write.csv(gcam_apg, file = "agyield_impact_canesm5_r1i1p1f1_w5e5_food-95.csv")






# CanESM5: Food-ref----

pr_projection_file <- "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-ref_pr_global_monthly_2015_2100.nc"
tas_projection_file <- "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-ref_tas_global_monthly_2015_2100.nc"


gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = pr_projection_file,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = tas_projection_file,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/CanESM5/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'CanESM5',
  climate_scenario = 'Food-ref',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = FALSE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

## step-by-step ----
gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = NULL,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = NULL,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/CanESM5/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'CanESM5',
  climate_scenario = 'Food-ref',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = TRUE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

out_yield_shock <- yield_shock_projection(use_default_coeff = FALSE,
                                          climate_model = 'CanESM5',
                                          climate_scenario = 'Food-ref',
                                          base_year = 2015,
                                          start_year = 2015,
                                          end_year = 2100,
                                          smooth_window = 20,
                                          diagnostics = FALSE,
                                          output_dir = output_dir)


gcam_apg <- gcam_agprodchange(data = out_yield_shock,
                              climate_model = 'CanESM5',
                              climate_scenario = 'Food-ref',
                              member = 'r1i1p1f1',
                              bias_adj = 'w5e5',
                              cfe = 'no-cfe',
                              gcam_version = 'gcam7',
                              diagnostics = FALSE,
                              output_dir = output_dir)


write.csv(gcam_apg, file = "agyield_impact_canesm5_r1i1p1f1_w5e5_food-ref.csv")






# MRI-ESM2-0: Food-ref----

pr_projection_file <- "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-ref_pr_global_monthly_2015_2100.nc"
tas_projection_file <- "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-ref_tas_global_monthly_2015_2100.nc"


gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = pr_projection_file,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = tas_projection_file,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/MRI-ES2-0/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'MRI-ES2-0',
  climate_scenario = 'Food-ref',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = FALSE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

## step-by-step ----
gaia::yield_impact(
  pr_hist_ncdf = NULL,                    # path to historical precipitation NetCDF file (must follow ISIMIP format); only if you wish to use your own historical precipitation observation
  tas_hist_ncdf = NULL,                   # path to historical temperature NetCDF file (must follow ISIMIP format); only if you wish to use your own historical temperature observation
  pr_proj_ncdf = NULL,      # path to future projected precipitation NetCDF file (must follow ISIMIP format)
  tas_proj_ncdf = NULL,    # path to future projected temperature NetCDF file (must follow ISIMIP format)
  timestep = 'monthly',                   # specify the time step of the NetCDF data (monthly or daily)
  climate_hist_dir = 'C:/Model/gaia/gaia_example/climate_hist',
  climate_impact_dir = 'C:/Model/gaia/impact/weighted_climate/MRI-ES2-0/',
  historical_periods = c(1960:2001),      # vector of historical years selected for fitting
  climate_model = 'MRI-ES2-0',
  climate_scenario = 'Food-ref',
  member = 'r1i1p1f1',                    # label of ensemble member name
  bias_adj = 'w5e5',                      # label of climate data for bias adjustment for the global climate model (GCM)
  cfe = 'no-cfe',                         # label of CO2 fertilization effect in the formula (default is no CFE)
  gcam_version = 'gcam7',                 # output is different depending on the GCAM version (gcam6 or gcam7)
  use_default_coeff = FALSE,              # set to TRUE when there is no historical climate data available
  base_year = 2015,                       # GCAM base year
  start_year = 2015,                      # start year of the projected climate data
  end_year = 2100,                        # end year of the projected climate data
  smooth_window = 20,                     # number of years as smoothing window
  co2_hist = NULL,                        # historical annual CO2 concentration. If NULL, will use default value
  co2_proj = NULL,                        # projected annual CO2 concentration. If NULL, will use default value
  crop_select = NULL,                     # set to NULL for the default crops
  diagnostics = TRUE,                     # set to TRUE to output diagnostic plots
  output_dir = output_dir                 # path to the output folder
)

out_yield_shock <- yield_shock_projection(use_default_coeff = FALSE,
                                          climate_model = 'MRI-ES2-0',
                                          climate_scenario = 'Food-ref',
                                          base_year = 2015,
                                          start_year = 2015,
                                          end_year = 2100,
                                          smooth_window = 20,
                                          diagnostics = FALSE,
                                          output_dir = output_dir)


gcam_apg <- gcam_agprodchange(data = out_yield_shock,
                              climate_model = 'MRI-ES2-0',
                              climate_scenario = 'Food-ref',
                              member = 'r1i1p1f1',
                              bias_adj = 'w5e5',
                              cfe = 'no-cfe',
                              gcam_version = 'gcam7',
                              diagnostics = FALSE,
                              output_dir = output_dir)


write.csv(gcam_apg, file = "agyield_impact_mri-esm2-0_r1i1p1f1_w5e5_food-ref.csv")





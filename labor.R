library(HELPS)

SOI_LIST <- SECTOR_ALL
# "OFIB_I", "BANA_I", "PLNT_I", "CITR_I", "TROF_I", "TEMF_I", "TOMA_I", "ONIO_I", "VEGE_I",
# SOI_LIST <- c("ORTS_I",
#               "OPUL_I", "CNUT_I", "OILP_I", "OOIL_I", "COFF_I",
#               "RCOF_I", "COCO_I", "RUBB_I", "TEAS_I", "TOBA_I", "REST_I", "WHEA_R", "RICE_R", "MAIZ_R", "SOYB_R",
#               "BARL_R", "MILL_R", "PMIL_R", "SORG_R", "OCER_R", "POTA_R", "SWPO_R", "YAMS_R", "CASS_R", "BEAN_R",
#               "CHIC_R", "COWP_R", "PIGE_R", "LENT_R", "GROU_R", "SUNF_R", "RAPE_R", "SESA_R", "SUGC_R", "SUGB_R",
#               "COTT_R", "OFIB_R", "BANA_R", "PLNT_R", "CITR_R", "TROF_R", "TEMF_R", "TOMA_R", "ONIO_R", "VEGE_R",
#               "ORTS_R", "OPUL_R", "CNUT_R", "OILP_R", "OOIL_R", "COFF_R", "RCOF_R", "COCO_R", "RUBB_R", "TEAS_R",
#               "TOBA_R", "REST_R", "NONCROP")


# store the climate impact files into a path: in this script, the files are stored in `C:/Model/Claudia T/`. Change accordingly.
# run crop yield impact before running this script. Save crop yield impact in both csv and xml.

# Inputs ----

LHR_Dunne <- function(WBGT, workload = NULL){ # Dunne et al., 2013
  eta = ifelse(WBGT <= 25, 1, 1 - 0.25*((WBGT - 25)^(2/3)))
}

YEAR <- seq(2015, 2100, 1); YEAR # all future years
# YEAR <- seq(2015, 2100, 5); YEAR # GCAM model years

## Food-ref: CanESM5 ----
start_t = Sys.time()
for (s in 1:length(SOI_LIST)){
  ANNUAL_REG <- list()
  SOI <- SOI_LIST[[s]]
  print(SOI)
  for (i in 1:length(YEAR)){
    YEAR_INPUT <- YEAR[i]
    print(YEAR_INPUT)
    esi.mon <- cal_heat_stress(TempRes = "month", SECTOR = SOI, HS = WBGT_ESI, YEAR_INPUT = YEAR_INPUT,
                               "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-ref_hurs_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-ref_tas_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-ref_rsds_global_monthly_2015_2100.nc")

    pwc.mon <- cal_pwc(WBGT = esi.mon,  LHR = LHR_Dunne, workload = "high")
    rm(esi.mon)
    pwc.ann <- monthly_to_annual(input_rack = pwc.mon, SECTOR = SOI)
    rm(pwc.mon)
    # ANNUAL_GRID[[i]] <- pwc.foster.ann
    reg_pwc <- grid_to_region(grid_annual_value = pwc.ann, SECTOR = SOI, rast_boundary = reg_WB_raster)
    rm(pwc.ann)
    ANNUAL_REG[[i]] <- reg_pwc %>% dplyr::mutate(crop = SOI, year = YEAR_INPUT)
  }
  saveRDS(ANNUAL_REG, file = paste0("C:/Model/Claudia T/ANNUAL_WB_Food-ref_CanESM5_Dunne_", SOI,".rds"))
}
end_t = Sys.time()
end_t - start_t








## Food-95: CanESM5 ----
start_t = Sys.time()
for (s in 1:length(SOI_LIST)){
  ANNUAL_REG <- list()
  SOI <- SOI_LIST[[s]]
  print(SOI)
  for (i in 1:length(YEAR)){
    YEAR_INPUT <- YEAR[i]
    print(YEAR_INPUT)
    esi.mon <- cal_heat_stress(TempRes = "month", SECTOR = SOI, HS = WBGT_ESI, YEAR_INPUT = YEAR_INPUT,
                               "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-95_hurs_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-95_tas_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/CanESM5_STITCHES_W5E5v2_Food-95_rsds_global_monthly_2015_2100.nc")

    pwc.mon <- cal_pwc(WBGT = esi.mon,  LHR = LHR_Dunne, workload = "high")
    rm(esi.mon)
    pwc.ann <- monthly_to_annual(input_rack = pwc.mon, SECTOR = SOI)
    rm(pwc.mon)
    # ANNUAL_GRID[[i]] <- pwc.foster.ann
    reg_pwc <- grid_to_region(grid_annual_value = pwc.ann, SECTOR = SOI, rast_boundary = reg_WB_raster)
    rm(pwc.ann)
    ANNUAL_REG[[i]] <- reg_pwc %>% dplyr::mutate(crop = SOI, year = YEAR_INPUT)
  }
  saveRDS(ANNUAL_REG, file = paste0("C:/Model/Claudia T/ANNUAL_WB_Food-95_CanESM5_Dunne_", SOI,".rds"))
}
end_t = Sys.time()
end_t - start_t


## Food-95: MRI ----
start_t = Sys.time()
for (s in 1:length(SOI_LIST)){
  ANNUAL_REG <- list()
  SOI <- SOI_LIST[[s]]
  print(SOI)
  for (i in 1:length(YEAR)){
    YEAR_INPUT <- YEAR[i]
    print(YEAR_INPUT)
    esi.mon <- cal_heat_stress(TempRes = "month", SECTOR = SOI, HS = WBGT_ESI, YEAR_INPUT = YEAR_INPUT,
                               "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-95_hurs_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-95_tas_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-95_rsds_global_monthly_2015_2100.nc")

    pwc.mon <- cal_pwc(WBGT = esi.mon,  LHR = LHR_Dunne, workload = "high")
    rm(esi.mon)
    pwc.ann <- monthly_to_annual(input_rack = pwc.mon, SECTOR = SOI)
    rm(pwc.mon)
    # ANNUAL_GRID[[i]] <- pwc.foster.ann
    reg_pwc <- grid_to_region(grid_annual_value = pwc.ann, SECTOR = SOI, rast_boundary = reg_WB_raster)
    rm(pwc.ann)
    ANNUAL_REG[[i]] <- reg_pwc %>% dplyr::mutate(crop = SOI, year = YEAR_INPUT)
  }
  saveRDS(ANNUAL_REG, file = paste0("C:/Model/Claudia T/ANNUAL_WB_Food-95_MRI_Dunne_", SOI,".rds"))
}
end_t = Sys.time()
end_t - start_t








## Food-ref: MRI ----
start_t = Sys.time()
for (s in 1:length(SOI_LIST)){
  ANNUAL_REG <- list()
  SOI <- SOI_LIST[[s]]
  print(SOI)
  for (i in 1:length(YEAR)){
    YEAR_INPUT <- YEAR[i]
    print(YEAR_INPUT)
    esi.mon <- cal_heat_stress(TempRes = "month", SECTOR = SOI, HS = WBGT_ESI, YEAR_INPUT = YEAR_INPUT,
                               "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-ref_hurs_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-ref_tas_global_monthly_2015_2100.nc",
                               "C:/Model/Claudia T/MRI-ESM2-0_STITCHES_W5E5v2_Food-ref_rsds_global_monthly_2015_2100.nc")

    pwc.mon <- cal_pwc(WBGT = esi.mon,  LHR = LHR_Dunne, workload = "high")
    rm(esi.mon)
    pwc.ann <- monthly_to_annual(input_rack = pwc.mon, SECTOR = SOI)
    rm(pwc.mon)
    # ANNUAL_GRID[[i]] <- pwc.foster.ann
    reg_pwc <- grid_to_region(grid_annual_value = pwc.ann, SECTOR = SOI, rast_boundary = reg_WB_raster)
    rm(pwc.ann)
    ANNUAL_REG[[i]] <- reg_pwc %>% dplyr::mutate(crop = SOI, year = YEAR_INPUT)
  }
  saveRDS(ANNUAL_REG, file = paste0("C:/Model/Claudia T/ANNUAL_WB_Food-ref_MRI_Dunne_", SOI,".rds"))
}
end_t = Sys.time()
end_t - start_t



# XMLs ----
library(gcamdata)
# translate HELPS outputs to GCAM input
# read in required csv inputs
basin_to_country_mapping <- read.csv("inst/extdata/water/basin_to_country_mapping.csv", skip = 7, header = T)
# FAO_ag_items_PRODSTAT <- read.csv("inst/extdata/aglu/FAO/FAO_ag_items_PRODSTAT.csv", skip = 10, header = T)

L2082.AgCoef_laborcapital_bio_irr_mgmt <- load_from_cache("L2082.AgCoef_laborcapital_bio_irr_mgmt")[[1]]
L2082.AgCoef_laborcapital_ag_irr_mgmt <- load_from_cache("L2082.AgCoef_laborcapital_ag_irr_mgmt")[[1]]
# update labor IO under heat stress -------------------------------------

# get all subsectors ----
L2082.AgCoef_laborcapital_bio_irr_mgmt %>%
  filter(year == 2015, minicam.energy.input == "Labor_Ag") %>%
  separate(AgSupplySubsector, into = c("GCAM_Subsector", "waterbasin"), sep = "_") %>%
  select(AgSupplySector, GCAM_Subsector) %>% unique() -> biomass_sec

L2082.AgCoef_laborcapital_ag_irr_mgmt %>%
  filter(year == 2015, minicam.energy.input == "Labor_Ag") %>%
  separate(AgSupplySubsector, into = c("GCAM_Subsector", "waterbasin"), sep = "_") %>%
  select(AgSupplySector, GCAM_Subsector) %>% unique() -> crop_sec

biomass_sec %>% bind_rows(crop_sec) -> GCAM.sec; GCAM.sec

crop_mapping <- read.csv("C:/Model/HS_package/GCAM_SPAM_crop_mapping.csv") %>%
  select(SPAM, GCAM_Subsector) %>%
  mutate(SPAM = trimws(SPAM))

input.dir <- "C:/Model/Claudia T/"
output.dir <- "../../"


## define scenarios and GCM ----

### MRI-ESM2-0: Food-95 ----
GCM <- "MRI"
SCE <- "Food-95"

eta_list <- list.files(input.dir, pattern = paste0("ANNUAL_WB_",SCE,"_",GCM,"_Dunne"), full.names = FALSE); eta_list

eta_list <- eta_list[!grepl("NONCROP", eta_list)] # no NONCROP mapping for GCAM crops

yld_HS <- read.csv("C:/Model/gaia/agyield_impact_mri-esm2-0_r1i1p1f1_w5e5_food-95.csv") %>%
  select(-X)

### MRI-ESM2-0: Food-ref ----
GCM <- "MRI"
SCE <- "Food-ref"

eta_list <- list.files(input.dir, pattern = paste0("ANNUAL_WB_",SCE,"_",GCM,"_Dunne"), full.names = FALSE); eta_list
eta_list <- eta_list[!grepl("NONCROP", eta_list)] # no NONCROP mapping for GCAM crops

yld_HS <- read.csv("C:/Model/gaia/agyield_impact_mri-esm2-0_r1i1p1f1_w5e5_food-ref.csv") %>%
  select(-X)


### CanESM5: Food-ref ----
GCM <- "CanESM5"
SCE <- "Food-ref"

eta_list <- list.files(input.dir, pattern = paste0("ANNUAL_WB_",SCE,"_",GCM,"_Dunne"), full.names = FALSE); eta_list

eta_list <- eta_list[!grepl("NONCROP", eta_list)] # no NONCROP mapping for GCAM crops

yld_HS <- read.csv("C:/Model/gaia/agyield_impact_canesm5_r1i1p1f1_w5e5_food-ref.csv") %>%
  select(-X)

### CanESM5: Food-95 ----
GCM <- "CanESM5"
SCE <- "Food-95"

eta_list <- list.files(input.dir, pattern = paste0("ANNUAL_WB_",SCE,"_",GCM,"_Dunne"), full.names = FALSE); eta_list

eta_list <- eta_list[!grepl("NONCROP", eta_list)] # no NONCROP mapping for GCAM crops

yld_HS <- read.csv("C:/Model/gaia/agyield_impact_canesm5_r1i1p1f1_w5e5_food-95.csv") %>%
  select(-X)




## start processing----

ETA <- list()
for (i in 1:length(eta_list)){
  eta <- readRDS(paste0(input.dir, eta_list[[i]]))
  df.eta <- do.call(rbind, eta) %>% filter(!is.na(region_id))
  basin_to_country_mapping %>% full_join(df.eta %>% rename(GCAM_basin_ID = region_id),
                                         by = "GCAM_basin_ID") -> df.eta.WB
  df.eta.WB$index = i
  ETA[[i]] <- df.eta.WB
  df.eta.WB
  print(dim(df.eta.WB))
}

eta.all <- do.call(rbind, ETA) %>%
  separate(crop, into = c("SPAM", "IRR_RFD"), sep = "_") # if include NONCROP, then NA values for IRR_RFD

eta.all %>% filter(!is.na(SPAM)) ->
  #Guam, Micronesia, Federated States of, French Southern Territories missing values
  eta.all

eta.all %>% filter(!is.na(smw)) %>%
  filter(SPAM != "NONCROP") %>% #
  left_join(crop_mapping, by = "SPAM", relationship = "many-to-many") %>%
  left_join(GCAM.sec, by = c("GCAM_Subsector")) %>%
  mutate(AgSupplySubsector = paste(GCAM_Subsector, GLU_name, sep = aglu.CROP_GLU_DELIMITER)) %>%
  group_by(GCAM_basin_ID, AgSupplySector, AgSupplySubsector, IRR_RFD, year) %>%
  summarise(value = weighted.mean(value, smw, na.rm = T)) %>%
  as_tibble() %>%
  repeat_add_columns(tibble(MGMT = c("hi", "lo"))) %>%
  mutate(IRR_RFD = ifelse(IRR_RFD == "I", "IRR", IRR_RFD),
         IRR_RFD = ifelse(IRR_RFD == "R", "RFD", IRR_RFD),
         AgProductionTechnology = paste(paste(AgSupplySubsector, IRR_RFD, sep = aglu.IRR_DELIMITER),
                                        MGMT, sep = aglu.MGMT_DELIMITER)) ->
  HS_PC

## smooth the heat stress shock ----
## define the window for moving average ---
library(data.table)
library(zoo)
HS_PC %>% select(-MGMT, -IRR_RFD) %>%
  # filter(!is.na(value)) %>%
  as.data.table() ->
  dt
# define moving average time window
K = 20
### Centered Moving Average with Adaptive Edges ----
setkey(dt, GCAM_basin_ID, AgSupplySector, AgSupplySubsector, AgProductionTechnology, year)
dt[, y_ma := rollapply(value, width = K, FUN = mean, fill = NA, partial = TRUE, align = "center"),
   by = .(GCAM_basin_ID, AgSupplySector, AgSupplySubsector, AgProductionTechnology)]

dt %>% as.data.frame() %>%
  mutate(value = ifelse(is.na(y_ma), 1, y_ma)) %>%
  filter(year %in% c(MODEL_FINAL_BASE_YEAR, MODEL_FUTURE_YEARS)) %>%
  group_by(AgSupplySector, AgSupplySubsector, AgProductionTechnology) %>%
  mutate(IO_mult = value[year == MODEL_FINAL_BASE_YEAR]/value) %>%  # IO_t / IO_2015 = PWC_2015 / PWC_t; value = PWC
  select(-y_ma) ->
  HS_PC_GCAM

saveRDS(HS_PC_GCAM, file = paste0("C:/Model/Claudia T/HELPS_GCAM_mapping_Dunne_",SCE,"_",GCM,".rds"))

# build XMLS ----


# read in gaia impact for land scaling ----

inputs_of("module_aglu_ag_prodchange_ref_IRR_MGMT_xml") %>% load_from_cache() -> all_data
L2052.AgProdChange_ag_irr_ref <- get_data(all_data, "L2052.AgProdChange_ag_irr_ref")
L2052.AgProdChange_bio_irr_ref <- get_data(all_data, "L2052.AgProdChange_bio_irr_ref")

yld_ref <- L2052.AgProdChange_ag_irr_ref %>%
  bind_rows(L2052.AgProdChange_bio_irr_ref)

yld_ref %>% mutate(scenario = "Ref") %>%
  bind_rows(yld_HS %>% mutate(scenario = GCM)) %>%
  arrange(scenario, region, AgSupplySector, AgSupplySubsector, AgProductionTechnology, year) %>%
  mutate(inter_mult = (1+AgProdChange)^5) %>%
  group_by(scenario, region, AgSupplySector, AgSupplySubsector, AgProductionTechnology) %>%
  mutate(mult = cumprod(inter_mult)) %>%
  group_by(region, AgSupplySector, AgSupplySubsector, AgProductionTechnology, year) %>%
  mutate(yld_mult = mult / mult[scenario == "Ref"],
         AIO_mult = 1 / yld_mult) %>%
  filter(scenario == GCM) %>%
  select(-inter_mult, -mult, -scenario) ->
  yld_IO_gaia


# read in HELPS labor impact ----
HS_PC_GCAM <- readRDS(paste0("C:/Model/Claudia T/HELPS_GCAM_mapping_Dunne_",SCE,"_",GCM,".rds"))

# read in XML gcamdata inputs ----
inputs_of("module_aglu_ag_input_IRR_MGMT_xml") %>% load_from_cache() -> all_data

MODULE_INPUTS <-
  c("L2062.AgCoef_Fert_ag_irr_mgmt",
    "L2062.AgCoef_Fert_bio_irr_mgmt",
    "L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA",
    "L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA",
    "L2082.AgCoef_laborcapital_for_tfp_MA")

get_data_list(all_data, MODULE_INPUTS, strip_attributes = TRUE)


outputs_of("module_aglu_L2082.ag_laborcapital_irr_mgmt") %>% load_from_cache() -> all_data

MODULE_OUTPUTS <-
  c("L2082.AgCoef_laborcapital_ag_irr_mgmt",
    "L2082.AgCoef_laborcapital_bio_irr_mgmt",
    "L2082.AgCoef_laborcapital_for")

get_data_list(all_data, MODULE_OUTPUTS, strip_attributes = TRUE)

## ag ----

L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA %>%
  filter(year >= 2025) %>%
  left_join_error_no_match(HS_PC_GCAM, by = c("AgSupplySector","AgSupplySubsector","AgProductionTechnology","year")) %>%
  mutate(coefficient = ifelse(minicam.energy.input == "Labor_Ag", coefficient*IO_mult, coefficient)) %>% # apply HS multiplier to labor IO
  ungroup() %>%
  select(colnames(L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA)) %>%
  bind_rows(L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA %>% filter(year < 2025)) ->
  L2082.AgCoef_laborcapital_ag_irr_mgmt_HS

L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA %>%
  filter(year >= 2025) %>%
  left_join_error_no_match(yld_IO_gaia, by = c("region", "AgSupplySector","AgSupplySubsector","AgProductionTechnology","year")) %>%
  mutate(coefficient = coefficient*AIO_mult) %>% # apply yld multiplier to labor and capital IO
  ungroup() %>%
  select(colnames(L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA)) %>%
  bind_rows(L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA %>% filter(year < 2025)) ->
  L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia

L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia %>%
  filter(year >= 2025) %>%
  left_join_error_no_match(HS_PC_GCAM, by = c("AgSupplySector","AgSupplySubsector","AgProductionTechnology","year")) %>%
  mutate(coefficient = ifelse(minicam.energy.input == "Labor_Ag", coefficient*IO_mult, coefficient)) %>% # apply HS multiplier to labor IO
  ungroup() %>%
  select(colnames(L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA)) %>%
  bind_rows(L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia %>% filter(year < 2025)) ->
  L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia_HS

## bio ----
L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA %>%
  filter(year >= 2025) %>%
  left_join_error_no_match(HS_PC_GCAM, by = c("AgSupplySector","AgSupplySubsector","AgProductionTechnology","year")) %>%
  mutate(coefficient = ifelse(minicam.energy.input == "Labor_Ag", coefficient*IO_mult, coefficient)) %>% # apply HS multiplier to labor IO
  ungroup() %>%
  select(colnames(L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA)) %>%
  bind_rows(L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA %>% filter(year < 2025)) ->
  L2082.AgCoef_laborcapital_bio_irr_mgmt_HS

L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA %>%
  filter(year >= 2025) %>%
  left_join_error_no_match(yld_IO_gaia, by = c("region", "AgSupplySector","AgSupplySubsector","AgProductionTechnology","year")) %>%
  mutate(coefficient = coefficient*AIO_mult) %>% # apply yld multiplier to labor and capital IO
  ungroup() %>%
  select(colnames(L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA)) %>%
  bind_rows(L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA %>% filter(year < 2025)) ->
  L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia

L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia %>%
  filter(year >= 2025) %>%
  left_join_error_no_match(HS_PC_GCAM, by = c("AgSupplySector","AgSupplySubsector","AgProductionTechnology","year")) %>%
  mutate(coefficient = ifelse(minicam.energy.input == "Labor_Ag", coefficient*IO_mult, coefficient)) %>% # apply HS multiplier to labor IO
  ungroup() %>%
  select(colnames(L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA)) %>%
  bind_rows(L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia %>% filter(year < 2025)) ->
  L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia_HS


# Write XMLs ----

# make sure to include the price.unit.conversion into XML
# otherwise, will see negative profit error in calibration years

# labor only: Dunne
create_xml(paste0("C:/Model/Claudia T/ag_input_IRR_MGMT_HELPS_",SCE,"_",GCM,"_Dunne.xml")) %>%
  add_xml_data(L2062.AgCoef_Fert_ag_irr_mgmt, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_ag_irr_mgmt_HS, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_ag_irr_mgmt_HS, "AgPriceConversion") %>%
  add_xml_data(L2062.AgCoef_Fert_bio_irr_mgmt, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_bio_irr_mgmt_HS, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_bio_irr_mgmt_HS, "AgPriceConversion") %>%
  add_xml_data(L2082.AgCoef_laborcapital_for_tfp_MA, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_for_tfp_MA, "AgPriceConversion") %>%
  add_precursors("L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA",
                 "L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA",
                 "L2082.AgCoef_laborcapital_for_tfp_MA",
                 "L2062.AgCoef_Fert_ag_irr_mgmt",
                 "L2062.AgCoef_Fert_bio_irr_mgmt") ->
  ag_input_IRR_MGMT_HS.xml
ag_input_IRR_MGMT_HS.xml %>% gcamdata::run_xml_conversion()


# crop only: land-scaling
create_xml(paste0("C:/Model/Claudia T/ag_input_IRR_MGMT_gaia_",SCE,"_",GCM,"_Dunne_LS.xml")) %>%
  add_xml_data(L2062.AgCoef_Fert_ag_irr_mgmt, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia, "AgPriceConversion") %>%
  add_xml_data(L2062.AgCoef_Fert_bio_irr_mgmt, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia, "AgPriceConversion") %>%
  add_xml_data(L2082.AgCoef_laborcapital_for_tfp_MA, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_for_tfp_MA, "AgPriceConversion") %>%
  add_precursors("L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA",
                 "L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA",
                 "L2082.AgCoef_laborcapital_for_tfp_MA",
                 "L2062.AgCoef_Fert_ag_irr_mgmt",
                 "L2062.AgCoef_Fert_bio_irr_mgmt") ->
  ag_input_IRR_MGMT_HS.xml
ag_input_IRR_MGMT_HS.xml %>% gcamdata::run_xml_conversion()


# crop&labor: land-scaling
create_xml(paste0("C:/Model/Claudia T/ag_input_IRR_MGMT_HELPS_",SCE,"_",GCM,"_Dunne_LS.xml")) %>%
  add_xml_data(L2062.AgCoef_Fert_ag_irr_mgmt, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia_HS, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_ag_irr_mgmt_gaia_HS, "AgPriceConversion") %>%
  add_xml_data(L2062.AgCoef_Fert_bio_irr_mgmt, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia_HS, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_bio_irr_mgmt_gaia_HS, "AgPriceConversion") %>%
  add_xml_data(L2082.AgCoef_laborcapital_for_tfp_MA, "AgCoef") %>%
  add_xml_data(L2082.AgCoef_laborcapital_for_tfp_MA, "AgPriceConversion") %>%
  add_precursors("L2082.AgCoef_laborcapital_ag_irr_mgmt_tfp_MA",
                 "L2082.AgCoef_laborcapital_bio_irr_mgmt_tfp_MA",
                 "L2082.AgCoef_laborcapital_for_tfp_MA",
                 "L2062.AgCoef_Fert_ag_irr_mgmt",
                 "L2062.AgCoef_Fert_bio_irr_mgmt") ->
  ag_input_IRR_MGMT_HS.xml
ag_input_IRR_MGMT_HS.xml %>% gcamdata::run_xml_conversion()


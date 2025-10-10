# install.packages("remotes")
# remotes::install_github("JGCRI/gcamdata")

conv_MIL_BIL = 1000.0
# conv_75_90 = 2.129
conv_75_90 = 2.129173
conv_75_15 = 3.507477
CONV_90_15 <- conv_75_15 / conv_75_90

# Load libs ----
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsci)
library(scales)
library(dplyr)
library(gcamdata)
library(purrr)
library(patchwork)
# library(broom)
library(sf)

source("R/LoadPackagesFuncs.R")
source("R/GCAM_module_funcs.R")

DIR_DATA <- "data"
DIR_OUTPUT <- "output"
DIR_MODULE = "Climate"

Project <- "HeatStress"
Version <- "Vclimate"
Scenario <- Load_GCAM(projnm = Project, versionnm = Version, return_availscen = T); Scenario

MODEL_FUTURE_YEARS  <- seq(2020, 2100, 5); MODEL_FUTURE_YEARS

# Check availability
Load_GCAM(projnm = Project, return_availversion = T)
Load_GCAM(projnm = Project, versionnm = Version, return_availscen = T)
Load_GCAM(projnm = Project, versionnm = Version, return_availquery = T)


# Modify/customize read csv function ----
read_csv_bind <- function(.multiCSVpath){
  
  library(doParallel)
  myCluster <-
    makeCluster(4, # number of cores to use
                type = "PSOCK") # type of cluster
  #detectCores()
  registerDoParallel(myCluster)
  
  foreach(csvDir = .multiCSVpath,
          .combine=rbind,
          .packages = "dplyr" ,.errorhandling = "remove"
  ) %dopar% {
    readr::read_csv(csvDir, skip = 1)%>%
      select(-matches("^X|\\...")) %>%
      na.omit() %>%
      filter(scenario != "scenario") %>%
      mutate(scenario = gsub(",date.*$", "", scenario)) %>%
      gcamdata::gather_years() %>%
      mutate(ss = sub(".*/([^/]+)/.*", "\\1", csvDir))
  } -> df
  
  stopCluster(myCluster)
  return(df)
}

rm(ListVclimate)
# Load everything into lists ----
Load_GCAM(projnm = Project, versionnm = "Vclimate", outputlistnm = "ListVclimate")

# create a project data output folder and save data
# dir.create(file.path(DIR_OUTPUT, Project, "ProjectRDS"), showWarnings = F) # somehow not working
ListVclimate %>% saveRDS(file.path(DIR_OUTPUT, Project, "ProjectRDS", paste0("ListVclimate", ".RDS")))

# Load the list [when needed]
ListVclimate <- readRDS(file.path(DIR_OUTPUT, Project, "ProjectRDS", paste0("ListVclimate", ".RDS")))

SCENARIO <- Scenario[!grepl("food_95_open|Ref", Scenario)]; SCENARIO
# SCENARIO <- Scenario[grepl("food_ref|Ref", Scenario)]; SCENARIO

PluckBind <- function(.query ){
  ListVclimate %>% purrr::pluck(.query) %>%
    mutate(branch = scenario, scenario = ss) %>%
    filter(scenario %in% SCENARIO)
}

PluckBind_climate <- function(.query ){
  ListVclimate %>% purrr::pluck(.query) %>%
    mutate(branch = scenario, scenario = ss) %>%
    filter(scenario %in% Scenario)
}

SCE_NAME <- function(.data ){
  .data %>% 
    mutate(ESM = "Ref",
           ESM = ifelse(grepl("MRI", scenario), "MRI", ESM),
           ESM = ifelse(grepl("CanESM5", scenario), "CanESM5", ESM),
           Hector = ifelse(grepl("food_ref", scenario), "Hector", "Ref"),
           Hector = ifelse(grepl("food_95", scenario), "Hector-M95", Hector),
           scenario = gsub("food_", "", scenario),
           scenario = gsub("ref_", "", scenario),
           scenario = gsub("95_", "", scenario),
           scenario = gsub("MRI_", "", scenario),
           scenario = gsub("CanESM5_", "", scenario),
           impact = scenario,
           scenario = paste0(Hector, "_", ESM))
}

REG <- st_read("data/maps/region_boundaries_moirai_combined_3p1_0p5arcmin.shp") 


## theme1 ----
theme1 <- theme(axis.text.x = element_text(angle = 40, hjust = 0.9, vjust = 1), legend.text.align = 0,
                strip.background = element_rect(fill="grey99"),
                strip.text = element_text(size = 12),
                axis.text.x.bottom = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linetype = 2, color = "grey80", size = 0.3),
                panel.spacing.y = unit(0.5, "lines"),
                panel.spacing.x = unit(0.5, "lines"))

## theme2 ----
theme2 <- theme(axis.text.x = element_text(angle = 40, hjust = 0.9, vjust = 1), legend.text.align = 0,
                strip.background = element_rect(fill="grey99"),
                strip.text = element_text(size = 12),
                axis.text.x.bottom = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.spacing.y = unit(0.5, "lines"),
                panel.spacing.x = unit(0.5, "lines"))


themeds <- theme(
  # panel.border = element_rect(colour = "black", size=1),
  text = element_text(family= fontfamily, size = 15),
  axis.text.y = element_text(angle = 0, color = "black", size = 15, margin = margin(r = 10)),
  axis.text.x = element_text(angle = 90, color = "black", size = 15, margin = margin(t = 10), vjust= 0.5),
  axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(size = 15, margin = margin(t = 10, r = 0, b = 0, l = 0))
)


# GMT ----
PluckBind("GMT") %>% 
  select(scenario, region, temperature, Units, year, value) -> 
  df.GMT

# Climate Forcing ----
PluckBind("TotalClimateForcing") %>% 
  select(scenario, region, Units, year, value) -> 
  df.forcing

# CO2 emission ----
PluckBind("CO2Emission") -> 
  df.CO2

# df.GMT %>% 
#   SCE_NAME() %>% 
#   filter(grepl("95", scenario)) %>% 
#   mutate(scenario = gsub("food_95_open", "Reference", scenario)) %>% 
#   group_by(year) %>% 
#   # filter(scenario != "HELPS_ref_open") %>% 
#   mutate(index = 100*(value / value[scenario == "Reference"] -1),
#          delta = value - value[scenario == "Reference"]) %>% 
#   select(scenario, year, value = delta) %>% 
#   mutate(variable = "GMT (degrees C)") %>% 
#   bind_rows(df.CO2 %>% 
#               group_by(scenario, Units, year) %>% 
#               summarise(value = sum(value, na.rm = T)) %>% 
#               filter(grepl("95", scenario)) %>% 
#               mutate(scenario = gsub("food_95_open", "Reference", scenario)) %>% 
#               group_by(year) %>% 
#               filter(year >= 2015) %>% 
#               mutate(index = 100*(value / value[scenario == "Reference"] -1),
#                      delta = value - value[scenario == "Reference"]) %>% 
#               select(scenario, year, value = index) %>% 
#               mutate(variable = "CO2 Emission (%)")) %>% 
#   bind_rows(df.forcing %>% 
#               filter(grepl("95", scenario)) %>% 
#               mutate(scenario = gsub("food_95_open", "Reference", scenario)) %>%  
#               group_by(year) %>% 
#               mutate(index = 100*(value / value[scenario == "Reference"] -1),
#                      delta = value - value[scenario == "Reference"]) %>% 
#               select(scenario, year, value = delta) %>% 
#               mutate(variable = "Climate Forcing (W/m^2)")) %>% 
#   bind_rows(df.GMT %>% 
#               filter(grepl("ref", scenario)) %>% 
#               mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
#               group_by(year) %>% 
#               # filter(scenario != "HELPS_ref_open") %>% 
#               mutate(index = 100*(value / value[scenario == "Reference"] -1),
#                      delta = value - value[scenario == "Reference"]) %>% 
#               select(scenario, year, value = delta) %>% 
#               mutate(variable = "GMT (degrees C)") %>% 
#               bind_rows(df.CO2 %>% 
#                           group_by(scenario, Units, year) %>% 
#                           summarise(value = sum(value, na.rm = T)) %>% 
#                           filter(grepl("ref", scenario)) %>% 
#                           mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
#                           group_by(year) %>% 
#                           filter(year >= 2015) %>% 
#                           mutate(index = 100*(value / value[scenario == "Reference"] -1),
#                                  delta = value - value[scenario == "Reference"]) %>% 
#                           select(scenario, year, value = index) %>% 
#                           mutate(variable = "CO2 Emission (%)")) %>% 
#               bind_rows(df.forcing %>% 
#                           filter(grepl("ref", scenario)) %>% 
#                           mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>%  
#                           group_by(year) %>% 
#                           mutate(index = 100*(value / value[scenario == "Reference"] -1),
#                                  delta = value - value[scenario == "Reference"]) %>% 
#                           select(scenario, year, value = delta) %>% 
#                           mutate(variable = "Climate Forcing (W/m^2)"))) %>% 
#   filter(year >= 2015,
#          !scenario %in% c("Reference", "Ref")) %>% 
#   SCE_NAME() ->
#   df.plot.95
# 
# df.plot.95 %>% 
#   ggplot(aes(x = year, y = value, color = impact)) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
#   geom_line(linewidth = 1) +
#   facet_grid(variable ~ scenario, scales = "free_y") +
#   scale_color_brewer(palette = "Set2") +
#   labs(x = "", y = "Response", title = "Response to climate impacts across scenarios") +
#   theme_bw()  ->
#   ClimateSummary; ClimateSummary
# 
# Write_png(ClimateSummary, "ClimateSummary", DIR_MODULE, w = 8, h = 6, r = 300)




# for paper ----
### CO2 emission by resource production ----

PluckBind("EmissionByResource") %>% 
  mutate(agg_sec = "fuel supply") %>% 
  group_by(scenario, year, region, agg_sec) %>% 
  summarise(value = sum(value, na.rm = T)) ->
  r_agg_sec_32

### CO2 emission from LUC ----
PluckBind("LUCEmission") %>% 
  mutate(agg_sec = "land use") %>% 
  group_by(scenario, year, region, agg_sec) %>% 
  summarise(value = sum(value, na.rm = T)) ->
  luc_agg_sec_32


### CO2 emission by sector ----
PluckBind("EmissionBySubsector") %>% 
  mutate(agg_sec = sector,
         agg_sec = ifelse(grepl("H2", sector), "fuel supply", agg_sec),
         agg_sec = ifelse(sector %in% c("refining", "gas pipeline", "gas processing"), "fuel supply", agg_sec),
         agg_sec = ifelse(grepl("resid", sector), "buildings", agg_sec),
         agg_sec = ifelse(grepl("trn_", sector), "transport", agg_sec),
         agg_sec = ifelse(grepl("elec", sector), "electricity", agg_sec),
         agg_sec = ifelse(grepl("backup", sector), "electricity", agg_sec),
         agg_sec = ifelse(grepl("regional", sector), "bio", agg_sec),
         agg_sec = ifelse(grepl("paper", sector), "industrial", agg_sec),
         agg_sec = ifelse(grepl("comm", sector), "buildings", agg_sec),
         agg_sec = ifelse(grepl("industrial", sector), "industrial", agg_sec),
         agg_sec = ifelse(grepl("use", sector), "industrial", agg_sec),
         agg_sec = ifelse(grepl("process heat", sector), "industrial", agg_sec),
         agg_sec = ifelse(grepl("chemical", sector), "industrial", agg_sec),
         agg_sec = ifelse(sector %in% c("ammonia", "cement","iron and steel", "district heat","desalinated water","alumina"), "industrial", agg_sec)) %>% 
  group_by(scenario, year, region, agg_sec) %>% 
  summarise(value = sum(value, na.rm = T)) ->
  s_agg_sec_32

s_agg_sec_32 %>% 
  bind_rows(r_agg_sec_32) %>% 
  bind_rows(luc_agg_sec_32) %>% 
  group_by(scenario, year, region, agg_sec) %>% 
  summarise(value = sum(value, na.rm = T)) ->
  agg_sec_32

agg_sec_32 %>% 
  mutate(agg_sec = ifelse(agg_sec %in% c("land use"), "land use", agg_sec),
         agg_sec = ifelse(agg_sec %in% c("bio", "buildings", "electricity", "fuel supply", "industrial", 
                                         "transport"), "Industry/Energy", agg_sec)) %>% 
  group_by(scenario, year, region, agg_sec) %>% 
  summarise(value = sum(value, na.rm = T)) %>% 
  filter(agg_sec == "Industry/Energy") %>% 
  group_by(region, year) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(delta = value - value[scenario == "Reference"],
         index = 100 * (value / value[scenario == "Reference"] - 1)) %>% 
  SCE_NAME() %>% 
  filter(impact == "all", 
         year == 2100) %>% 
  rename(reg_nm = region)->
  CO2_32

# breaks <- seq(-60, 60, by = 10)
# colors <- scales::div_gradient_pal(low = "blue", mid = "yellow", high = "red")(seq(0, 1, length.out = length(breaks)))


REG %>% 
  left_join(CO2_32, by = "reg_nm") %>% 
  ggplot() +
  geom_sf(aes(fill = delta)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "MtC") +
  facet_wrap(~ scenario) +
  labs(title = "Industry/Energy CO2 Emissions relative to BwO: 2100") +
  coord_sf(datum = NA) +
  theme_bw() + theme0 + theme1 ->
IndustryCO2MtC

Write_png(IndustryCO2MtC, "IndustryCO2MtC", DIR_MODULE, w = 8, h = 4, r = 300)


REG %>% 
  left_join(CO2_32, by = "reg_nm") %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
  facet_wrap(~ scenario) +
  labs(title = "Industry/Energy CO2 Emissions relative to BwO: 2100") +
  coord_sf(datum = NA) +
  theme_bw() + theme0 + theme1 ->
IndustryCO2Pct

Write_png(IndustryCO2Pct, "IndustryCO2Pct", DIR_MODULE, w = 8, h = 4, r = 300)


#### one-at-a-time ----

agg_sec_32 %>% 
  mutate(agg_sec = ifelse(agg_sec %in% c("land use"), "land use", agg_sec),
         agg_sec = ifelse(agg_sec %in% c("bio", "buildings", "electricity", "fuel supply", "industrial", 
                                         "transport"), "Industry/Energy", agg_sec)) %>% 
  group_by(scenario, year, region, agg_sec) %>% 
  summarise(value = sum(value, na.rm = T)) %>% 
  group_by(region, year, agg_sec) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(delta = value - value[scenario == "Reference"],
         index = 100 * (value / value[scenario == "Reference"] - 1)) %>% 
  filter(year == 2100) %>% rename(reg_nm = region) %>% 
  SCE_NAME() ->
  CO2_32_sce


REG %>% left_join(CO2_32_sce %>% filter(impact == "crop", agg_sec == "Industry/Energy"), by = "reg_nm") %>%
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "hdcd", agg_sec == "Industry/Energy"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "hydro", agg_sec == "Industry/Energy"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "labor", agg_sec == "Industry/Energy"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "runoff", agg_sec == "Industry/Energy"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "all", agg_sec == "Industry/Energy"), by = "reg_nm")) %>% 
  ggplot() +
  geom_sf(aes(fill = delta)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "MtC") +
  labs(title = "Industry/Energy CO2 Emissions relative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(impact~scenario) +
  theme_bw() + theme0 + theme1 ->
IndustryCO2MtC_individual

Write_png(IndustryCO2MtC_individual, "IndustryCO2MtC_individual", DIR_MODULE, w = 10, h = 6, r = 300)


## LUC emission ----
agg_sec_32 %>% 
  mutate(agg_sec = ifelse(agg_sec %in% c("land use"), "land use", agg_sec),
         agg_sec = ifelse(agg_sec %in% c( "bio","buildings", "electricity", "fuel supply", "industrial", 
                                          "transport"), "Industry/Energy", agg_sec)) %>% 
  group_by(scenario, year, region, agg_sec) %>% 
  summarise(value = sum(value, na.rm = T)) %>% 
  filter(agg_sec == "land use") %>% 
  group_by(region, year) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(delta = value - value[scenario == "Reference"],
         index = 100 * (value / value[scenario == "Reference"] - 1)) %>% 
  SCE_NAME() %>% 
  filter(impact == "all", year == 2100) %>% rename(reg_nm = region)->
  CO2_32

REG %>% 
  left_join(CO2_32, by = "reg_nm") %>% 
  ggplot() +
  geom_sf(aes(fill = delta)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "MtC") +
  facet_wrap(~ scenario) +
  labs(title = "LUC CO2 Emissions relative to BwO: 2100") +
  coord_sf(datum = NA) +
  theme_bw() + theme0 + theme1 ->
LUCCO2MtC

Write_png(LUCCO2MtC, "LUCCO2MtC", DIR_MODULE, w = 8, h = 4, r = 300)


REG %>% 
  left_join(CO2_32, by = "reg_nm") %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
  facet_wrap(~ scenario) +
  labs(title = "LUC CO2 Emissions relative to BwO: 2100") +
  coord_sf(datum = NA) +
  theme_bw() + theme0 + theme1 ->
LUCCO2Pct

Write_png(LUCCO2Pct, "LUCCO2Pct", DIR_MODULE, w = 8, h = 4, r = 300)


REG %>% left_join(CO2_32_sce %>% filter(impact == "crop", agg_sec == "land use"), by = "reg_nm") %>%
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "hdcd", agg_sec == "land use"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "hydro", agg_sec == "land use"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "labor", agg_sec == "land use"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "runoff", agg_sec == "land use"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(CO2_32_sce %>% filter(impact == "all", agg_sec == "land use"), by = "reg_nm")) %>% 
  ggplot() +
  geom_sf(aes(fill = delta)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "MtC") +
  labs(title = "LUC CO2 Emissions relative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(impact~scenario) +
  theme_bw() + theme0 + theme1 ->
  LUCCO2MtC_individual

Write_png(LUCCO2MtC_individual, "LUCCO2MtC_individual", DIR_MODULE, w = 10, h = 6, r = 300)

# land ----

PluckBind("Aggland") %>% filter(year >= 2015) %>% 
  left_join_error_no_match(LandMapping %>% select(LandLeaf, land = LandCover4), by = "LandLeaf") %>%
  group_by(scenario, region, land, year, branch) %>%
  # to Mha
  summarise(value = sum(value)/10, .groups = "drop") %>%
  Agg_reg(land, region) %>% filter(year >= 2015) %>% 
  filter(!grepl("Fixed", land)) %>%
  mutate(land = ifelse(grepl("Cropland", land), "Cropland" , land),
         land = ifelse(grepl("Forest", land), "Forest" , land),
         land = ifelse(grepl("Other Natural|Pasture", land), "Other Natural" , land)) %>% 
  group_by_at(vars(-value)) %>% summarise(value = sum(value), .groups = "drop") -> Pland

Pland %>% group_by(land, region, year) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(delta = value - value[scenario == "Reference"],
         index = 100 * (value / value[scenario == "Reference"] - 1)) %>% 
  filter(year == 2100) %>% 
  rename(reg_nm = region) ->
  df.ag.land

REG %>% left_join(df.ag.land, by = "reg_nm") %>%
  bind_rows(REG %>% left_join(df.ag.land %>% filter(scenario == "Forest"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(df.ag.land %>% filter(scenario == "Other Natural"), by = "reg_nm")) %>% 
  filter(!scenario %in% c("Ref", "Reference")) %>% 
  SCE_NAME() %>% na.omit() ->
  delta.ag.land

delta.ag.land %>% 
  filter(impact == "all") %>% 
  ggplot() +
  geom_sf(aes(fill = delta)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "Mha") +
  labs(title = "Land use relative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(land ~ scenario) +
  theme_bw() + theme0 + theme1 ->
LandUse

Write_png(LandUse, "LandUse", DIR_MODULE, w = 9, h = 4, r = 300)

# residential energy demand ----

PluckBind("ResidFinalEnergy") %>% 
  filter(grepl('(resid cooling|resid heating)', sector)) %>% 
  mutate(DD = if_else(grepl('resid cooling', sector), 'cooling', 'heating')) ->
  df.HDCD

df.HDCD %>% 
  group_by(scenario, year, region, DD) %>% 
  summarise(value = sum(value, na.rm = T)) %>% 
  group_by(year, region, DD) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(delta = value - value[scenario == "Reference"],
         index = 100 * (value / value[scenario == "Reference"] - 1)) %>% 
  filter(year == 2100) %>% rename(reg_nm = region)->
  HDCD_32

REG %>% left_join(HDCD_32 %>% filter(DD == "heating"), by = "reg_nm") %>% 
  bind_rows(REG %>% left_join(HDCD_32 %>% filter(DD == "cooling"), by = "reg_nm")) %>% 
  filter(!scenario %in% c("Ref", "Reference")) %>% 
  SCE_NAME() %>% na.omit() ->
  delta.HDCD

delta.HDCD %>% 
  filter(impact == "all") %>% 
  ggplot() +
  geom_sf(aes(fill = delta)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "EJ") +
  labs(title = "Residential energy demand\nrelative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(DD ~ scenario) +
  theme_bw() + theme0 + theme1 ->
HDCDEJ

Write_png(HDCDEJ, "HDCDEJ", DIR_MODULE, w = 8, h = 4, r = 300)


# index has an outlier, so the maps are pretty much white
# need a better way to present the results
# delta.HDCD %>% 
#   ggplot() +
#   geom_sf(aes(fill = index)) +
#   scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
#   labs(title = "Residential energy demand\nrelative to BwO: 2100") +
#   coord_sf(datum = NA) +
#   facet_grid(DD~scenario) +
#   theme_bw() + theme0 + theme1 ->
# HDCDPct
# 
# Write_png(HDCDPct, "HDCDPct", DIR_MODULE, w = 8, h = 4, r = 300)



PluckBind("SAM_NA") %>% 
  # select(scenario, region, Account, year, value) %>% 
  filter(Account== "gdp-per-capita") -> 
  df.income

df.income %>% select(scenario, region, year, GDPpc = value) %>% head

# df.HDCD %>% select(scenario, region, sector, DD, year, EJ = value) %>% 
#   SCE_NAME() %>% 
#   filter(year >= 2015) %>% 
#   filter(impact =="all") %>% 
#   filter(DD == "cooling") %>% 
#   mutate(sector = factor(sector, levels = paste0("resid cooling modern_d", 1:10))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = EJ, fill = sector), stat = "identity", position = "stack") +
#   facet_wrap(~ region, ncol = 8) +
#   scale_fill_brewer(palette = "RdYlGn") +
#   theme_bw()

# df.HDCD %>% select(scenario, region, sector, DD, year, EJ = value) %>% 
#   filter(year >= 2015) %>% 
#   filter(scenario == "all_crop_LK") %>%
#   filter(DD == "heating") %>% 
#   mutate(source = ifelse(grepl("TradBio", sector), "TradBio", sector),
#          source = ifelse(grepl("modern", sector), "modern", source),
#          source = ifelse(grepl("coal", sector), "coal", source)) %>% 
#   filter(source == "coal") %>% 
#   mutate(sector = gsub("resid heating ", "", sector),
#          sector = factor(sector, levels = paste0("coal_d", 1:10))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = EJ, fill = sector), stat = "identity", position = "stack") +
#   facet_wrap(~ region, ncol = 8) +
#   scale_fill_brewer(palette = "RdYlGn") +
#   theme_bw()
# 
# df.HDCD %>% select(scenario, region, sector, DD, year, EJ = value) %>% 
#   filter(year >= 2015) %>% 
#   filter(scenario == "all_crop_LK") %>%
#   filter(DD == "heating") %>% 
#   mutate(source = ifelse(grepl("TradBio", sector), "TradBio", sector),
#          source = ifelse(grepl("modern", sector), "modern", source),
#          source = ifelse(grepl("coal", sector), "coal", source)) %>% 
#   filter(source == "modern") %>% 
#   mutate(sector = gsub("resid heating ", "", sector),
#          sector = factor(sector, levels = paste0("modern_d", 1:10))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = EJ, fill = sector), stat = "identity", position = "stack") +
#   facet_wrap(~ region, ncol = 8) +
#   scale_fill_brewer(palette = "RdYlGn") +
#   theme_bw()
# 
# df.HDCD %>% select(scenario, region, sector, DD, year, EJ = value) %>% 
#   filter(year >= 2015) %>% 
#   filter(scenario == "all_crop_LK") %>%
#   filter(DD == "heating") %>% 
#   mutate(source = ifelse(grepl("TradBio", sector), "TradBio", sector),
#          source = ifelse(grepl("modern", sector), "modern", source),
#          source = ifelse(grepl("coal", sector), "coal", source)) %>% 
#   filter(source == "TradBio") %>% 
#   mutate(sector = gsub("resid heating ", "", sector),
#          sector = factor(sector, levels = paste0("TradBio_d", 1:10))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = EJ, fill = sector), stat = "identity", position = "stack") +
#   facet_wrap(~ region, ncol = 8) +
#   scale_fill_brewer(palette = "RdYlGn") +
#   theme_bw()



# price index ----

PluckBind("Agprices") %>% 
  group_by(region, sector, year) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(index = 100 * value / value[scenario == "Reference"] - 100) %>% 
  filter(year == 2100) %>% rename(reg_nm = region) ->
  Ag.Price.32

REG %>% left_join(Ag.Price.32 %>% filter(sector == "Corn"), by = "reg_nm") %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% filter(sector == "Wheat"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% filter(sector == "Soybean"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% filter(sector == "Rice"), by = "reg_nm")) %>% 
  filter(!scenario %in% c("Ref", "Reference")) %>% 
  SCE_NAME() %>% na.omit() ->
  delta.Ag.Price

delta.Ag.Price %>% 
  filter(impact == "all") %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
  labs(title = "Major crop price\nrelative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(sector ~ scenario) +
  theme_bw() + theme0 + theme1 ->
  Ag.Price

Write_png(Ag.Price, "Ag.Price", DIR_MODULE, w = 8, h = 6, r = 300)

### one-at-a-time ----
sec_list <- c("Corn", "Rice", "Vegetables", "Wheat", "Soybean")

for (f in sec_list){ 

REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "crop", sector == f), by = "reg_nm") %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "hdcd", sector == f), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "hydro", sector == f), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "labor", sector == f), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "runoff", sector == f), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "all", sector == f), by = "reg_nm")) %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
  labs(title = paste0(f, " price relative to BwO: 2100")) +
  coord_sf(datum = NA) +
  facet_grid(impact~scenario) +
  theme_bw() + theme0 + theme1 ->
  price_individual

Write_png(price_individual, paste0(f, "_price_individual"), DIR_MODULE, w = 10, h = 6, r = 300)

}



PluckBind("Meatprices") %>% 
  group_by(region, sector, year) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(index = 100 * value / value[scenario == "Reference"] - 100) %>% 
  filter(year == 2100) %>% rename(reg_nm = region) ->
  Mt.Price.32

REG %>% left_join(Mt.Price.32 %>% filter(sector == "Beef"), by = "reg_nm") %>% 
  bind_rows(REG %>% left_join(Mt.Price.32 %>% filter(sector == "Dairy"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Mt.Price.32 %>% filter(sector == "Pork"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Mt.Price.32 %>% filter(sector == "Poultry"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Mt.Price.32 %>% filter(sector == "SheepGoat"), by = "reg_nm")) %>% 
  filter(!scenario %in% c("Ref", "Reference")) %>% 
  SCE_NAME() %>% na.omit() ->
  delta.Mt.Price

delta.Ag.Price %>% 
  filter(impact == "all") %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
  labs(title = "Major livestock price\nrelative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(sector ~ scenario) +
  theme_bw() + theme0 + theme1 ->
  Mt.Price

Write_png(Mt.Price, "Mt.Price", DIR_MODULE, w = 8, h = 6, r = 300)

### one-at-a-time ----
sec_list <- c("Beef", "Dairy", "Pork", "Poultry", "SheepGoat")

for (f in sec_list){ 
  
  REG %>% left_join(Mt.Price.32 %>% SCE_NAME() %>% filter(impact == "crop", sector == f), by = "reg_nm") %>% 
    bind_rows(REG %>% left_join(Mt.Price.32 %>% SCE_NAME() %>% filter(impact == "hdcd", sector == f), by = "reg_nm")) %>% 
    bind_rows(REG %>% left_join(Mt.Price.32 %>% SCE_NAME() %>% filter(impact == "hydro", sector == f), by = "reg_nm")) %>% 
    bind_rows(REG %>% left_join(Mt.Price.32 %>% SCE_NAME() %>% filter(impact == "labor", sector == f), by = "reg_nm")) %>% 
    bind_rows(REG %>% left_join(Mt.Price.32 %>% SCE_NAME() %>% filter(impact == "runoff", sector == f), by = "reg_nm")) %>% 
    bind_rows(REG %>% left_join(Mt.Price.32 %>% SCE_NAME() %>% filter(impact == "all", sector == f), by = "reg_nm")) %>% 
    ggplot() +
    geom_sf(aes(fill = index)) +
    scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
    labs(title = paste0(f, " price relative to BwO: 2100")) +
    coord_sf(datum = NA) +
    facet_grid(impact~scenario) +
    theme_bw() + theme0 + theme1 ->
    price_individual
  
  Write_png(price_individual, paste0(f, "_price_individual"), DIR_MODULE, w = 10, h = 6, r = 300)
  
}



REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "crop", sector == "Corn"), by = "reg_nm") %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "hdcd", sector == "Corn"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "hydro", sector == "Corn"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "labor", sector == "Corn"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "runoff", sector == "Corn"), by = "reg_nm")) %>% 
  bind_rows(REG %>% left_join(Ag.Price.32 %>% SCE_NAME() %>% filter(impact == "all", sector == "Corn"), by = "reg_nm")) %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
  labs(title = "Corn price relative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(impact~scenario) +
  theme_bw() + theme0 + theme1 ->
  price_individual

Write_png(price_individual, "Corn_price_individual", DIR_MODULE, w = 10, h = 6, r = 300)



PluckBind("ElecPriceBySector") %>% 
  group_by(region, fuel, year) %>%
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(index = 100 * value / value[scenario == "Reference"] - 100) %>% 
  filter(year == 2100) %>% rename(reg_nm = region) ->
  Elec.Price.32

REG %>% left_join(Elec.Price.32 %>% filter(fuel == "electricity"), by = "reg_nm") %>% 
  # bind_rows(REG %>% left_join(Ag.Price.32 %>% filter(fuel == "elect_td_ind"), by = "reg_nm")) %>% 
  # bind_rows(REG %>% left_join(Ag.Price.32 %>% filter(fuel == "elect_td_trn"), by = "reg_nm")) %>% 
  # bind_rows(REG %>% left_join(Ag.Price.32 %>% filter(fuel == "elect_td_bld"), by = "reg_nm")) %>% 
  filter(!scenario %in% c("Ref", "Reference")) %>% 
  SCE_NAME() %>% na.omit() ->
  delta.Elec.Price

delta.Elec.Price %>% 
  filter(impact == "all") %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "red", high = "blue", midpoint = 0, name = "%") +
  labs(title = "Electricity price\nrelative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_wrap( ~ scenario) +
  theme_bw() + theme0 + theme1 ->
  ELec.Price

Write_png(ELec.Price, "ELec.Price", DIR_MODULE, w = 6, h = 4, r = 300)


# STOP HERE !!!! -----

# non-CO2 Emissions ----
PluckBind("NonCO2Emission") %>% 
  filter(grepl('CH4|SO2', GHG)) %>% 
  mutate(GHG = ifelse(grepl("CH4", GHG), "CH4", "SO2")) %>% 
  group_by(scenario, region, GHG, Units, year) %>% 
  summarise(value = sum(value, na.rm = T)) %>% 
  group_by(region, GHG, Units, year) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>% 
  mutate(delta = value - value[scenario == "Reference"],
         index = 100 * (value / value[scenario == "Reference"] - 1)) %>% 
  filter(year == 2100) %>% rename(reg_nm = region)->
  NonCO2_32

REG %>% left_join(NonCO2_32 %>% filter(GHG == "CH4"), by = "reg_nm") %>% 
  bind_rows(REG %>% left_join(NonCO2_32 %>% filter(GHG == "SO2"), by = "reg_nm")) %>% 
  filter(!scenario %in% c("Ref", "Reference")) %>% 
  SCE_NAME() %>% na.omit() ->
  delta.NonCO2

delta.NonCO2 %>% 
  filter(impact == "all") %>% 
  ggplot() +
  geom_sf(aes(fill = delta)) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, name = "Tg") +
  labs(title = "NonCO2 Emissions\nrelative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(GHG ~ scenario) +
  theme_bw() + theme0 + theme1 ->
  NonCO2Tg

Write_png(NonCO2Tg, "NonCO2Tg", DIR_MODULE, w = 8, h = 4, r = 300)


delta.NonCO2 %>% 
  filter(impact == "all") %>% 
  ggplot() +
  geom_sf(aes(fill = index)) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, name = "%") +
  labs(title = "NonCO2 Emissions \nrelative to BwO: 2100") +
  coord_sf(datum = NA) +
  facet_grid(GHG ~ scenario) +
  theme_bw() + theme0 + theme1 ->
  NonCO2Pct

Write_png(NonCO2Pct, "NonCO2Pct", DIR_MODULE, w = 8, h = 4, r = 300)


PluckBind("NonCO2Emission") %>% 
  filter(grepl('CH4|SO2', GHG)) %>% 
  mutate(GHG = ifelse(grepl("CH4", GHG), "CH4", "SO2")) %>% 
  group_by(scenario,GHG, Units, year) %>% 
  summarise(value = sum(value, na.rm = T)) %>% 
  group_by(GHG, Units, year) %>% 
  mutate(scenario = gsub("food_ref_open", "Reference", scenario)) %>%
  filter(scenario != "food_95_open") %>% 
  mutate(delta = value - value[scenario == "Reference"],
         index = 100 * (value / value[scenario == "Reference"] - 1)) ->
  NonCO2_glb

NonCO2_glb %>% 
  SCE_NAME() %>% 
  filter(impact != "Reference") %>% 
  filter(year >= 2015) %>% 
  ggplot() +
  geom_line(aes(x = year, y = index, color = impact)) +
  facet_grid(GHG ~ scenario) +
  theme_bw() +
  labs(x = "", y = "%", title = "Global NonCO2 Emission Chnages due to Climate Change")

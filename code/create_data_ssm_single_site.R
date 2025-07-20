create_data_ssm_single_site <- function(cwac_data, steps = c("missing", "model"), ...){

  varargs <- list(...)

  sp_code <- cwac_data |>
    dplyr::filter(!is.na(SppRef)) |>
    dplyr::pull(SppRef) |>
    unique()

  # Add missing counts ------------------------------------------------------

  message("Adding missing counts and guessing dates for them")

  # Figure out which counts are missing and which are zero for the species of interest
  site_data_w_miss <- CWAC::addMissingCwacCounts(cwac_data, min(cwac_data$Year):max(cwac_data$Year))

  # Give missing surveys a date based on the dates from other surveys
  month_summer <- cwac_data %>%
    dplyr::mutate(month = lubridate::month(StartDate)) %>%
    dplyr::filter(Season == "S", !is.na(month)) %>%
    dplyr::count(month) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::pull(month)

  if(length(month_summer) == 0) month_summer <- 1

  month_winter <- cwac_data %>%
    dplyr::mutate(month = lubridate::month(StartDate)) %>%
    dplyr::filter(Season == "W", !is.na(month)) %>%
    dplyr::count(month) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::pull(month)

  if(length(month_winter) == 0) month_winter <- 7

  site_data_w_miss <- site_data_w_miss %>%
    dplyr::mutate(StartDate = dplyr::case_when(is.na(StartDate) & Season == "S" ~ as.Date(paste(Year, month_summer, "01", sep = "-")),
                                               is.na(StartDate) & Season == "W" ~ as.Date(paste(Year, month_winter, "01", sep = "-")),
                                               TRUE ~ StartDate)) %>%
    dplyr::arrange(StartDate, Season)

  # Prepare counts for the species and site of interest
  counts <- prepSsmData(counts = site_data_w_miss,
                        spp_sel = sp_code,
                        keep = c("Country", "LocationCode", "CountCondition", "X", "Y", "Year"))


  # Prepare data for modelling ----------------------------------------------

  message("Preparing data for modelling: removing seasons other than S and W")

  # Remove seasons other than summer and winter
  counts_mod <- counts %>%
    dplyr::filter(Season %in% c("S", "W"))

  # Create other useful variables
  counts_mod <- counts_mod %>%
    dplyr::mutate(season_id = dplyr::case_when(Season == "S" ~ 1,
                                               Season == "W" ~ 2,
                                               TRUE ~ 3),
                  date = lubridate::date(StartDate)) %>%
    dplyr::group_by(LocationCode) %>%
    dplyr::mutate(site_id = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(year_id = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site_id, year_id, season_id) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(visit_id = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(site_id, year_id, season_id, visit_id) %>%
    dplyr::mutate(spp_id = sp_code)

  # Order data by: site, year, season
  counts_mod <- counts_mod %>%
    dplyr::arrange(site_id, year_id, season_id)

  return(counts_mod)

}

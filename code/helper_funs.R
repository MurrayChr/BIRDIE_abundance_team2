getCwacModelData_sp <- function(
    site_code,
    genus,
    specific_epithet){
  raw_dat <- getCwacSiteCounts(site_code) %>%
    filter(Genus == genus & Species == sp_epithet)
  out <- create_data_ssm_single_site(raw_dat)
  return(out)
}
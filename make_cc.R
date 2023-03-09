# Define case-crossover dataset generating function -----------------------
# t0: original time-points (days)
# z0: original exposure (days)
# t1: vector of case time-points (may include repeats if multiple cases on same day)
# x1: matrix of covariates with same number of rows as t1
make_cc <- function(t0, z0, t1, x1){
  
  # Create dataset containing all strata
  all_strata <- data.frame(Date = t0, z0) %>%
    mutate(Year = lubridate::year(Date),
           Month = lubridate::month(Date),
           DoW = weekdays(Date),
           Strata = as.numeric(interaction(DoW, Month)))
  
  CC <- data.frame(Case.Date = t1, x1) %>%
    mutate(ID = row_number()) %>%
    left_join(select(all_strata, Date, Strata), by = c('Case.Date' = 'Date')) %>% # Determine the strata to which the case belongs to
    left_join(all_strata, by = join_by(Strata), multiple = 'all') %>% # Add records for the remaining observations in the case strata
    filter(!is.na(Case.Date)) %>% # Remove strata that don't have cases
    mutate(Y = ifelse(Date == Case.Date, 1, 0)) %>% # Identify the record for which the case occurred
    select(-c(Case.Date, Strata))
  
  return(CC)
}

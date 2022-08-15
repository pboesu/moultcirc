#package simulated data
library(dplyr)
sim_data <- readRDS('../../2022_moult_methods_paper/data/simulated_datasets/simulated_no_recaptures_start_duration_temp_no_no_linear_set100.rds') %>% ungroup()

sim_data_small <- slice_sample(sim_data, n = 500) %>% mutate(yday_shifted = (yday + 150) %% 365,
                                                             circday = yday*(2*pi / 365 ) - pi,
                                                             circday_shifted = ifelse(circday + 2.5 < pi, circday + 2.5, circday + 2.5 - 2*pi)) %>% select(yday, yday_shifted,moult_score)

usethis::use_data(sim_data_small, overwrite = TRUE)

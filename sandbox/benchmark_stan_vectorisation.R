#\\benchmark vectorisation
library(moultcirc)
library(microbenchmark)

data("sim_data_small")

#test_that("model runs through for moult contiguous in [-pi,pi]", {

mbm = microbenchmark(optT = uz2_circ("moult_score", "yday", data = sim_data_small[1:100,], chains = 2, parallel_chains = 2, refresh = 1000),
                optF = uz2_circ("moult_score", "yday", data = sim_data_small[1:100,], chains = 2, parallel_chains = 2, refresh = 1000, opt = FALSE),
                optT1k = uz2_circ("moult_score", "yday", data = sim_data_small, chains = 2, parallel_chains = 2, refresh = 1000),
                optF1k = uz2_circ("moult_score", "yday", data = sim_data_small, chains = 2, parallel_chains = 2, refresh = 1000, opt = FALSE),
                optTs = uz2_circ("moult_score", "yday_shifted", data = sim_data_small[1:100,], chains = 2, parallel_chains = 2, refresh = 1000),
                optFs = uz2_circ("moult_score", "yday_shifted", data = sim_data_small[1:100,], chains = 2, parallel_chains = 2, refresh = 1000, opt = FALSE),
                times = 10)

saveRDS(mbm, 'sandbox/mbm_uz2_circ.rds')
pdf('sandbox/mbm_uz2_circ.pdf')
boxplot(mbm)
dev.off()

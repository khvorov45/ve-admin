# Summarising simulation results
# Arseniy Khvorov
# Created 2019/11/19
# Last edit 2019/11/19

library(tidyverse)
library(tools)

sim_dir <- "sim"
sim_summ_dir <- "sim-summary"

# Functions ===================================================================

read_res <- function(filepath) {
  read_csv(filepath, col_types = cols())
}

summ_res <- function(res) {
  nongroup <- c("ve_est", "n_study", "seed", "index")
  grouping_vars <- names(res)[!names(res) %in% nongroup]
  res %>%
    group_by(!!!syms(grouping_vars)) %>%
    summarise(
      ve_est_mean = mean(ve_est),
      ve_est_sd = sd(ve_est),
      n_study_mean = mean(n_study)
    ) %>%
    select(
      name, vary_name, type, ve_est_mean, ve_est_sd, n_study_mean, everything()
    )
}

save_summ <- function(summ, name, folder) {
  write_csv(summ, file.path(folder, paste0(name, ".csv")))
}

# Script ======================================================================

all_res_nms <- list_files_with_exts(sim_dir, "csv")
all_res <- map(all_res_nms, read_res)
names(all_res) <- str_replace(basename(all_res_nms), ".csv", "")

all_summs <- map(all_res, summ_res)

iwalk(all_summs, save_summ, sim_summ_dir)

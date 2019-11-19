# Simulations of VE studies
# Arseniy Khvorov
# Created 2019/11/13
# Last edit 2019/11/13

library(tidyverse)
library(rlang)
library(furrr) # Like purrr but with parallel support
library(extraDistr) # For the categorical and bernoulli disributions

plan(multiprocess, workers = 20) # Parallel setup

# Directory names that will be used later
sim_dir <- "sim"
pars_dir <- "pars"

# Settings ====================================================================

nsim <- 5 # Number of simulations
pars_dict <- read_csv(file.path(pars_dir, "pars.csv")) # Parameter values

# These will be varied one at a time in a population that consists of one group
# Results are saved to one-<nsim>sims.csv
vary_table <- list(
  pvac = seq(0.05, 0.5, 0.05),
  pflu = seq(0.05, 0.15, 0.01),
  ve = seq(0.1, 0.9, 0.1),
  pnonflu = seq(0.1, 0.3, 0.02),
  sympt = seq(0.1, 0.9, 0.1),
  clin = seq(0.1, 0.9, 0.1),
  test_clin = seq(0.1, 0.9, 0.1),
  test_nonclin = seq(0, 0.3, 0.05),
  sens_vac = seq(0.9, 1, 0.01),
  spec_vac = seq(0.5, 1, 0.05),
  sens_flu = seq(0.5, 1, 0.05),
  spec_flu = seq(0.9, 1, 0.01)
)

# These will be used for simulations of populations with multiple groups.
# Pattern is (low, mid, high). All groups are set to one of 6 patterns:
# all low, all mid, all high, 1 low 2 high (3 groups, so 3 of these).
# Results are saved to mult-<nsim>.csv
vary_table_mult <- list(
  prop = c(0.33, 0.7, 0.15),
  pvac = c(0.05, 0.3, 0.5),
  pflu = c(0.05, 0.1, 0.15),
  ve = c(0.15, 0.33, 0.7),
  pnonflu = c(0.1, 0.15, 0.3),
  sympt = c(0.1, 0.5, 0.9),
  clin = c(0.1, 0.5, 0.9),
  test_clin = c(0.1, 0.5, 0.9),
  test_nonclin = c(0, 0.15, 0.3),
  sens_vac = c(0.9, 0.95, 1),
  spec_vac = c(0.5, 0.75, 1),
  sens_flu = c(0.5, 0.75, 1),
  spec_flu = c(0.9, 0.95, 1)
)

# These will be varied simultaneously in a population with one group in it to
# investigate the interaction between ve and pvac
# in presence of misclassification by running it on different populations
# (parameter sets) which have different miscalssification profiles.
vary_table_ve <- list(
  ve = seq(0.1, 0.9, 0.1),
  pvac = seq(0.1, 0.9, 0.1)
)

# Functions ===================================================================

# Population simulation
simulate_pop <- function(nsam = 1000,
                         pvac = 0.5,
                         pflu = 0.1,
                         ve = 0.5,
                         pnonflu = 0.1,
                         sympt = 1,
                         clin = 1,
                         test_clin = 1,
                         test_nonclin = 1,
                         sens_vac = 1,
                         spec_vac = 1,
                         sens_flu = 1,
                         spec_flu = 1,
                         seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)

  # Work out infection status
  calc_status <- function(n, flu, nonflu) {
    healthy <- 1 - flu - nonflu
    probs <- matrix(c(flu, nonflu, healthy), ncol = 3)
    rcat(n, probs) %>% recode("1" = "Flu", "2" = "Nonflu", "3" = "Healthy")
  }

  # Simulate data
  pop <- tibble(
    vaccinated = rbern(nsam, pvac),
    vac_record_prob = if_else(vaccinated == 1, sens_vac, 1 - spec_vac),
    vac_record = rbern(nsam, vac_record_prob),
    flu_prob = if_else(vaccinated == 1, pflu * (1 - ve), pflu),
    nonflu_prob = pnonflu,
    status = calc_status(nsam, flu_prob, nonflu_prob),
    symptom_prob = if_else(status == "Healthy", 0, sympt),
    symptom = rbern(nsam, symptom_prob),
    clinic_prob = if_else(symptom == 0, 0, clin),
    clinic = rbern(nsam, clinic_prob),
    tested_prob = if_else(clinic == 0, test_nonclin, test_clin),
    tested = rbern(nsam, tested_prob),
    test_result_prob = case_when(
      tested == 0 ~ NA_real_,
      status == "Flu" ~ sens_flu,
      TRUE ~ 1 - spec_flu
    ),
    test_result = suppressWarnings(rbern(nsam, test_result_prob))
  )
  attr(pop, "true_vals") <- tibble(
    nsam, pvac, pflu, ve, pnonflu, sympt, clin,
    test_clin, test_nonclin, sens_vac, spec_vac, sens_flu, spec_flu, seed
  )
  pop
}

# Count cases/controls
summarise_pop <- function(pop, type = "surv") {
  type <- arg_match(type, c("surv", "admin"))
  if (type == "surv") pop <- pop %>%
      filter(clinic == 1, !is.na(test_result)) %>%
      mutate(case = test_result == 1)
  else pop <- pop %>%
      filter(!is.na(test_result)) %>%
      mutate(case = test_result == 1)
  summ <- pop %>%
    count(case, vac_record) %>%
    mutate(type = type)
  attr(summ, "true_vals") <- attr(pop, "true_vals")
  summ
}

# Calculate ve from case/control counts
calc_ve <- function(n, case, vac_record) {
  if (length(n) < 4) {
    warn("cannot calculate ve")
    return(NA_real_)
  }
  cv <- n[case & vac_record]
  cu <- n[case & !vac_record]
  nv <- n[!case & vac_record]
  nu <- n[!case & !vac_record]
  1  - (cv / cu) / (nv / nu)
}

# Simulate a TN study (both surveillance and administrative data)
simulate_study <- function(..., name = "default") {
  pop <- simulate_pop(...)
  study <- map_dfr(c("surv", "admin"), ~ summarise_pop(pop, .x)) %>%
    mutate(name = name)
  attr(study, "true_vals") <- attr(pop, "true_vals") %>% mutate(name = name)
  study
}

# Simualate a TN study from a parameter dictionary
# Assumes that nms subset from pars_dict are all in one population
simulate_study_pd <- function(nsam, nms, pars_dict,
                              init_seed = sample.int(.Machine$integer.max, 1)) {
  pars <- pars_dict %>%
    filter(name %in% nms) %>%
    mutate(
      prop = prop / sum(prop),
      nsam = round(nsam * prop, 0),
      seed = init_seed + row_number() - 1
    )
  if (nrow(pars) != length(nms)) abort("one row per name in pars_dict")
  true_vals <- pars
  if (length(nms) > 1) {
    pars_total <- pars %>%
      select(-name, -prop, -nsam, -seed) %>%
      summarise_all(function(vec) sum(vec * pars$prop)) %>%
      mutate(name = "total", prop = 1, nsam = nsam, seed = NA_integer_)
    true_vals <- bind_rows(true_vals, pars_total)
  }
  study <- pmap_dfr(select(pars, -prop), function(..., name) {
    do.call(simulate_study, c(list(...), list(name = name)))
  })
  if (length(nms) > 1) {
    total <- study %>%
      group_by(case, vac_record, type) %>%
      summarise(n = sum(n)) %>%
      ungroup() %>%
      mutate(name = "total")
    study <- bind_rows(study, total)
  }
  attr(study, "true_vals") <- true_vals
  study
}

# Calculate VE, attach true values
summarise_study <- function(study) {
  study %>%
    group_by(name, type) %>%
    summarise(ve_est = calc_ve(n, case, vac_record), n_study = sum(n)) %>%
    ungroup() %>%
    left_join(attr(study, "true_vals"), by = "name")
}

# Simulate and summarise many studies
many_studies <- function(nsim, nsam, nms, pars_dict,
                         init_seed = sample.int(.Machine$integer.max, 1)) {
  per_sim <- length(nms)
  future_map_dfr(1:nsim, function(ind) {
    simulate_study_pd(nsam, nms, pars_dict, init_seed + (ind - 1) * per_sim) %>%
      summarise_study() %>%
      mutate(index = ind)
  })
}

# Create a table with all possible combinations from var_list
create_all_combos <- function(var_list, pars_dict) {
  add_combos <- function(curr, add_vals, add_name) {
    curr_rows <- nrow(curr)
    curr %>%
      slice(rep(1:n(), each = length(add_vals))) %>%
      mutate(!!sym(add_name) := rep(add_vals, times = curr_rows))
  }
  reduce2(
    var_list, names(var_list), add_combos, .init = tibble(.rows = 1)
  )
}

# Vary (multiple) parameters in a population with one group in it
vary_pars <- function(var_names, vary_table, nsim, nsam, set_name, pars_dict,
                     init_seed = sample.int(.Machine$integer.max, 1)) {
  if (length(set_name) > 1) abort("population has to consist of one group")
  if (any(!var_names %in% names(pars_dict)))
    abort("some var_names not in pars_dict")
  needed_vars <- vary_table[names(vary_table) %in% var_names] %>%
    create_all_combos()
  pars <- pars_dict %>%
    mutate(prop = 1, clin = 1) %>% # Only have meaning in multi-group context
    filter(name == set_name) %>%
    select(-var_names) %>%
    group_split(name) %>%
    map_dfr(function(.tbl) {
      if (nrow(.tbl) > 1) abort("one row per name in pars_dict")
      .tbl %>% slice(rep(1, nrow(needed_vars))) %>% bind_cols(needed_vars)
    })
  map_dfr(
    1:nrow(pars), function(ind) many_studies(
      nsim, nsam, set_name, slice(pars, ind), init_seed + (ind - 1) * nsim
    )
  ) %>% mutate(vary_name = paste(var_names, collapse = "-"))
}

# Create parameter combinations for multi-group variation
create_mult <- function(vary_list, nms, pars_dict) {
  if (length(nms) != 3) abort("nms should be of length 3")
  pars <- pars_dict %>% filter(name %in% nms)
  mult <- tibble()
  for (var_name in names(vary_list)) {
    if (!var_name %in% names(pars_dict))
      abort(paste0(var_name, " not in pars_dict"))
    vals <- vary_list[[var_name]]
    if (length(vals) != 3) abort("values should be of length 3")
    entry <- pars
    entry$vary_par <- var_name
    eqs <- tibble()
    lbls <- c("low", "med", "high")
    for (i in seq_along(lbls)) {
      all_eq <- mutate(
        entry, !!sym(var_name) := vals[[i]],
        vary_type = paste0("all-", lbls[[i]])
      )
      eqs <- bind_rows(eqs, all_eq)
    }
    lows <- tibble()
    for (i in 1:3) {
      one_low <- entry %>%
        mutate(
          !!sym(var_name) := if_else(name == nms[[i]], vals[[1]], vals[[3]]),
          vary_type = paste0(nms[[i]], "-low")
        )
      lows <- bind_rows(lows, one_low)
    }
    entry <- bind_rows(all_eq, lows)
    mult <- bind_rows(mult, entry)
  }
  mult
}

# Vary parameters in a multi-group context
vary_mult <- function(var_names, vary_table, nsim, nsam, set_names, pars_dict,
                      init_seed = sample.int(.Machine$integer.max, 1)) {
  pars_dicts <- create_mult(vary_table, set_names, pars_dict)
  pars_dicts <- filter(pars_dicts, vary_par %in% var_names)
  nms <- unique(pars_dicts$name)
  pars_split <- group_split(pars_dicts, vary_type, vary_par)
  imap_dfr(
    pars_split, function(pars, ind) {
      vary_type <- unique(pars$vary_type)
      pars$vary_type <- NULL
      vary_par <- unique(pars$vary_par)
      pars$vary_par <- NULL
      many_studies(
        nsim, nsam, nms, pars, init_seed + (ind - 1) * nsim * length(nms)
      ) %>% mutate(vary_name = paste0(vary_par, "--", vary_type))
    }
  )
}

# Vary one parameter at a time in a population with one group in it
vary_pars_1aat <- function(par_names, set_names, nsim, nsam,
                           vary_table, pars_dict,
                           init_seed = sample.int(.Machine$integer.max, 1)) {
  sims_per_par <- map_dbl(vary_table[par_names], function(.x) length(.x) * nsim)
  sims_per_set <- sum(sims_per_par)
  imap_dfr(
    set_names,
    function(set_name, set_ind) {
      imap_dfr(
        par_names,
        function(par_name, par_ind) {
          if (par_ind == 1) sims_per_prev_par <- 0
          else sims_per_prev_par <- sims_per_par[[par_ind - 1]]
          vary_pars(
            par_name, vary_table, nsim, nsam, set_name, pars_dict,
            init_seed = init_seed + sims_per_prev_par * (par_ind - 1) +
              sims_per_set * (set_ind - 1)
          )
        }
      )
    }
  )
}

# Vary multiple parameters at a time in a population with one group in it
vary_pars_maat <- function(par_names, set_names, nsim, nsam,
                           vary_table, pars_dict,
                           init_seed = sample.int(.Machine$integer.max, 1)) {
  sims_per_set <- nsim * nrow(
    create_all_combos(vary_table[names(vary_table) %in% par_names], pars_dict)
  )
  imap_dfr(
    set_names,
    function(set_name, set_ind) {
      vary_pars(
        par_names, vary_table, nsim, nsam, set_name, pars_dict,
        init_seed = init_seed + sims_per_set * (set_ind - 1)
      )
    }
  )
}

# Save results
save_res <- function(res, name, nsim, folder) {
  write_csv(
    res, file.path(folder, paste0(name, "-", nsim, "sims.csv"))
  )
}

# Script ======================================================================
# Uncomment to regenerate results

# Multiple groups

# sims_mult <- map_dfr(
#   names(vary_table_mult),
#   vary_mult,
#   vary_table_mult, nsim, 5e5, c("children", "adults", "elderly"),
#   pars_dict, 20191118
# )
# save_res(sims_mult, "mult", nsim, sim_dir)

# Single group, single parameter

# sims_one <- vary_pars_1aat(
#   names(vary_table), pars_dict$name,
#   nsim, 5e5, vary_table, pars_dict, 20191118
# )
# save_res(sims_one, "one", nsim, sim_dir)

# VE investigation

# sims_ve <- vary_pars_maat(
#   names(vary_table_ve),
#   c(
#     "special_no", "special_se_f",
#     "special_sp_f", "special_se_v", "special_sp_v"
#   ),
#   nsim, 5e5, vary_table, pars_dict, 20191118
# )
# save_res(sims_ve, "veinv", nsim, sim_dir)

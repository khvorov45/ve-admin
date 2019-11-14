# Simulations of VE studies
# Arseniy Khvorov
# Created 2019/11/13
# Last edit 2019/11/13

library(tidyverse)
library(rlang)
library(furrr) # Like purrr but with parallel support
library(extraDistr) # For the categorical and bernouilli disributions

plan(multiprocess) # Parallel setup

sim_dir <- "sim"
pars_dir <- "pars"

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
simulate_study_pd <- function(nsam, nms, pars_dict,
                              init_seed = sample.int(.Machine$integer.max, 1)) {
  pars <- pars_dict %>%
    filter(name %in% nms) %>%
    mutate(
      prop = prop / sum(prop),
      nsam = round(nsam * prop, 0),
      seed = init_seed + row_number() - 1
    )
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

summarise_study <- function(study) {
  study %>%
    group_by(name, type) %>%
    summarise(ve_est = calc_ve(n, case, vac_record),) %>%
    ungroup() %>%
    left_join(attr(study, "true_vals"), by = "name")
}

# Script ======================================================================

stud <- simulate_study_pd(1e5, c("children", "adults", "elderly"), pars_dict)
stud <- simulate_study()
summarise_study(stud)


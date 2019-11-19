# Tests of simulations of VE studies
# Arseniy Khvorov
# Created 2019/11/13
# Last edit 2019/11/13

library(testthat)

test_that("simulation function works", {

  # Sample size
  pop <- simulate_pop(
    nsam = 30, pvac = 0.5, pflu = 0.1, ve = 0.5, pnonflu = 0.1, sympt = 0.9,
    clin = 0.8, test_clin = 0.7, test_nonclin = 0.6,
    sens_vac = 0.6, spec_vac = 0.7, sens_flu = 0.8, spec_flu = 0.9,
  )
  expect_equal(nrow(pop), 30)

  # Vaccinated proportion
  n_large <- 1e6
  pop <- simulate_pop(nsam = n_large, pvac = 0.6)
  expect_equal(sum(pop$vaccinated) / nrow(pop), 0.6, tol = 0.001)

  # Infection status (and VE)
  pop <- simulate_pop(n_large, pvac = 0.5, pflu = 0.2, pnonflu = 0.2, ve = 0.5)
  pop_sum_temp <- pop %>%
    mutate(case = status == "Flu") %>%
    count(vaccinated, case, name = "n_cell") %>%
    mutate(vaccinated = if_else(vaccinated == 1, "vac", "unvac")) %>%
    pivot_wider(names_from = vaccinated, values_from = n_cell)
  r_vac <- pop_sum_temp$vac[pop_sum_temp$case] / sum(pop_sum_temp$vac)
  r_unvac <- pop_sum_temp$unvac[pop_sum_temp$case] / sum(pop_sum_temp$unvac)
  expect_equal(r_vac, 0.2 * 0.5, tol = 0.001)
  expect_equal(r_unvac, 0.2, tol = 0.001)
  expect_equal(sum(pop$status == "Nonflu") / nrow(pop), 0.2, tol = 0.001)

  # Vaccination record
  pop <- simulate_pop(n_large, pvac = 0.5, sens_vac = 0.8, spec_vac = 0.9)
  exp_rec <- sum(pop$vaccinated) * 0.8 + sum(!pop$vaccinated) * (1 - 0.9)
  expect_equal(exp_rec / n_large, sum(pop$vac_record) / n_large, tol = 0.001)

  # Flu test
  pop <- simulate_pop(n_large, pflu = 0.2, sens_flu = 0.8, spec_flu = 0.9)
  exp_pos <- sum(pop$status == "Flu") * 0.8 +
    sum(pop$status != "Flu") * (1 - 0.9)
  expect_equal(exp_pos / n_large, sum(pop$test_result) / n_large, tol = 0.001)

  # Symptomatic
  pop <- simulate_pop(n_large, sympt = 0.5)
  expect_equal(
    sum(pop$status != "Healthy") / n_large * 0.5,
    sum(pop$symptom) / n_large,
    tol = 0.001
  )

  # Clinical
  pop <- simulate_pop(n_large, clin = 0.5)
  expect_equal(
    sum(pop$status != "Healthy") / n_large * 0.5,
    sum(pop$clinic) / n_large,
    tol = 0.001
  )

  # Tested probability
  pop <- simulate_pop(n_large, test_clin = 0.5, test_nonclin = 0.2)
  expect_equal(
    sum(pop$tested & pop$status != "Healthy") / sum(pop$status != "Healthy"),
    0.5, tol = 0.001
  )
  expect_equal(
    sum(pop$tested & pop$status == "Healthy") / sum(pop$status == "Healthy"),
    0.2, tol = 0.001
  )
})

test_that("population summary works", {
  pop <- simulate_pop(1e3, pflu = 0.3, pnonflu = 0.3)
  summ_surv <- summarise_pop(pop, type = "surv")
  expect_equal(
    summ_surv$n[summ_surv$case & summ_surv$vac_record],
    nrow(filter(pop, vac_record == 1, test_result == 1, clinic == 1))
  )
  expect_equal(
    summ_surv$n[summ_surv$case & !summ_surv$vac_record],
    nrow(filter(pop, vac_record == 0, test_result == 1, clinic == 1))
  )
  expect_equal(
    summ_surv$n[!summ_surv$case & summ_surv$vac_record],
    nrow(filter(pop, vac_record == 1, test_result == 0, clinic == 1))
  )
  expect_equal(
    summ_surv$n[!summ_surv$case & !summ_surv$vac_record],
    nrow(filter(pop, vac_record == 0, test_result == 0, clinic == 1))
  )
  summ_admin <- summarise_pop(pop, type = "admin")
  expect_equal(
    summ_admin$n[summ_admin$case & summ_admin$vac_record],
    nrow(filter(pop, vac_record == 1, test_result == 1))
  )
  expect_equal(
    summ_admin$n[summ_admin$case & !summ_admin$vac_record],
    nrow(filter(pop, vac_record == 0, test_result == 1))
  )
  expect_equal(
    summ_admin$n[!summ_admin$case & summ_admin$vac_record],
    nrow(filter(pop, vac_record == 1, test_result == 0))
  )
  expect_equal(
    summ_admin$n[!summ_admin$case & !summ_admin$vac_record],
    nrow(filter(pop, vac_record == 0, test_result == 0))
  )
})

test_that("ve caclulations work", {
  pop <- simulate_pop(2e6, pflu = 0.3, pnonflu = 0.3, ve = 0.4)
  pop_sum_surv <- summarise_pop(pop, type = "surv")
  expect_equal(
    calc_ve(pop_sum_surv$n, pop_sum_surv$case, pop_sum_surv$vac_record),
    0.4, tol = 0.001
  )
})

test_that("study simulation and summary work", {
  stud_one <- simulate_study()
  stud_one_summ <- summarise_study(stud_one)
  stud_many <- simulate_study_pd(
    1e5, c("children", "adults", "elderly"), pars_dict
  )
  stud_many_summ <- summarise_study(stud_many)
  expect_true("prop" %in% names(stud_many_summ))
  expect_equal(
    sort(names(stud_many_summ)[names(stud_many_summ) != "prop"]),
    sort(names(stud_one_summ))
  )
})

test_that("many studies are simulated", {
  studs <- many_studies(
    10, 1e5, c("children", "adults", "elderly"), pars_dict, init_seed = 1
  )
  expect_true(all(studs$seed[!is.na(studs$seed)] %in% 1:30))
})

test_that("parameters can be varied", {
  vary_table_light <- list(pvac = c(0.05, 0.1, 0.2), sens_vac = c(0.9, 0.95))
  expect_equal(nrow(create_all_combos(vary_table_light)), 4)
  sims <- vary_pars(
    c("pvac", "sens_vac"), vary_table_light, 5, 1e5, "children", pars_dict,
    init_seed = 1
  )
  expect_equal(
    sims$pvac %>% sort() %>% unique(), c(0.05, 0.1)
  )
  expect_equal(
    sims$sens_vac %>% sort() %>% unique(), c(0.9, 0.95)
  )
  expect_equal(
    sims$seed %>% sort() %>% unique(), 1:20
  )
  vary_table_mult_light <- list(
    ve = c(0.15, 0.33, 0.7),
    pvac = c(0.05, 0.3, 0.5)
  )
  sims_mult <- vary_mult(
    "pvac", vary_table_mult_light, 5, 1e5, c("children", "adults", "elderly"),
    pars_dict, 1
  )
  expect_equal(sort(names(sims_mult)), sort(names(sims)))
})

test_that("one parameter at a time in one group works", {
  vary_table_light <- list(pvac = c(0.05, 0.1, 0.2), sens_vac = c(0.9, 0.95))
  sims <- vary_pars_1aat(
    c("pvac", "sens_vac"), c("children", "adults"),
    5, 1e5, vary_table_light, pars_dict,
    init_seed = 10
  )
  expect_equal(length(unique(sims$seed)), 2 * (3 * 5 + 2 * 5))
})

test_that("multiple parameters at a time in one group works", {
  vary_table_light <- list(pvac = c(0.05, 0.1, 0.2), sens_vac = c(0.9, 0.95))
  sims <- vary_pars_maat(
    c("pvac", "sens_vac"), c("children", "adults"),
    5, 1e5, vary_table_light, pars_dict,
    init_seed = 10
  )
  expect_equal(length(unique(sims$seed)), 2 * (6 * 5))
})


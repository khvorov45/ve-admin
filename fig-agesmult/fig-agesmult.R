# Graphs of summaries of multiple group simulations
# Arseniy Khvorov
# Created 2019/11/21
# Last edit 2019/11/21

library(tools)
library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(ggrepel)

# Directory names that will be used later
fig_agesmult_dir <- "fig-agesmult"
sim_summ_dir <- "sim-summary"

# Settings ====================================================================

par_lbls <- list(
  pvac = bquote(italic(v)),
  pflu = bquote(italic(v)),
  ve = bquote(italic(e)),
  pnonflu = bquote(italic(l)),
  sympt = bquote(italic(p)),
  clin = bquote(italic(c)), # c
  test_clin = bquote(italic(t)[italic(a)]),
  test_nonclin = bquote(italic(t)[italic(n)]),
  sens_vac = bquote(italic(s)[italic(e)][","][italic(v)]),
  spec_vac = bquote(italic(s)[italic(p)][","][italic(v)]),
  sens_flu = bquote(italic(s)[italic(e)]),
  spec_flu = bquote(italic(s)[italic(e)])
)

agesmult_theme <- theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(0, "null"),
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    legend.margin = margin(0, 0, 0, 0)
  )

# Functions ===================================================================

plot_mult <- function(mult_onepar, par_lbls) {
  true_line <- function(vary_name) {
    if (vary_name == "ve") geom_abline(slope = 1, intercept = 0, lty = "1111")
    else geom_hline(aes(yintercept = ve), lty = "1111")
  }
  vary_name <- unique(mult_onepar$vary_par)
  ind <- mult_onepar %>%
    filter(name != "Total") %>%
    group_by(name, !!sym(vary_name), type, ve) %>%
    summarise(ve_est_mean := mean(ve_est_mean))
  total <- filter(mult_onepar, name == "Total")
  mod_dat <- bind_rows(total, ind)
  pl <- mod_dat %>%
    ggplot(aes(!!sym(vary_name), ve_est_mean)) +
    dark_mode(agesmult_theme) +
    facet_wrap(~name, nrow = 1) +
    scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
    scale_linetype_discrete("Data source") +
    scale_shape_manual("Data source", values = c(17, 18)) +
    xlab(par_lbls[[vary_name]]) +
    ylab("Estimated VE") +
    true_line(vary_name) +
    geom_hline(yintercept = 0, lty = "3111") +
    geom_line(aes(lty = type)) +
    geom_point(aes(shape = type)) +
    geom_text_repel(
      data = filter(mod_dat, type == "Surveillance", name == "Total"),
      mapping = aes(label = vary_patt_ind), col = "gray50"
    )
  attr(pl, "plotname") <- vary_name
  pl
}

save_agesmult <- function(pl, folder) {
  ggsave_dark(
    file.path(folder, paste0("agesmult-", attr(pl, "plotname"), ".pdf")),
    pl, dark = FALSE,
    width = 15, height = 7, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

agesmult <- read_csv(file.path(sim_summ_dir, "mult-500sims.csv")) %>%
  mutate(
    name = toTitleCase(name) %>%
      factor(levels = c("Children", "Adults", "Elderly", "Total")),
    type = recode(type, "admin" = "Administrative", "surv" = "Surveillance"),
    vary_par = str_extract(vary_name, ".*--") %>% str_replace("--", ""),
    vary_patt = str_extract(vary_name, "--.*") %>% str_replace("--", ""),
    vary_patt_ind = recode(
      vary_patt, "all-low" = "l", "all-med" = "m", "all-high" = "h",
      "children-high" = "ch", "adults-high" = "ah", "elderly-high" = "eh"
    )
  )

mult_pls <- agesmult %>%
  group_split(vary_par) %>%
  map(plot_mult, par_lbls)

walk(mult_pls, save_agesmult, fig_agesmult_dir)

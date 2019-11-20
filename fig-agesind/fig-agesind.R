# Figures of summaries of individual age group group simulations
# Arseniy Khvorov
# Created 2019/11/20
# Last edit 2019/11/21

library(tools)
library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directory names that will be used later
fig_dir <- "fig-agesind"
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

onevar_plot_theme <- theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(0, "null"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90),
    panel.grid.minor = element_blank(),
    legend.box.spacing = unit(0, "null"),
    legend.margin = margin(0, 0, 0, 0)
  )

# Functions ===================================================================

# Plots variation of one parameter across multiple groups
onevar_plot <- function(onevar_onepar, par_lbls, theme) {
  vary_name <- unique(onevar_onepar$vary_name)
  xlbls <- unique(onevar_onepar[[vary_name]])
  pl <- onevar_onepar %>%
    ggplot(aes(!!sym(vary_name), ve_est_mean)) +
    dark_mode(theme) +
    scale_x_continuous(breaks = xlbls) +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    scale_linetype_discrete("Data source") +
    scale_shape_manual("Data source", values = c(17, 18)) +
    xlab(par_lbls[[vary_name]]) +
    ylab("Estimated VE") +
    facet_wrap(~name, nrow = 1) +
    geom_hline(aes(yintercept = ve), lty = "1111") +
    geom_point(aes(shape = type)) +
    geom_line(aes(lty = type))
  attr(pl, "plotname") <- vary_name
  pl
}

save_plot <- function(pl, folder) {
  ggsave_dark(
    file.path(folder, paste0("agesind-", attr(pl, "plotname"), ".pdf")),
    pl, dark = FALSE,
    width = 15, height = 6.5, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

# Age-representative groups, one parameter varied
onevar <- read_csv(file.path(sim_summ_dir, "one-500sims.csv"))
ages  <- filter(onevar, name %in% c("children", "adults", "elderly"))

ages_plots <- ages %>%
  mutate(
    name = toTitleCase(name) %>%
      factor(levels = c("Children", "Adults", "Elderly")),
    type = recode(type, "admin" = "Administrative", "surv" = "Surveillance")
  ) %>%
  group_split(vary_name) %>%
  map(onevar_plot, par_lbls, onevar_plot_theme)

walk(ages_plots, save_plot, fig_dir)


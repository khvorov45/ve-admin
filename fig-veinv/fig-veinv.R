# Figures of VE investion summary
# Arseniy Khvorov
# Created 2019/11/21
# Last edit 2019/11/21

library(tools)
library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(latex2exp)

# Directory names that will be used later
fig_veinv_dir <- "fig-veinv"
sim_summ_dir <- "sim-summary"

# Settings ====================================================================

veinv_theme <- theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(0, "null"),
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    legend.margin = margin(0, 0, 0, 0)
  )

# Functions ===================================================================

veinv_plot <- function(veinv, theme) {
  label_facets <- function(lbls) {
    if (names(lbls) == "pvac")
      lbls$pvac <- paste0("$\\textit{v} = ", lbls$pvac, "$")
    lbls[[1]] <- TeX(lbls[[1]])
    label_parsed(lbls)
  }
  veinv %>%
    ggplot(aes(ve, ve_est_mean)) +
    dark_mode(theme) +
    facet_grid(name_lbls ~ pvac, labeller = label_facets) +
    scale_x_continuous(breaks = seq(0.2, 0.8, 0.2)) +
    scale_y_continuous(breaks = seq(0.2, 0.8, 0.2)) +
    scale_linetype_discrete("Data source") +
    scale_shape_manual("Data source", values = c(17, 18)) +
    xlab(TeX("\\textit{e}")) +
    ylab("Estimated VE") +
    geom_abline(slope = 1, intercept = 0, lty = "1111") +
    geom_point(aes(shape = type)) +
    geom_line(aes(lty = type))
}

save_veinv <- function(pl, n_hor, n_vertical, name, folder) {
  ggsave_dark(
    file.path(folder, paste0("veinv-", name, ".pdf")),
    pl, dark = FALSE,
    width = 2.2 * n_hor, height = 3.5 * n_vertical, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

veinv <- read_csv(file.path(sim_summ_dir, "veinv-500sims.csv")) %>%
  mutate(
    name_lbls = recode(
      name,
      "special_no" = "No mis.",
      "special_se_f" = "$\\textit{s}_e = 0.95$",
      "special_sp_f" = "$\\textit{s}_p = 0.95$",
      "special_se_v" = "$\\textit{s}_{e,v} = 0.95$",
      "special_sp_v" = "$\\textit{s}_{p,v} = 0.95$"
    ),
    type = recode(type, "admin" = "Administrative", "surv" = "Surveillance")
  )

full <- veinv_plot(veinv, veinv_theme)
save_veinv(full, 10, 5, "full", fig_veinv_dir)

limited <- veinv %>%
  filter(
    name %in% c("special_no", "special_sp_v"),
    round(pvac, 1) %in% c(0.1, 0.3, 0.5, 0.7, 0.9)
  ) %>%
  veinv_plot(veinv_theme)
save_veinv(limited, 5, 2, "limited", fig_veinv_dir)


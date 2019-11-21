# Explanation of the effect of spec_vac
# Arseniy Khvorov
# Created 2019/04/10
# Last edit 2019/11/21

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(latex2exp)

# Directory names that will be used later
fig_sfinv_dir <- "fig-sfinv"

# Settings ====================================================================

sfinv_theme <- theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.spacing = unit(0, "null"),
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    legend.margin = margin(0, 0, 0, 0)
  )

# Functions ===================================================================

ORs <- function(f, e, l, s) (f * (1 - e) + l * (1 - s)) / (f + l * (1 - s))

ORa <- function(f, e, s) (
  ((1 - s * (1 - f * (1 - e))) * (s * (1 - f))) /
    ((1 - s * (1 - f)) * (s * (1 - f * (1 - e))))
)

label_facets_sfinv <- function(lbls) {
  lbls$varied <- recode(
    lbls$varied, "f" = "$\\textit{f}$", "Sp" = "$\\textit{s}_p$"
  ) %>% TeX()
  label_parsed(lbls)
}

# Script ======================================================================

df_f <- tibble(
  f = seq(0.05, 0.15, 10^-4),
  e = 0.5,
  l = 0.15,
  s = 0.95,
  ORad = ORa(f, e, s),
  ORsu = ORs(f, e, l, s),
  ORtrue = 1 - e,
  varied = "f"
)

df_s <- tibble(
  f = 0.1,
  e = 0.5,
  l = 0.15,
  s = seq(0.9, 1, 10^-4),
  ORad = ORa(f, e, s),
  ORsu = ORs(f, e, l, s),
  ORtrue = 1 - e,
  varied = "Sp"
)

dfg <- bind_rows(df_f, df_s) %>%
  pivot_longer(
    c(ORad, ORsu), names_to = "ORtype", values_to = "OR"
  ) %>%
  mutate(
    ORtype = recode(ORtype, "ORad" = "Administrative", "ORsu" = "Surveillance")
  )

dfg %>%
  ggplot(aes(y = OR, lty = ORtype)) +
  dark_mode(sfinv_theme) +
  geom_hline(aes(yintercept = ORtrue), linetype = "1111") +
  geom_line(data = subset(dfg, varied == "f"), aes(x = f), lwd = 1) +
  geom_line(data = subset(dfg, varied == "Sp"), aes(x = s), lwd = 1) +
  facet_wrap(vars(varied), scales = "free_x", labeller = label_facets_sfinv) +
  scale_linetype_discrete(name = "Data Source") +
  xlab("Value")

ggsave_dark(
  filename = file.path(fig_sfinv_dir, "sfinv.pdf"), dark = FALSE,
  width = 12, height = 6,
  units = "cm", dpi = 500
)


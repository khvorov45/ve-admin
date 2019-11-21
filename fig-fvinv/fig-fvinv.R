# Effects of IP_flu and p_vac
# Arseniy Khvorov
# Created 2019/03/26
# Last edit 2019/11/22

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directory names that will be used later
fig_fvinv_dir <- "fig-fvinv"

# Functions ===================================================================

# Overall OR in terms of p_vac and IP_flu in the first group
ORov <- function(v1, f1) {
  (1 - e) * ((v1 * f1 + 2 * v * f) * (3 - v1 - 2 * v)) /
    (((1 - v1) * f1 + 2 * (1 - v) * f) * (v1 + 2 * v))
}

# Script ======================================================================

e <- 0.5 # True VE
v <- 0.5 # p_vac in the other 2 groups
f <- 0.1 # IP_flu in the other 2 groups

v1 <- rep(seq(0.1, 0.9, length.out = 100), 100)
f1 <- rep(seq(0.01, 0.3, length.out = 100), 100)

df <- tibble(
  f1 = sort(f1), v1, VE = 1 - ORov(v1, f1), Bias = VE - e,
  RelFluProp = f1 / f, RelVacProp = v1 / v
)
df$VE <- 1 - fun(df$v1, df$f1)
df$Bias <- df$VE - e
df$RelFluProp <- df$f1 / f
df$RelVacProp <- df$v1 / v

ggplot(df, aes(RelFluProp, RelVacProp)) +
  dark_theme_bw(verbose = FALSE) +
  geom_tile(aes(col = Bias, fill = Bias)) +
  scale_color_gradient2(low = "blue", high = "red") +
  scale_fill_gradient2(low = "blue", high = "red") +
  geom_hline(yintercept = 1, color = "gray") +
  geom_vline(xintercept = 1, color = "gray") +
  ylab("Relative vaccinated proportion") + xlab("Relative flu incidence")
ggsave_dark(
  file.path(fig_fvinv_dir, "fig-fvinv.pdf"), dark = FALSE,
  width = 16, height = 10, units = "cm"
)

# Graph of how a subgroup may impact the overall estimate
# Arseniy Khvorov
# Created 2019/03/20
# Last edit 2019/11/21

library(ggplot2)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directory names that will be used later
fig_ortn_dir <- "fig-ortn"

# Functions ===================================================================

ORbad <- function(v1) {
  (
    fv * (
      n1 * v1 + n2 * v2 + n3 * v3
    ) * (
      n1 * (1 - v1) * hu + l * (
        n1 * (1 - v1) + n2 * (1 - v2) + n3 * (1 - v3)
      )
    )
  ) /
  (
    fu * (
      n1 * (1 - v1) + n2 * (1 - v2) + n3 * (1 - v3)
    ) * (
      n1 * v1 * hv + l * (
        n1 * v1 + n2 * v2 + n3 * v3
      )
    )
  )
}

# Script ======================================================================

v1 <- seq(0.1, 0.9, 0.0001)

v2 <- 0.5
v3 <- 0.5

n1 <- 1000
n2 <- 1000
n3 <- 1000

fv <- 0.025
fu <- 0.05

l <- 0.1

hv <- 1 - l - fv
hu <- 1 - l - fu

df <- data.frame(v1, ORbad = ORbad(v1))

pl <- ggplot(df, aes(v1, ORbad)) +
  dark_theme_bw(verbose = FALSE) +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = fv / fu, linetype = "3333") +
  geom_hline(yintercept = 1, linetype = "3111") +
  geom_line() +
  scale_x_continuous(
    name = "Vaccinated proportion in subgroup 1", breaks = seq(0,1,0.1)
  ) +
  scale_y_continuous(
    name = "Overall population OR", breaks = seq(0, 2, 0.1)
  )
pl
ggsave_dark(
  file.path(fig_ortn_dir, "fig-ortn.pdf"), pl, dark = FALSE,
  width = 12, height = 7, units = "cm"
)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gpindex)

# Start with some data on prices and quantities for 6 products
# over 5 periods
price6
quantity6

# We'll only need prices and quantities for a few periods
p0 <- price6[[1]]
p1 <- price6[[2]]
p2 <- price6[[3]]
q0 <- price6[[1]]
q1 <- price6[[2]]

s0 <- p0 * q0
s1 <- p1 * q1

# Laspeyres index
arithmetic_mean(p1 / p0, s0)

## -----------------------------------------------------------------------------
# Lloyd-Moulton index (elasticity of substitution -1)
quadratic_mean <- generalized_mean(2)
quadratic_mean(p1 / p0, s0)

# Palgrave index
arithmetic_mean(p1 / p0, s1)

## -----------------------------------------------------------------------------
# Fisher index
fisher_mean(p1 / p0, s0, s1)

# Drobisch index
drobisch_mean <- nested_mean(1, c(1, 1))
drobisch_mean(p1 / p0, s0, s1)

# Geometric AG mean index (elasticity of substitution 0.25)
ag_mean <- nested_mean(0, c(0, 1), c(0.25, 0.75))
ag_mean(p1 / p0, s0, s0)

## -----------------------------------------------------------------------------
# Laspeyres index, again
arithmetic_mean(p1 / p0, index_weights("Laspeyres")(p0, q0))

laspeyres_index(p1, p0, q0)

## -----------------------------------------------------------------------------
quadratic_decomposition <- transmute_weights(2, 1)

arithmetic_mean(p1 / p0, quadratic_decomposition(p1 / p0, s0))
quadratic_mean(p1 / p0, s0)

quadratic_contributions <- contributions(2)
quadratic_contributions(p1 / p0, s0)

quadratic_update <- factor_weights(2)

quadratic_mean(p2 / p0, s0)
quadratic_mean(p2 / p1, quadratic_update(p1 / p0, s0)) *
  quadratic_mean(p1 / p0, s0)

## -----------------------------------------------------------------------------
ag_decomposition <- nested_transmute(0, c(0, 1), 1, c(0.25, 0.75))

ag_mean(p1 / p0, s0, s0)
arithmetic_mean(p1/ p0, ag_decomposition(p1 / p0, s0, s0))

ag_contributions <- nested_contributions(0, c(0, 1), c(0.25, 0.75))
ag_contributions(p1 / p0, s0, s0)


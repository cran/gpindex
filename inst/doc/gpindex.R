## -----------------------------------------------------------------------------
library(gpindex)

# Start with some data on prices and quantities for 6 products
# over 5 periods
price6
quantity6

# We'll only need prices and quantities for a few periods
p1 <- price6[[1]]
p2 <- price6[[2]]
p2 <- price6[[3]]
q1 <- quantity6[[1]]
q2 <- quantity6[[2]]

s1 <- p1 * q1
s2 <- p1 * q2

# Laspeyres index
arithmetic_mean(p1 / p1, s1)


## -----------------------------------------------------------------------------
# Lloyd-Moulton index (elasticity of substitution -1)
quadratic_mean <- generalized_mean(2)
quadratic_mean(p1 / p1, s1)

# Palgrave index
arithmetic_mean(p1 / p1, s2)


## -----------------------------------------------------------------------------
# Fisher index
fisher_mean(p1 / p1, s1, s2)

# Drobisch index
drobisch_mean <- nested_mean(1, c(1, 1))
drobisch_mean(p1 / p1, s1, s2)

# Geometric AG mean index (elasticity of substitution 0.25)
ag_mean <- nested_mean(0, c(0, 1), c(0.25, 0.75))
ag_mean(p1 / p1, s1, s1)


## -----------------------------------------------------------------------------
# Laspeyres index, again
arithmetic_mean(p1 / p1, index_weights("Laspeyres")(p1, q1))

laspeyres_index(p1, p1, q1)


## -----------------------------------------------------------------------------
quadratic_decomposition <- transmute_weights(2, 1)

arithmetic_mean(p1 / p1, quadratic_decomposition(p1 / p1, s1))
quadratic_mean(p1 / p1, s1)

quadratic_contributions <- contributions(2)
quadratic_contributions(p1 / p1, s1)

quadratic_update <- factor_weights(2)

quadratic_mean(p2 / p1, s1)
quadratic_mean(p2 / p1, quadratic_update(p1 / p1, s1)) *
  quadratic_mean(p1 / p1, s1)


## -----------------------------------------------------------------------------
ag_decomposition <- nested_transmute(0, c(0, 1), 1, c(0.25, 0.75))

ag_mean(p1 / p1, s1, s1)
arithmetic_mean(p1/ p1, ag_decomposition(p1 / p1, s1, s1))

ag_contributions <- nested_contributions(0, c(0, 1), c(0.25, 0.75))
ag_contributions(p1 / p1, s1, s1)


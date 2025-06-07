## -----------------------------------------------------------------------------
library(gpindex)

p2 <- price6[[2]]
p1 <- price6[[1]]
q2 <- quantity6[[2]]
q1 <- quantity6[[1]]

rel <- p2 / p1

s1 <- scale_weights(p1 * q1)
s2 <- scale_weights(p2 * q2)

quadratic_mean <- generalized_mean(2)
quadratic_decomposition <- transmute_weights(2, 1)

v <- quadratic_decomposition(rel, s1)

all.equal(
  quadratic_mean(rel, s1),
  arithmetic_mean(rel, v)
)


## -----------------------------------------------------------------------------
(v - s1)[order(rel)]


## -----------------------------------------------------------------------------
all.equal(
  transmute_weights(1, 0)(rel, v),
  transmute_weights(2, 0)(rel, s1)
)


## -----------------------------------------------------------------------------
v1 <- nested_transmute(0, c(1, -1), 1)(rel, s1, s2)

all.equal(fisher_mean(rel, s1, s2), arithmetic_mean(rel, v1))

all.equal(
  v1,
  quadratic_decomposition(rel, nested_transmute(0, c(1, -1), 2)(rel, s1, s2))
)


## -----------------------------------------------------------------------------
v2 <- nested_transmute2(0, c(1, -1), 1)(rel, s1, s2)

all.equal(fisher_mean(rel, s1, s2), arithmetic_mean(rel, v2))

all.equal(
  v2,
  quadratic_decomposition(rel, nested_transmute2(0, c(1, -1), 2)(rel, s1, s2))
)


## -----------------------------------------------------------------------------
summary(v1 - v2)


## -----------------------------------------------------------------------------
group <- rep(c("a", "b"), each = 3)

s1_by_group <- split(s1, group)
rel_by_group <- split(rel, group)

index_a <- quadratic_mean(rel_by_group$a, s1_by_group$a)
index_b <- geometric_mean(rel_by_group$b)

quadratic_mean(c(index_a, index_b), sapply(s1_by_group, sum))

decomp_a <- quadratic_decomposition(rel_by_group$a, s1_by_group$a)
decomp_b <- transmute_weights(0, 1)(rel_by_group$b)

v <- Map(
  `*`,
  quadratic_decomposition(c(index_a, index_b), sapply(s1_by_group, sum)),
  list(decomp_a, decomp_b)
) |>
  unlist()

arithmetic_mean(rel, v)


## -----------------------------------------------------------------------------
V <- sum(p2 * q2) / sum(p1 * q1)

v <- nested_transmute2(0, c(1, -1), -1)(rel, s1, s2)

all.equal(
  arithmetic_mean(V / rel, v),
  fisher_mean(q2 / q1, s1, s2)
)

V / rel * v


## -----------------------------------------------------------------------------
# Arithmetic hybrid index
all.equal(
  arithmetic_mean(p2 / p1, p2 * q1),
  contraharmonic_mean(p2 / p1, p1 * q1)
)

# Palgrave index
all.equal(
  arithmetic_mean(p2 / p1, p2 * q2),
  contraharmonic_mean(p2 / p1, p1 * q2)
)


# -----------------------------------------------------------------------------
# simulate_data.R
# Generates a SYNTHETIC attitudinal survey dataset: four latent classes with
# different positions on three attitude dimensions, ordinal items drawn around
# those positions, realistic non-response codes and a survey weight.
#
# The point is that the class structure is KNOWN, so the methods in
# lca_blueprint.Rmd can be checked against a ground truth. No real survey
# responses are involved anywhere in this file.
# -----------------------------------------------------------------------------

simulate_survey <- function(n = 1200, seed = 2026) {
  set.seed(seed)

  # Four latent classes with different mean positions on three attitude
  # dimensions (attitude, policy support, institutional norms)
  class <- sample(1:4, n, replace = TRUE, prob = c(0.35, 0.30, 0.20, 0.15))
  centre <- rbind(
    c(0.8, 0.8, 0.8),    # broadly supportive on all three
    c(0.3, 0.6, 0.7),    # low on attitude, high on norms
    c(-0.3, 0.2, -0.4),  # mixed
    c(-0.8, -0.6, -0.8)  # broadly low on all three
  )
  z <- centre[class, ]

  # helper: draw an ordinal item from a latent position
  ordinal <- function(pos, min, max, reverse = FALSE, noise = 0.7) {
    latent <- pos + rnorm(length(pos), 0, noise)
    if (reverse) latent <- -latent
    cuts <- qnorm(seq(0, 1, length.out = max - min + 2))[-c(1, max - min + 2)]
    as.numeric(cut(latent, c(-Inf, cuts, Inf))) + min - 1
  }
  add_missing <- function(x, code, rate = 0.03) {
    x[runif(length(x)) < rate] <- code
    x
  }

  df <- data.frame(ID = seq_len(n))

  # Summated-scale items, 1-7 (no don't-know option)
  for (i in 1:3) df[[paste0("scale_item_", i)]] <- ordinal(z[, 1], 1, 7)

  # Further attitude items, 1-5; items 1, 5, 6, 8, 9 are negatively worded
  for (i in 1:9) {
    df[[paste0("attitude_", i)]] <- ordinal(z[, 1], 1, 5, reverse = i %in% c(1, 5, 6, 8, 9))
  }
  df$attitude_10 <- ordinal(z[, 1], 1, 5, reverse = TRUE)
  df$attitude_11 <- ordinal(z[, 1], 1, 5, reverse = TRUE)
  df$attitude_12 <- ordinal(z[, 1], 1, 5, reverse = TRUE)
  df$attitude_13 <- ordinal(z[, 1], 1, 5, reverse = TRUE)

  # Items with a 'don't know' option coded 6
  df$attitude_14 <- add_missing(ordinal(z[, 1], 1, 5, reverse = TRUE), 6)
  df$attitude_15 <- add_missing(ordinal(z[, 1], 1, 5), 6)
  df$policy_7   <- add_missing(ordinal(z[, 2], 1, 5), 6)
  for (i in 1:7) df[[paste0("norms_", i)]] <- add_missing(ordinal(z[, 3], 1, 5, reverse = TRUE), 6)
  df$norms_8 <- add_missing(ordinal(z[, 3], 1, 5, reverse = TRUE), 6)
  df$norms_9 <- add_missing(ordinal(z[, 3], 1, 5), 6)
  df$norms_10 <- add_missing(ordinal(z[, 3], 1, 5), 6)
  df$norms_11 <- add_missing(ordinal(z[, 3], 1, 5, reverse = TRUE), 6)
  df$norms_12 <- add_missing(ordinal(z[, 3], 1, 5), 7)

  # Policy-support items, 0-10, don't know = 999
  for (i in 1:6) df[[paste0("policy_", i)]] <- add_missing(ordinal(z[, 2], 0, 10), 999)

  # Covariates
  df$sex <- haven::labelled(sample(1:2, n, replace = TRUE),
                               c(Male = 1, Female = 2))
  df$age <- round(pmin(pmax(rnorm(n, 48 + 6 * (class == 4), 16), 16), 90))
  df$income_raw <- sample(c(1:16, 99), n, replace = TRUE)
  df$education_raw <- sample(1:10, n, replace = TRUE)
  p_group <- function(k) switch(k,
    c(.25, .30, .05, .20, .15, .05),
    c(.35, .25, .15, .08, .12, .05),
    c(.25, .15, .40, .05, .05, .10),
    c(.10, .10, .65, .02, .03, .10))
  df$vote <- vapply(class, function(k)
    sample(c(1:5, 999), 1, prob = p_group(k)), numeric(1))
  df$group_var <- sample(c(1:10, 997), n, replace = TRUE)

  # Post-stratification-style weight, mean 1
  w <- exp(rnorm(n, 0, 0.35))
  df$weight <- w / mean(w)

  df
}

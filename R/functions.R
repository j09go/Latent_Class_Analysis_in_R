# -----------------------------------------------------------------------------
# functions.R  -  helper functions for lca_blueprint.Rmd
# -----------------------------------------------------------------------------

# ---- Recoding ---------------------------------------------------------------

#' Set survey non-response codes (e.g. "don't know" = 6, 999) to NA
recode_to_na <- function(x, codes) {
  x <- as.numeric(x)
  x[x %in% codes] <- NA
  x
}

#' Reverse a bounded scale so that high values point in the same direction
#' for all items. Values outside [min, max] are set to NA with a warning
#' instead of being silently reversed into the valid range.
reverse_scale <- function(x, min = 1, max = 5) {
  x <- as.numeric(x)
  bad <- !is.na(x) & (x < min | x > max)
  if (any(bad)) {
    warning(sum(bad), " value(s) outside ", min, "-", max,
            " set to NA (check for unrecoded non-response codes)", call. = FALSE)
    x[bad] <- NA
  }
  (min + max) - x
}

#' Collapse an ordinal item into 3 categories (1 = low, 2 = medium, 3 = high)
#' using upper bounds for the low and medium categories.
collapse_3 <- function(x, low_max, mid_max, max) {
  x <- as.numeric(x)
  bad <- !is.na(x) & x > max
  if (any(bad)) {
    warning(sum(bad), " value(s) above ", max, " set to NA", call. = FALSE)
    x[bad] <- NA
  }
  dplyr::case_when(
    is.na(x)      ~ NA_real_,
    x <= low_max  ~ 1,
    x <= mid_max  ~ 2,
    TRUE          ~ 3
  )
}

# ---- Measurement checks -----------------------------------------------------

#' Internal consistency (Cronbach's alpha) and the largest pairwise
#' correlation for a block of items
reliability_check <- function(df, items, collinearity_threshold = 0.7) {
  x <- as.data.frame(lapply(df[items], as.numeric))
  alpha <- suppressMessages(suppressWarnings(psych::alpha(x)))$total$raw_alpha
  r <- cor(x, use = "pairwise.complete.obs")
  r[lower.tri(r, diag = TRUE)] <- NA
  idx <- which(abs(r) > collinearity_threshold, arr.ind = TRUE)
  high <- if (nrow(idx) == 0) "none" else
    paste0(rownames(r)[idx[, 1]], " ~ ", colnames(r)[idx[, 2]],
           " (r = ", round(r[idx], 2), ")", collapse = "; ")
  tibble::tibble(
    n_items = length(items),
    cronbach_alpha = round(alpha, 2),
    max_abs_r = round(max(abs(r), na.rm = TRUE), 2),
    pairs_above_threshold = high
  )
}

#' Response-category shares per item, used to justify collapsing sparse
#' categories (rule of thumb: categories below ~10% are hard for LCA to use)
category_shares <- function(df, items, sparse_threshold = 10) {
  purrr::map_dfr(items, function(v) {
    tab <- table(df[[v]], useNA = "no")
    tibble::tibble(item = v,
                   category = names(tab),
                   n = as.vector(tab),
                   pct = round(100 * as.vector(tab) / sum(tab), 1))
  }) |>
    dplyr::group_by(item) |>
    dplyr::summarise(
      n_categories = dplyr::n(),
      min_pct = min(pct),
      n_sparse = sum(pct < sparse_threshold),
      .groups = "drop"
    )
}

#' Weighted shares of a categorical variable, excluding non-response codes
weighted_shares <- function(df, var, weight = "weight", exclude = NULL) {
  df |>
    dplyr::filter(!is.na(.data[[var]]), !(.data[[var]] %in% exclude)) |>
    dplyr::group_by(category = haven::as_factor(.data[[var]])) |>
    dplyr::summarise(weighted_n = sum(.data[[weight]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(share = round(100 * weighted_n / sum(weighted_n), 1))
}

# ---- Latent class analysis --------------------------------------------------

lca_formula <- function(indicators) {
  as.formula(paste0("cbind(", paste(indicators, collapse = ", "), ") ~ 1"))
}

#' Relative entropy (Ramaswamy et al., 1993): 1 = perfect class separation,
#' 0 = none. pmax() avoids log(0) for posteriors that are exactly zero.
relative_entropy <- function(posterior) {
  p <- pmax(posterior, 1e-12)
  1 - sum(-p * log(p)) / (nrow(p) * log(ncol(p)))
}

#' Fit poLCA models for a range of K and collect fit statistics.
#'
#' na.rm = FALSE keeps respondents with item non-response in the model
#' (full-information ML under MAR) instead of dropping them listwise.
#'
#' llik_scale / n_eff are only needed for the row-replication weighting
#' workaround (see the robustness section): the replicated data inflate the
#' log-likelihood and N, so information criteria are recomputed on the scale
#' of the original sample.
fit_lca_range <- function(data, indicators, ks = 2:6, nrep = 10,
                          maxiter = 5000, llik_scale = 1,
                          n_eff = nrow(data)) {
  f <- lca_formula(indicators)
  purrr::map(ks, function(k) {
    m <- poLCA::poLCA(f, data = data, nclass = k, nrep = nrep,
                      maxiter = maxiter, na.rm = FALSE, verbose = FALSE)
    ll <- m$llik * llik_scale
    list(
      model = m,
      stats = tibble::tibble(
        K = k,
        LL = round(ll, 1),
        npar = m$npar,
        BIC = round(-2 * ll + m$npar * log(n_eff), 1),
        AIC = round(-2 * ll + 2 * m$npar, 1),
        entropy = round(relative_entropy(m$posterior), 3),
        min_class_pct = round(100 * min(m$P), 1)
      )
    )
  }) |>
    purrr::set_names(paste0("K", ks))
}

fit_table <- function(fits) purrr::map_dfr(fits, "stats")

#' Conditional item-response probabilities in long format
tidy_cond_probs <- function(model,
                            levels = c("Low", "Medium", "High")) {
  purrr::imap_dfr(model$probs, function(mat, item) {
    tibble::as_tibble(mat, rownames = "class") |>
      tidyr::pivot_longer(-class, names_to = "category",
                          values_to = "probability") |>
      dplyr::mutate(
        item = item,
        class = as.integer(gsub("\\D", "", class)),
        category = factor(levels[as.integer(gsub("\\D", "", category))],
                          levels = levels)
      )
  }) |>
    dplyr::select(item, class, category, probability)
}

#' Heuristic check of local independence: polychoric correlations between
#' indicators within each modally assigned class. Large residual correlations
#' suggest the class structure does not fully explain the item association.
check_local_independence <- function(data, indicators, class_var,
                                     threshold = 0.7) {
  classes <- sort(unique(stats::na.omit(data[[class_var]])))
  purrr::map_dfr(classes, function(k) {
    x <- data[!is.na(data[[class_var]]) & data[[class_var]] == k, indicators]
    # polychoric() needs variation in every item
    x <- x[, vapply(x, function(v) length(unique(stats::na.omit(v))) > 1,
                    logical(1)), drop = FALSE]
    rho <- tryCatch(
      suppressWarnings(psych::polychoric(as.data.frame(x))$rho),
      error = function(e) NULL)
    if (is.null(rho)) {
      return(tibble::tibble(class = k, item_1 = NA, item_2 = NA, rho = NA))
    }
    idx <- which(abs(rho) > threshold & upper.tri(rho), arr.ind = TRUE)
    if (nrow(idx) == 0) {
      return(tibble::tibble(class = k, item_1 = "none above threshold",
                            item_2 = NA, rho = NA))
    }
    tibble::tibble(class = k,
                   item_1 = rownames(rho)[idx[, 1]],
                   item_2 = colnames(rho)[idx[, 2]],
                   rho = round(rho[idx], 2))
  })
}

#' Match the classes of two solutions on the same respondents (Hungarian
#' algorithm, maximising overlap of modal assignments) and report agreement.
align_classes <- function(class_a, class_b) {
  tab <- table(class_a, class_b)
  m <- matrix(0, max(dim(tab)), max(dim(tab)))
  m[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  assignment <- clue::solve_LSAP(m, maximum = TRUE)
  matched <- sum(m[cbind(seq_len(nrow(m)), as.integer(assignment))])
  list(
    crosstab = tab,
    mapping = tibble::tibble(
      class_a = rownames(tab),
      class_b = colnames(tab)[as.integer(assignment)[seq_len(nrow(tab))]]
    ),
    agreement_pct = round(100 * matched / sum(tab), 1)
  )
}

# ---- Figures ----------------------------------------------------------------

theme_report <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())
}

plot_cond_probs <- function(cond_probs, item_labels = NULL) {
  d <- cond_probs
  if (!is.null(item_labels)) {
    d$item <- dplyr::recode(d$item, !!!item_labels)
  }
  # keep the indicator order of the model (top to bottom)
  d$item <- factor(d$item, levels = rev(unique(d$item)))
  ggplot2::ggplot(d, ggplot2::aes(x = category, y = item, fill = probability)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", probability)),
                       size = 2.6) +
    ggplot2::facet_wrap(~ paste("Class", class), nrow = 1) +
    ggplot2::scale_fill_gradient(low = "white", high = "#1f4e79",
                                 limits = c(0, 1)) +
    ggplot2::labs(title = "Conditional item-response probabilities by class",
                  x = "Endorsement", y = NULL, fill = "Probability") +
    theme_report(11)
}

plot_group_shares <- function(df, group_labels, group_colours,
                              class_labels = NULL, weight = "weight") {
  d <- df |>
    dplyr::filter(!is.na(class),
                  vote %in% as.numeric(names(group_labels))) |>
    dplyr::mutate(group = factor(group_labels[as.character(vote)],
                                 levels = group_labels)) |>
    dplyr::group_by(class, group) |>
    dplyr::summarise(w = sum(.data[[weight]], na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(class) |>
    dplyr::mutate(pct = 100 * w / sum(w)) |>
    dplyr::ungroup()
  if (!is.null(class_labels)) {
    d$class <- factor(class_labels[as.character(d$class)],
                      levels = rev(class_labels))
  } else {
    d$class <- factor(paste("Class", d$class),
                      levels = rev(paste("Class", sort(unique(d$class)))))
  }
  ggplot2::ggplot(d, ggplot2::aes(x = pct, y = class, fill = group)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(pct >= 4,
                                                  sprintf("%.0f%%", pct), "")),
                       position = ggplot2::position_stack(vjust = 0.5),
                       colour = "white", size = 3) +
    ggplot2::scale_fill_manual(values = group_colours) +
    ggplot2::labs(title = "Covariate composition by latent class",
                  subtitle = "Weighted shares within each class",
                  x = NULL, y = NULL, fill = "Group") +
    theme_report() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}

plot_class_profiles <- function(df, weight = "weight") {
  d <- df |>
    dplyr::filter(!is.na(class)) |>
    dplyr::group_by(class) |>
    dplyr::summarise(
      `Tertiary education` = weighted.mean(education_cat == "Tertiary education",
                                           .data[[weight]], na.rm = TRUE),
      Female = weighted.mean(sex_f == "Female", .data[[weight]], na.rm = TRUE),
      `Binary covariate` = weighted.mean(binary_cov == 1, .data[[weight]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(-class, names_to = "variable", values_to = "share")
  # grouped bars rather than lines: the three variables are separate
  # categories, so connecting them would suggest a continuum that isn't there
  ggplot2::ggplot(d, ggplot2::aes(variable, share, fill = factor(class))) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.75, colour = "grey20") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * share)),
                       position = ggplot2::position_dodge(width = 0.8),
                       vjust = -0.4, size = 3) +
    ggplot2::scale_fill_grey(start = 0.25, end = 0.85) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(title = "Socio-demographic profiles by latent class",
                  x = NULL, y = "Weighted share", fill = "Class") +
    theme_report()
}

plot_age_by_class <- function(df, weight = "weight") {
  d <- dplyr::filter(df, !is.na(class)) |>
    dplyr::mutate(class = factor(class), w = .data[[weight]])
  means <- d |>
    dplyr::group_by(class) |>
    dplyr::summarise(age = weighted.mean(age, w, na.rm = TRUE))
  # densities are weighted by replicating respondents in proportion to their
  # weight (works across ggplot2 versions, unlike the weight aesthetic)
  d_rep <- tidyr::uncount(d, pmax(round(w * 10), 1))
  ggplot2::ggplot(d_rep, ggplot2::aes(class, age)) +
    ggplot2::geom_violin(fill = "grey80", colour = "grey30") +
    ggplot2::geom_point(data = means, shape = 21, fill = "white", size = 3) +
    ggplot2::labs(title = "Weighted age distribution by latent class",
                  subtitle = "Points show weighted means",
                  x = "Latent class", y = "Age (years)") +
    theme_report()
}

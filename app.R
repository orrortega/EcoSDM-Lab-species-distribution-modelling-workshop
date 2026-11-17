# EcoSDM Lab | species distribution modelling workshop
# Simulated data by default; optionally analyse uploaded presence/absence CSV files.
# Run with shiny::runApp("path/to/eco_sdm_workshop")

required_packages <- c("shiny", "bslib", "ggplot2", "mgcv", "ranger")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}

library(shiny)
library(ggplot2)
library(mgcv)
options(shiny.maxRequestSize = 20 * 1024^2)  # 20 MB CSV uploads

sim_landscape <- read.csv("data/simulated_landscape.csv")
sim_observations <- read.csv("data/simulated_observations.csv")
sim_labels <- c(
  grassland = "Semi-natural grassland (%)",
  temperature = "Mean temperature (°C)",
  rainfall = "Annual rainfall (mm)",
  elevation = "Elevation (m)",
  distance_road = "Distance to road (km)"
)
model_names <- c(glm = "GLM", gam = "GAM", rf = "Random Forest")
model_colors <- c(glm = "#296b92", gam = "#c17c26", rf = "#277653")
map_colours <- c("#f3eee2", "#c9d9a9", "#6baf82", "#286b56", "#103e36")

candidate_targets <- data.frame(
  site = c("A · Western valley", "B · Central meadow", "C · Eastern uplands"),
  east = c(20, 50, 80), north = c(25, 45, 55)
)
# Candidate sites stay at the same locations, but use the active simulation truth.
get_candidate_sites <- function(landscape) {
  do.call(rbind, lapply(seq_len(nrow(candidate_targets)), function(i) {
    target <- candidate_targets[i, ]
    k <- which.min((landscape$x_km - target$east)^2 +
                     (landscape$y_km - target$north)^2)
    data.frame(site = target$site, short_label = LETTERS[i],
               landscape[k, , drop = FALSE], row.names = NULL)
  }))
}

# Students change ecological niches, sampling locations, and sampling intensity.
# Calibrating the intercept to a chosen mean probability makes extreme niche
# settings feasible without hiding the effects of the other ecological controls.
simulate_own_species <- function(landscape, n_surveys, grassland_effect,
                                 temperature_opt, temperature_width,
                                 rainfall_opt, rainfall_width, prevalence,
                                 sampling, seed) {
  d <- landscape
  grass <- (d$grassland - 50) / 25
  temp <- (d$temperature - temperature_opt) / temperature_width
  rain <- (d$rainfall - rainfall_opt) / rainfall_width
  ecological_score <- grassland_effect * grass - 0.5 * temp^2 -
    0.5 * rain^2 + 0.15 * grass * rain

  target <- prevalence / 100
  objective <- function(a) mean(stats::plogis(ecological_score + a)) - target
  # This bracket guarantees a root even for very narrow ecological niches.
  bracket <- c(-max(ecological_score) - 30, -min(ecological_score) + 30)
  intercept <- stats::uniroot(objective, interval = bracket)$root
  d$true_probability <- stats::plogis(ecological_score + intercept)

  # Road-biased sampling makes visits more likely near roads, not more likely
  # wherever the species is present. Presence remains a Bernoulli outcome.
  weights <- if (identical(sampling, "road"))
    exp(-d$distance_road / 5) + 0.03 else NULL
  set.seed(as.integer(seed))
  ids <- sample.int(nrow(d), size = as.integer(n_surveys), replace = FALSE,
                    prob = weights)
  observations <- d[ids, setdiff(names(d), "true_probability"), drop = FALSE]
  observations$site_id <- sprintf("C%04d", seq_along(ids))
  observations$presence <- stats::rbinom(length(ids), size = 1,
                                          prob = d$true_probability[ids])
  rownames(observations) <- NULL
  list(landscape = d, observations = observations,
       expected_prevalence = mean(d$true_probability),
       observed_prevalence = mean(observations$presence))
}

# -------------------- Data preparation / input checking ---------------------
read_upload <- function(file) {
  tryCatch(
    read.csv(file$datapath, check.names = TRUE, stringsAsFactors = FALSE),
    error = function(e) validate(need(FALSE, paste("Cannot read CSV:", conditionMessage(e))))
  )
}

prepare_data <- function(surveys, grid = NULL, response_col, x_col, y_col,
                         simulated = FALSE) {
  needed <- c(response_col, x_col, y_col)
  validate(need(all(nzchar(needed)) && length(unique(needed)) == 3,
                "Choose three different columns for presence (0/1), X and Y."))
  validate(need(all(needed %in% names(surveys)),
                "One of the chosen response or coordinate columns is missing from the survey CSV."))
  validate(need(nrow(surveys) >= 40,
                "Please provide at least 40 survey rows for this teaching app."))

  # Use only numeric predictor columns, with reasonable variation.
  candidates <- setdiff(names(surveys), c(needed, "presence", "x_km", "y_km", "true_probability"))
  candidates <- candidates[vapply(surveys[candidates], is.numeric, logical(1))]
  candidates <- candidates[vapply(surveys[candidates], function(v) {
    length(unique(v[is.finite(v)])) >= 5L
  }, logical(1))]
  if (!is.null(grid)) {
    # Extra numeric survey fields need not appear in a separate prediction grid.
    candidates <- candidates[candidates %in% names(grid)]
    candidates <- candidates[vapply(grid[candidates], is.numeric, logical(1))]
  }
  validate(need(length(candidates) >= 1,
                "No compatible numeric predictors were found. Each predictor must have at least 5 distinct survey values and, if a prediction CSV is provided, appear as a numeric column there too."))
  validate(need(length(candidates) <= 35,
                "This teaching app supports up to 35 numeric predictors. Remove unrelated numeric columns first."))

  response_text <- trimws(as.character(surveys[[response_col]]))
  empty_response <- is.na(response_text) | response_text == ""
  validate(need(all(empty_response | response_text %in% c("0", "1")),
                "Presence/absence must use 0 and 1 (not yes/no, counts or probabilities)."))
  response <- suppressWarnings(as.numeric(response_text))
  validate(need(is.numeric(surveys[[x_col]]) && is.numeric(surveys[[y_col]]),
                "X and Y survey coordinates must be numeric columns."))

  source_rows <- nrow(surveys)
  surveys$presence <- response
  surveys$x_km <- surveys[[x_col]]
  surveys$y_km <- surveys[[y_col]]
  obs_keep <- complete.cases(surveys[, c("presence", "x_km", "y_km", candidates), drop = FALSE]) &
    apply(surveys[, c("x_km", "y_km", candidates), drop = FALSE], 1L,
          function(v) all(is.finite(v)))
  surveys <- surveys[obs_keep, , drop = FALSE]
  validate(need(nrow(surveys) >= 40,
                "Fewer than 40 complete surveys remain after removing missing or non-finite values."))
  counts <- table(factor(surveys$presence, levels = c(0, 1)))
  validate(need(all(counts >= 8),
                "At least 8 absences and 8 presences are needed for the workshop's model comparisons."))
  validate(need(length(unique(surveys$x_km)) >= 4L && length(unique(surveys$y_km)) >= 4L,
                "Survey coordinates need variation in both X and Y to support spatial diagnostics."))
  validate(need(nrow(surveys) <= 10000L,
                "For responsiveness, limit the survey CSV to 10,000 rows."))

  has_grid <- !is.null(grid)
  grid_rows <- if (has_grid) nrow(grid) else NA_integer_
  if (has_grid) {
    grid_required <- c(x_col, y_col, candidates)
    validate(need(all(grid_required %in% names(grid)),
                  paste("Prediction CSV is missing columns:",
                        paste(setdiff(grid_required, names(grid)), collapse = ", "))))
    validate(need(all(vapply(grid[grid_required], is.numeric, logical(1))),
                  "Prediction coordinates and predictor columns must all contain numeric values."))
    grid$x_km <- grid[[x_col]]
    grid$y_km <- grid[[y_col]]
    grid_keep <- complete.cases(grid[, c("x_km", "y_km", candidates), drop = FALSE]) &
      apply(grid[, c("x_km", "y_km", candidates), drop = FALSE], 1L,
            function(v) all(is.finite(v)))
    grid <- grid[grid_keep, , drop = FALSE]
    validate(need(nrow(grid) >= 1L,
                  "No complete prediction locations remain in the prediction CSV."))
    validate(need(nrow(grid) <= 25000L,
                  "For responsiveness, limit the prediction CSV to 25,000 locations."))
  } else {
    # No spatial predictions beyond the surveyed locations without covariate data.
    grid <- surveys
  }

  labels <- if (simulated) sim_labels[candidates] else stats::setNames(candidates, candidates)
  list(observations = surveys, landscape = grid, predictors = candidates,
       labels = labels, simulated = simulated, has_grid = has_grid,
       points_only = !simulated, x_col = x_col, y_col = y_col,
       n_original = source_rows, n_survey_dropped = source_rows - nrow(surveys),
       n_grid_original = grid_rows,
       n_grid_dropped = if (has_grid) grid_rows - nrow(grid) else 0L)
}

# --------------------------- Model functions --------------------------------
make_split <- function(strategy, observations) {
  if (identical(strategy, "spatial")) {
    cutoff <- as.numeric(stats::quantile(observations$x_km, 0.75))
    test_ids <- which(observations$x_km >= cutoff)
    validate(need(length(test_ids) >= 2L && length(test_ids) < nrow(observations) - 2L,
                  "Spatial split could not create a useful test set; try random holdout."))
  } else {
    set.seed(42)
    # Stratify by observed outcome so a random test set includes both classes.
    test_ids <- unlist(lapply(c(0, 1), function(y) {
      ids <- which(observations$presence == y)
      sample(ids, size = max(1L, min(length(ids) - 1L, round(length(ids) * 0.25))))
    }), use.names = FALSE)
  }
  train <- observations[-test_ids, , drop = FALSE]
  test <- observations[test_ids, , drop = FALSE]
  validate(need(length(unique(train$presence)) == 2L,
                "Training sites contain only one presence class with this split; try random holdout."))
  list(train = train, test = test)
}

fit_sdm <- function(type, vars, train, gam_k = 6, trees = 250, leaf = 8,
                    gam_smooth_vars = NULL) {
  validate(need(length(vars) >= 1L, "Select at least one predictor in this model's sidebar."))
  validate(need(all(vars %in% names(train)), "Selected predictor is unavailable in the dataset."))
  # With the standalone helper, NULL means smooth all predictors (backwards
  # compatible). Shiny passes a character vector, including character(0), so
  # students can explicitly choose an entirely linear GAM for comparison.
  if (identical(type, "gam")) {
    if (is.null(gam_smooth_vars)) gam_smooth_vars <- vars
    gam_smooth_vars <- intersect(vars, gam_smooth_vars)
  }
  fit <- tryCatch({
    if (identical(type, "glm")) {
      f <- stats::reformulate(vars, response = "presence")
      stats::glm(f, data = train, family = stats::binomial())
    } else if (identical(type, "gam")) {
      # A GAM can mix parametric linear terms and penalised regression splines.
      ks <- vapply(gam_smooth_vars, function(v) {
        min(as.integer(gam_k), length(unique(train[[v]])) - 1L)
      }, integer(1))
      if (any(ks < 4L)) {
        stop("A selected smooth has too few unique training values. Switch it to linear or choose another predictor.")
      }
      terms <- vapply(vars, function(v) {
        if (v %in% gam_smooth_vars) {
          sprintf("s(%s, k = %d, bs = 'cr')", v, ks[[v]])
        } else v
      }, character(1))
      f <- stats::as.formula(paste("presence ~", paste(terms, collapse = " + ")))
      mgcv::gam(f, data = train, family = stats::binomial(), method = "REML")
    } else if (identical(type, "rf")) {
      dat <- data.frame(
        presence_factor = factor(ifelse(train$presence == 1, "present", "absent"),
                                 levels = c("absent", "present")),
        train[, vars, drop = FALSE], check.names = FALSE
      )
      ranger::ranger(
        presence_factor ~ ., data = dat, probability = TRUE,
        num.trees = as.integer(trees), min.node.size = as.integer(leaf),
        importance = "permutation", seed = 42, num.threads = 1
      )
    } else stop("Unknown model")
  }, error = function(e) validate(need(FALSE, paste(
    "Could not fit this model:", conditionMessage(e),
    "Try fewer predictors or a larger dataset."))))
  list(type = type, vars = vars, fit = fit,
       smooth_vars = if (identical(type, "gam")) gam_smooth_vars else character(0),
       gam_k = if (identical(type, "gam")) gam_k else NA_integer_)
}

predict_sdm <- function(model, newdata) {
  p <- if (identical(model$type, "rf")) {
    predict(model$fit, data = newdata[, model$vars, drop = FALSE])$predictions[, "present"]
  } else {
    stats::predict(model$fit, newdata = newdata, type = "response")
  }
  pmin(1, pmax(0, as.numeric(p)))
}

auc_rank <- function(y, p) {
  pos <- sum(y == 1)
  neg <- sum(y == 0)
  if (pos == 0 || neg == 0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - pos * (pos + 1) / 2) / (pos * neg)
}

scores <- function(y, p, threshold = 0.5) {
  decision <- as.integer(p >= threshold)
  actual_pos <- sum(y == 1)
  actual_neg <- sum(y == 0)
  safe_p <- pmin(1 - 1e-12, pmax(1e-12, p))
  c(AUC = auc_rank(y, p), Brier = mean((y - p)^2),
    LogLoss = -mean(y * log(safe_p) + (1 - y) * log(1 - safe_p)),
    Sensitivity = if (actual_pos) sum(y == 1 & decision == 1) / actual_pos else NA_real_,
    Specificity = if (actual_neg) sum(y == 0 & decision == 0) / actual_neg else NA_real_)
}

roc_points <- function(y, p, name) {
  if (!sum(y == 1) || !sum(y == 0)) {
    return(data.frame(model = character(), fpr = numeric(), tpr = numeric()))
  }
  thresholds <- seq(1, 0, length.out = 151)
  data.frame(model = name,
             fpr = vapply(thresholds, function(t) mean(p[y == 0] >= t), numeric(1)),
             tpr = vapply(thresholds, function(t) mean(p[y == 1] >= t), numeric(1)))
}

calibration_points <- function(y, p, name) {
  bin <- cut(p, breaks = seq(0, 1, length.out = 9),
             include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(sort(unique(bin[!is.na(bin)])), function(i) {
    ids <- which(bin == i)
    data.frame(model = name, predicted = mean(p[ids]),
               observed = mean(y[ids]), n = length(ids))
  }))
}

response_data <- function(model, variable, train, grid, n_grid = 37) {
  values <- seq(min(grid[[variable]]), max(grid[[variable]]), length.out = n_grid)
  rows <- unique(round(seq(1, nrow(train), length.out = min(250, nrow(train)))))
  reference <- train[rows, , drop = FALSE]
  estimate <- vapply(values, function(value) {
    modified <- reference
    modified[[variable]] <- value
    mean(predict_sdm(model, modified))
  }, numeric(1))
  data.frame(x = values, estimate = estimate)
}

plot_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", colour = "#143d35", size = 14),
          plot.subtitle = element_text(colour = "#64766d", size = 10),
          axis.title = element_text(colour = "#385248"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

spatial_layer <- function(data, variable, simulated = FALSE, limits = NULL,
                          legend_title = "Probability") {
  # For arbitrary uploaded points, do NOT pretend there is a regular raster.
  p <- ggplot(data, aes(x_km, y_km))
  if (isTRUE(simulated)) {
    p <- p + geom_raster(aes(fill = .data[[variable]]), interpolate = FALSE) +
      scale_fill_gradientn(colours = map_colours, limits = limits, name = legend_title)
  } else {
    p <- p + geom_point(aes(colour = .data[[variable]]), shape = 15,
                        size = if (nrow(data) > 10000) 0.65 else 1.8, alpha = 0.9) +
      scale_colour_gradientn(colours = map_colours, limits = limits, name = legend_title)
  }
  p + coord_equal(expand = FALSE) + plot_theme() +
    theme(panel.grid = element_blank())
}

# Signed probability difference between GAM and GLM across prediction locations.
# This is a difference in model predictions, not a residual or uncertainty interval.
glm_gam_difference <- function(preds) {
  stopifnot(all(c("glm", "gam") %in% names(preds)),
            length(preds$glm) == length(preds$gam))
  as.numeric(preds$gam) - as.numeric(preds$glm)
}

plot_glm_gam_difference <- function(ds, preds) {
  d <- ds$landscape
  d$probability_difference <- glm_gam_difference(preds)
  stopifnot(nrow(d) == length(d$probability_difference))
  max_difference <- max(abs(d$probability_difference), na.rm = TRUE)
  # A symmetric, data-dependent colour scale preserves the sign and makes small
  # differences visible. The legend displays its actual numeric range.
  bound <- if (is.finite(max_difference) && max_difference > 1e-8)
    max_difference else 0.01
  colour_args <- list(low = "#b65342", mid = "#f9f8f3", high = "#247e78",
                      midpoint = 0, limits = c(-bound, bound),
                      breaks = c(-bound, 0, bound),
                      labels = sprintf("%.3g", c(-bound, 0, bound)),
                      name = "GAM − GLM\nprobability")
  p <- ggplot(d, aes(x_km, y_km))
  if (isTRUE(ds$simulated)) {
    p <- p + geom_raster(aes(fill = probability_difference), interpolate = FALSE) +
      do.call(scale_fill_gradient2, colour_args)
  } else {
    p <- p + geom_point(aes(colour = probability_difference), shape = 15,
                        size = if (nrow(d) > 10000) 0.65 else 1.8,
                        alpha = 0.9) +
      do.call(scale_colour_gradient2, colour_args)
  }
  p + coord_equal(expand = FALSE) +
    labs(title = "GAM − GLM · predicted probability",
         subtitle = "Teal: GAM predicts higher | Rust: GLM predicts higher",
         x = if (ds$simulated) "Fictional km east" else "X coordinate",
         y = if (ds$simulated) "Fictional km north" else "Y coordinate") +
    plot_theme() + theme(panel.grid = element_blank())
}

plot_prediction <- function(model, show_observations, train, ds) {
  grid <- ds$landscape
  grid$estimate <- predict_sdm(model, grid)
  p <- spatial_layer(grid, "estimate", ds$simulated, limits = c(0, 1),
                     legend_title = "Predicted\noccurrence") +
    labs(title = paste0(model_names[model$type], " · predicted distribution"),
         subtitle = if (ds$simulated) "Fictional landscape grid"
                    else if (ds$has_grid) "Predictions at uploaded environmental locations"
                    else "Predictions ONLY at surveyed locations — no prediction grid uploaded",
         x = if (ds$simulated) "Easting (fictional km)" else "X coordinate",
         y = if (ds$simulated) "Northing (fictional km)" else "Y coordinate")
  if (isTRUE(show_observations)) {
    found <- train[train$presence == 1, , drop = FALSE]
    p <- p + geom_point(data = found, aes(x_km, y_km), inherit.aes = FALSE,
                        shape = 21, size = 1.7, fill = "white",
                        colour = "#203c37", stroke = 0.5, alpha = 0.85)
  }
  p
}

plot_response <- function(model, variable, train, ds, baseline = NULL) {
  d <- response_data(model, variable, train, ds$landscape)
  omitted <- !variable %in% model$vars
  p <- ggplot(d, aes(x, estimate)) +
    geom_line(linewidth = 1.25, colour = model_colors[[model$type]])
  if (!is.null(baseline)) {
    b <- response_data(baseline, variable, train, ds$landscape)
    p <- p + geom_line(data = b, aes(x, estimate), inherit.aes = FALSE,
                       colour = model_colors[["glm"]], linetype = "dashed", linewidth = 1.1)
  }
  p + scale_y_continuous(limits = c(0, 1)) +
    labs(title = if (is.null(baseline)) "Average model response" else "GAM (solid) vs GLM (dashed)",
         subtitle = if (omitted) "Not in this model: no direct effect"
                    else "Vary one predictor, averaging across surveyed environments",
         x = ds$labels[[variable]], y = "Mean predicted occurrence") + plot_theme()
}

# On the log-odds (link) scale, the shape is unambiguous. Unlike the average
# probability response above, this shows the selected model term itself.
plot_gam_partial <- function(model, variable, train, ds, n_grid = 110) {
  validate(need(model$type == "gam", "Partial smooths are only available for GAMs."))
  validate(need(variable %in% model$vars,
                "Select a predictor that is included in your GAM to see its term."))
  values <- seq(min(train[[variable]]), max(train[[variable]]), length.out = n_grid)
  newdata <- train[rep(1L, n_grid), , drop = FALSE]
  for (v in model$vars) newdata[[v]] <- stats::median(train[[v]])
  newdata[[variable]] <- values
  is_smooth <- variable %in% model$smooth_vars
  term <- if (is_smooth) paste0("s(", variable, ")") else variable
  pred <- stats::predict(model$fit, newdata = newdata, type = "terms",
                         terms = term, se.fit = TRUE)
  validate(need(term %in% colnames(pred$fit),
                "The chosen model term could not be extracted."))
  effect <- as.numeric(pred$fit[, term])
  se <- as.numeric(pred$se.fit[, term])
  d <- data.frame(x = values, effect = effect,
                  lower = effect - 1.96 * se, upper = effect + 1.96 * se)
  st <- as.matrix(summary(model$fit)$s.table)
  edf <- if (is_smooth && nrow(st) && term %in% rownames(st)) st[term, "edf"] else 1
  ggplot(d, aes(x, effect)) +
    geom_hline(yintercept = 0, colour = "#a6b8af", linetype = "dotted") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = model_colors[["gam"]], alpha = 0.15) +
    geom_line(colour = model_colors[["gam"]], linewidth = 1.2) +
    labs(title = if (is_smooth) sprintf("Smooth term · EDF %.2f", edf)
                     else "Linear term · EDF 1",
         subtitle = "Approximate 95% pointwise interval; other predictors held fixed",
         x = ds$labels[[variable]], y = "Contribution to log-odds (centred)") + plot_theme()
}

model_explanation <- function(model, labels) {
  flag <- function(p) ifelse(is.na(p), "Not estimable",
                            ifelse(p < 0.05, "Significant", "Not significant"))
  if (model$type == "glm") {
    ct <- stats::coef(summary(model$fit))
    p <- ct[, "Pr(>|z|)"]
    status <- flag(p)
    status[rownames(ct) == "(Intercept)"] <- "Intercept"
    data.frame(Term = rownames(ct),
               `Estimate (log-odds)` = round(ct[, "Estimate"], 3),
               `Std. error` = round(ct[, "Std. Error"], 3),
               `p-value` = format.pval(p, digits = 3, eps = 0.001),
               `At 5% level` = status,
               check.names = FALSE, row.names = NULL)
  } else if (model$type == "gam") {
    sm <- summary(model$fit)
    pt <- as.matrix(sm$p.table)
    st <- as.matrix(sm$s.table)
    parametric <- if (nrow(pt)) data.frame(
      Term = rownames(pt), Type = "Linear",
      EDF = ifelse(rownames(pt) == "(Intercept)", NA_real_, 1),
      `Estimate (log-odds)` = round(pt[, "Estimate"], 3),
      `p-value` = format.pval(pt[, "Pr(>|z|)"], digits = 3, eps = 0.001),
      `At 5% level` = ifelse(rownames(pt) == "(Intercept)", "Intercept", flag(pt[, "Pr(>|z|)"])),
      check.names = FALSE, row.names = NULL) else NULL
    smooths <- if (nrow(st)) data.frame(
      Term = rownames(st), Type = "Spline",
      EDF = round(st[, "edf"], 2),
      `Estimate (log-odds)` = NA_real_,
      `p-value` = format.pval(st[, "p-value"], digits = 3, eps = 0.001),
      `At 5% level` = flag(st[, "p-value"]),
      check.names = FALSE, row.names = NULL) else NULL
    rbind(parametric, smooths)
  } else {
    vi <- sort(model$fit$variable.importance, decreasing = TRUE)
    data.frame(Variable = unname(labels[names(vi)]),
               `Permutation importance` = round(unname(vi), 4),
               `p-value` = "Not available", check.names = FALSE)
  }
}

significant_predictors <- function(model) {
  if (identical(model$type, "glm")) {
    ct <- stats::coef(summary(model$fit))
    p <- ct[, "Pr(>|z|)"]
    return(model$vars[!is.na(p[model$vars]) & p[model$vars] < 0.05])
  }
  if (identical(model$type, "gam")) {
    sm <- summary(model$fit)
    pt <- as.matrix(sm$p.table)
    st <- as.matrix(sm$s.table)
    result <- vapply(model$vars, function(v) {
      if (v %in% model$smooth_vars) {
        term <- paste0("s(", v, ")")
        p <- if (term %in% rownames(st)) st[term, "p-value"] else NA_real_
      } else {
        p <- if (v %in% rownames(pt)) pt[v, "Pr(>|z|)"] else NA_real_
      }
      !is.na(p) && p < 0.05
    }, logical(1))
    return(model$vars[result])
  }
  character(0)
}

correlation_plot <- function(ds, vars, method = "pearson") {
  validate(need(length(vars) >= 2L,
                "Select at least two predictors to see their correlation."))
  validate(need(length(vars) <= 12L,
                "For a readable matrix, select at most 12 predictors."))
  validate(need(all(vars %in% ds$predictors), "Choose available predictors."))
  mat <- stats::cor(ds$observations[, vars, drop = FALSE], method = method)
  display_labels <- unname(ds$labels[vars])
  rownames(mat) <- display_labels
  colnames(mat) <- display_labels
  d <- as.data.frame(as.table(mat), responseName = "r")
  names(d)[1:2] <- c("Variable1", "Variable2")
  ggplot(d, aes(Variable1, Variable2, fill = r)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = sprintf("%.2f", r),
                  colour = abs(r) >= 0.62), size = if (length(vars) > 8) 2.7 else 3.3) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#263b34"), guide = "none") +
    scale_fill_gradient2(low = "#bf5848", mid = "#f7f8f3", high = "#257c68",
                         midpoint = 0, limits = c(-1, 1), name = "Correlation") +
    coord_equal() +
    labs(title = paste(if (method == "pearson") "Pearson" else "Spearman",
                       "correlation across surveyed sites"), x = NULL, y = NULL) +
    plot_theme() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9),
          panel.grid = element_blank())
}

model_formulation <- function(model) {
  if (model$type == "glm") {
    paste0("logit(P(presence)) = intercept + ", paste(model$vars, collapse = " + "))
  } else if (model$type == "gam") {
    # Print the ACTUAL fitted formula, including spline terms and any linear terms.
    paste(deparse(stats::formula(model$fit), width.cutoff = 500L), collapse = "")
  } else paste0("Probability forest using: ", paste(model$vars, collapse = ", "))
}

# General mathematical expressions, rather than equations fitted to any
# particular combination of selected predictors.
equation_card <- function(type) {
  eq <- switch(type,
    data = "$$Y_i\\mid p_i \\sim \\mathrm{Bernoulli}(p_i),\\qquad p_i=P(Y_i=1\\mid \\mathbf{x}_i)$$",
    glm = "$$Y_i\\sim\\mathrm{Bernoulli}(p_i),\\qquad\\mathrm{logit}(p_i)=\\beta_0+\\sum_{j=1}^{q}\\beta_j x_{ij}$$",
    gam = "$$Y_i\\sim\\mathrm{Bernoulli}(p_i),\\qquad\\mathrm{logit}(p_i)=\\beta_0+\\sum_{j\\in S} f_j(x_{ij})+\\sum_{j\\in L}\\beta_j x_{ij}$$",
    rf = "$$\\hat p(\\mathbf{x})=\\frac{1}{B}\\sum_{b=1}^{B}\\hat p_b(\\mathbf{x})$$",
    compare = "$$\\mathrm{Brier}=\\frac{1}{n_{test}}\\sum_{i=1}^{n_{test}}(y_i-\\hat p_i)^2$$")
  interpretation <- switch(type,
    data = "A survey is a Bernoulli observation: 1 = presence, 0 = absence. This teaching example assumes perfect detection.",
    glm = "Each selected predictor has a constant slope on the log-odds scale. Coefficients are conditional on the other selected predictors.",
    gam = "Tick a predictor to use a cubic regression spline; unticked predictors remain linear. Smoothness is penalised by REML. Smooth p-values are approximate and do not test nonlinearity.",
    rf = "Each tree provides an estimated probability; the forest averages predictions. There are no standard per-predictor coefficient p-values.",
    compare = "The Brier score is the mean squared error of predicted occurrence probabilities on held-out surveys. Smaller is better.")
  div(class = "equation-card",
      div(class = "equation-heading", "THE MODEL IN ONE EQUATION"),
      withMathJax(HTML(eq)),
      div(class = "hint", interpretation))
}

# ---------------------------------- UI --------------------------------------
model_tab <- function(type, title, description, question) {
  tabPanel(
    title = title, value = type,
    div(class = "model-intro", h3(title), p(description)),
    equation_card(type),
    fluidRow(
      column(3,
             div(class = "panel-card controls", h4("Design your model"),
                 checkboxGroupInput(paste0("vars_", type), "Environmental predictors",
                                    choices = stats::setNames(names(sim_labels), sim_labels),
                                    selected = c("grassland", "temperature", "rainfall")),
                 if (type %in% c("glm", "gam")) tagList(
                   actionButton(paste0("keep_significant_", type),
                                "Keep p < 0.05 predictors", class = "btn-outline-primary"),
                   div(class = "hint", "Exploratory shortcut: uses training-data p-values. ",
                       "Repeated selection makes those same p-values unsuitable for ",
                       "confirmatory inference."),
                   tags$hr()),
                 actionButton(paste0("restore_all_", type), "Select all predictors",
                              class = "btn-outline-secondary"),
                 tags$hr(),
                 if (type == "gam") tagList(
                   checkboxGroupInput("gam_smooth_vars", "Fit cubic regression splines to",
                                      choices = stats::setNames(names(sim_labels), sim_labels),
                                      selected = c("temperature", "rainfall")),
                   div(class = "hint", "Tick variables to fit s(x) smooths; unticked selected ",
                       "predictors stay linear. At least one tick shows the spline mechanism."),
                   sliderInput("gam_k", "Maximum spline basis dimension (k)",
                               min = 4, max = 12, value = 8, step = 2),
                   div(class = "hint", "REML may shrink a smooth back to EDF ≈ 1. ",
                       "More basis functions allow flexibility but do not guarantee curves.")),
                 if (type == "rf") tagList(
                   sliderInput("rf_trees", "Number of trees", min = 100, max = 500,
                               value = 250, step = 50),
                   sliderInput("rf_leaf", "Minimum observations per leaf", min = 3,
                               max = 25, value = 8, step = 1)),
                  tags$hr(),
                  actionButton(paste0("fit_", type),
                               paste("Fit", model_names[[type]]),
                               class = "btn-primary"),
                  uiOutput(paste0("fit_status_", type)),
                 tags$hr(),
                 checkboxInput(paste0("show_obs_", type), "Overlay training presences"),
                 selectInput(paste0("effect_", type), "Explore the response to",
                             choices = stats::setNames(names(sim_labels), sim_labels),
                             selected = if (type == "gam") "temperature" else "grassland"),
                  div(class = "hint", "Changes to predictors and model settings do not refit ",
                      "the model until you press Fit. Plot display controls remain interactive. ",
                      "All probability maps share a 0–1 colour scale.")),
             div(class = "panel-card learning-card", h4("Think about..."), p(question))),
      column(9,
             uiOutput(paste0("stats_", type)),
             div(class = "panel-card", plotOutput(paste0("map_", type), height = "365px")),
             fluidRow(
               column(6, div(class = "panel-card",
                             plotOutput(paste0("response_", type), height = "300px"))),
               column(6, div(class = "panel-card", h4("Inside the model"),
                             textOutput(paste0("formula_", type)),
                             tableOutput(paste0("summary_", type)),
                             div(class = "hint", if (type == "rf")
                               "Random Forest importance ranks predictors; it does not provide conventional coefficient significance tests."
                               else if (type == "gam") "EDF near 1 is effectively linear. A significant smooth is not necessarily nonlinear: its p-value tests an effect, not curvature."
                               else "Significance uses Wald tests on training data. p < 0.05 is not an ecological importance threshold.")))),
             if (type == "gam") div(class = "panel-card",
               h4("Inspect the fitted GAM term (on the log-odds scale)"),
               plotOutput("gam_partial_plot", height = "315px"),
               p(class = "hint", "The solid orange curve is the selected term; the band ",
                 "is its approximate pointwise 95% interval. EDF around 1 means ",
                 "the fitted smooth is nearly linear. For a visible ecological optimum, ",
                 "generate a species with narrow temperature tolerance, then include ",
                 "temperature as a smooth."))
    )
  )
  )
}

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#24735c"),
  tags$head(
    tags$title("EcoSDM Lab · Species distribution models"),
    tags$style(HTML("\n      body { background:#f3f6f2; color:#173d33; }\n      .container-fluid { max-width: 1480px; }\n      .hero { background:linear-gradient(118deg,#103e35,#2c765a); color:white;\n              padding:24px 28px; margin:16px 0 18px; border-radius:17px; }\n      .hero h1 { margin:0 0 5px; font-weight:750; font-size:30px; }\n      .hero p { color:#e0eee6; margin:0; max-width:1000px; font-size:15px; }\n      .hero .eyebrow { text-transform:uppercase; letter-spacing:1.4px; font-size:11px;\n                        font-weight:800; color:#a8e2bb; margin-bottom:8px; }\n      .panel-card { background:white; border:1px solid #e0e9e1; border-radius:13px;\n                    padding:17px 18px; margin:0 0 14px; box-shadow:0 2px 12px #183e3408; }\n      .model-intro { margin:20px 0 16px; }\n      .model-intro h3 { margin:0 0 5px; font-weight:750; }\n      .model-intro p { margin:0; color:#51685f; }\n      .controls .checkbox label { padding-bottom:6px; font-weight:450; }\n      .controls h4, .panel-card h4 { font-size:16px; font-weight:750; margin:0 0 13px; }\n      .hint { color:#60796e; font-size:12px; line-height:1.5; }\n      .learning-card { border-left:4px solid #c59c42; background:#fffcf4; }\n      .learning-card p { margin-bottom:0; }\n      .stat-row { display:flex; gap:12px; flex-wrap:wrap; margin:0 0 14px; }\n      .stat-pill { background:white; border:1px solid #dce7df; border-radius:12px;\n                   padding:11px 17px; min-width:145px; flex:1; }\n      .stat-pill .value { font-size:23px; line-height:1.2; font-weight:800; color:#17604b; }\n      .stat-pill .label { margin-top:2px; font-size:12px; color:#597369; }\n      .nav-tabs { margin-top:15px; border-bottom:1px solid #d0dfd3; }\n      .nav-tabs > li > a { font-weight:700; }\n      .tab-content { padding-bottom:28px; }\n      table { font-size:12px; }\n      .table th { background:#f4f7f2; }\n      .shiny-input-container { width:100% !important; }\n      .footer-note { color:#64766e; font-size:12px; margin:12px 0 30px; }\n      .prompt-line { background:#edf6ed; border-radius:10px; padding:12px 14px; }
      .equation-card { padding:14px 20px; margin:0 0 17px; border:1px solid #cfe0d1;
                       background:#f8fcf8; border-radius:12px; overflow-x:auto; }
      .equation-card .equation-heading { font-size:11px; letter-spacing:1.2px;
                                        color:#527362; font-weight:800; }
      .equation-card .MathJax_Display { margin:10px 0 8px; }
      .btn-outline-primary, .btn-outline-secondary { margin-bottom:7px; }
       .fit-ready { margin:9px 0; color:#27674e; font-weight:650; }
       .fit-stale { margin:9px 0; color:#9a5a14; font-weight:650; }
\n    "))
  ),
  div(class = "hero",
      div(class = "eyebrow", "Ecological restoration · interactive workshop"),
      h1("EcoSDM Lab"),
      p("Where should we restore habitat for a fictional grassland butterfly? ",
        "Explore a teaching dataset, design your own simulated species or upload ",
        "binary surveys. Fit three models and compare their predictions.")),
  div(class = "prompt-line", tags$b("Your mission: "),
      "Explore species–environment relationships and ask why models disagree. ",
      "Uploaded data are processed in your Shiny session (on the host server if deployed); ",
      "the built-in example is entirely simulated."),
  tabsetPanel(
    id = "main_tabs", type = "tabs",
    tabPanel(
      title = "1 · Data & upload", value = "data",
      div(class = "model-intro", h3("Understand the data before modelling"),
          p("Start with the fully simulated dataset, or load your own site-level ",
            "presence/absence records. The chosen dataset is used in ALL modelling tabs.")),
      equation_card("data"),
      fluidRow(
        column(4,
               div(class = "panel-card controls", h4("Choose a dataset"),
                   radioButtons("data_source", "Data source",
                                choices = c("Original butterfly example" = "simulation",
                                            "Design a simulated species" = "custom",
                                            "My own CSV files" = "upload"),
                                selected = "simulation"),
                   conditionalPanel(
                     condition = "input.data_source == 'custom'",
                     tags$hr(), h4("Design your species"),
                     sliderInput("sim_grassland", "Grassland association",
                                 min = -1.5, max = 2.5, value = 1.1, step = 0.1),
                     sliderInput("sim_temperature", "Preferred temperature (°C)",
                                 min = 10, max = 17, value = 13, step = 0.5),
                     sliderInput("sim_temp_width", "Temperature tolerance (°C)",
                                 min = 0.75, max = 4, value = 2, step = 0.25),
                     sliderInput("sim_rainfall", "Preferred rainfall (mm)",
                                 min = 800, max = 1350, value = 1050, step = 25),
                     sliderInput("sim_rain_width", "Rainfall tolerance (mm)",
                                 min = 100, max = 450, value = 250, step = 25),
                     sliderInput("sim_prevalence", "Mean occurrence probability (%)",
                                 min = 15, max = 65, value = 35, step = 5),
                     sliderInput("sim_n", "Number of surveyed locations",
                                 min = 200, max = 1200, value = 750, step = 50),
                     selectInput("sim_sampling", "Where do surveys happen?",
                                 choices = c("Across the landscape (random)" = "uniform",
                                             "More surveys near roads (biased)" = "road")),
                     numericInput("sim_seed", "Random seed", value = 2026,
                                  min = 1, max = 999999, step = 1),
                     actionButton("generate_species", "Generate species & surveys",
                                  class = "btn-primary"),
                     div(class = "hint", "After changing settings, press Generate. ",
                          "The app redraws the surveys; click Fit in each model tab, ",
                          "or Fit all in Compare & diagnose, to update model results. ",
                         "Road bias changes sampling locations, not the true occurrence mechanism."),
                     tags$hr()),
                   conditionalPanel(
                     condition = "input.data_source == 'upload'",
                     fileInput("survey_upload", "1. Survey observations (required)",
                               accept = c(".csv", "text/csv")),
                     fileInput("grid_upload", "2. Prediction locations (optional)",
                               accept = c(".csv", "text/csv")),
                     selectInput("response_col", "Response (presence = 1 / absence = 0)",
                                 choices = "presence", selected = "presence"),
                     selectInput("x_col", "X / easting coordinate", choices = "x_km",
                                 selected = "x_km"),
                     selectInput("y_col", "Y / northing coordinate", choices = "y_km",
                                 selected = "y_km"),
                     div(class = "hint", "CSV headers are made into valid R column names ",
                         "(e.g., 'forest cover' becomes 'forest.cover'). ",
                         "Coordinates must be numeric and named the same in both files. ",
                         "All environmental predictors must be numeric.")),
                   div(class = "hint", "You can switch between the original example, ",
                       "your generated species and uploaded surveys at any time.")),
               div(class = "panel-card", h4("Download CSV examples"),
                   downloadButton("download_survey_template", "Example surveys (.csv)"),
                   tags$div(style = "height:10px;"),
                   downloadButton("download_grid_template", "Example prediction grid (.csv)"),
                   p(class = "hint", "The examples also act as templates for column names ",
                     "and data structure. You can inspect them in Excel or R.")),
               div(class = "panel-card learning-card", h4("Important distinction"),
                   p("Survey rows have a known species outcome (0 or 1). Prediction ",
                     "locations only need environmental predictors and coordinates. ",
                     "Without a prediction CSV, the app maps predictions at survey points, ",
                     "not across unsampled space."))),
        column(8,
               uiOutput("data_status"),
               div(class = "panel-card", h4("How the data were created"),
                   uiOutput("data_generation_description")),
               div(class = "panel-card", h4("Explore your active survey data"),
                   selectInput("data_plot_var", "Environmental variable",
                               choices = stats::setNames(names(sim_labels), sim_labels),
                               selected = "grassland"),
                   plotOutput("data_scatter", height = "295px"),
                   tableOutput("data_columns")),
               div(class = "panel-card", h4("Environmental correlation matrix"),
                   selectInput("cor_method", "Correlation type",
                               choices = c("Pearson (linear)" = "pearson",
                                           "Spearman (rank)" = "spearman")),
                   selectizeInput("cor_vars", "Predictors to correlate (max 12)",
                                  choices = stats::setNames(names(sim_labels), sim_labels),
                                  selected = names(sim_labels), multiple = TRUE),
                   plotOutput("correlation_plot", height = "460px"),
                   p(class = "hint", "Correlations are calculated over surveyed locations, ",
                     "not across the full grid. Strongly correlated covariates can destabilise ",
                     "GLM coefficients and complicate interpretation; a large correlation ",
                     "does not prove one variable should be removed.")),
               div(class = "panel-card", h4("First 8 survey records"),
                   tableOutput("data_preview"),
                   div(class = "hint", "Missing or non-finite rows are removed consistently ",
                       "for the selected dataset before all models are fitted."))
        )
      )
    ),
    model_tab("glm", "2 · GLM",
              "A transparent baseline: each predictor has a linear effect on log-odds.",
              "What if the species prefers intermediate conditions?"),
    model_tab("gam", "3 · GAM",
              "Smooth functions let environmental responses curve without making them opaque.",
              "Which responses curve? Does the more flexible fit generalise?"),
    model_tab("rf", "4 · Random Forest",
              "Trees can capture nonlinear responses and interactions between predictors.",
              "Does strong prediction automatically imply an ecological explanation?"),
    tabPanel(
      title = "5 · Compare & diagnose", value = "compare",
      div(class = "model-intro", h3("Compare models on the same unseen sites"),
          p("AUC assesses discrimination; Brier and calibration assess probability estimates. ",
             "After changing the validation method, click Fit all models.")),
      equation_card("compare"),
      fluidRow(
        column(3,
               div(class = "panel-card controls", h4("Validation experiment"),
                   radioButtons("validation", "Choose withheld observations",
                                choices = c("Stratified random ~25%" = "random",
                                            "Spatial holdout: eastern ~25%" = "spatial"),
                                selected = "random"),
                   div(class = "hint", "Spatial holdout tests transfer to a different area. ",
                       "It is a single geographic holdout, not spatial cross-validation."),
                    tags$hr(),
                    actionButton("fit_all_models", "Fit all three models",
                                 class = "btn-primary"),
                    uiOutput("compare_fit_status"),
                   tags$hr(),
                   sliderInput("threshold", "Classification threshold",
                               min = 0.1, max = 0.9, value = 0.5, step = 0.05),
                   selectInput("res_model", "Residuals for",
                               choices = stats::setNames(names(model_names), unname(model_names)),
                               selected = "gam"),
                   conditionalPanel(
                     condition = "input.data_source == 'simulation' || input.data_source == 'custom'",
                     checkboxInput("reveal_truth", "Reveal simulated true suitability"))),
               div(class = "panel-card learning-card", h4("What matters?"),
                   p("Best AUC, calibration, interpretability or ecological plausibility?")),
               div(class = "panel-card", h4("Download active data"),
                   downloadButton("download_observations", "Survey records (.csv)"),
                   tags$div(style = "height:9px;"),
                   downloadButton("download_predictions", "Model predictions (.csv)"),
                   p(class = "hint", "Outputs reflect the currently selected dataset."))
        ),
        column(9,
               uiOutput("validation_note"),
               div(class = "panel-card", h4("Performance on unseen observations"),
                   tableOutput("comparison_table"),
                   div(class = "hint", "Higher AUC, sensitivity and specificity are better; ",
                       "lower Brier and log loss are better. Sensitivity and specificity ",
                       "change with the threshold; AUC and scoring rules do not.")),
               fluidRow(
                 column(6, div(class = "panel-card", plotOutput("roc_plot", height = "305px"))),
                 column(6, div(class = "panel-card", plotOutput("calibration_plot", height = "305px")))
               ),
               fluidRow(
                 column(6, div(class = "panel-card", plotOutput("split_plot", height = "320px"))),
                 column(6, div(class = "panel-card", plotOutput("residual_plot", height = "320px")))
               ),
               fluidRow(
                 column(6, div(class = "panel-card",
                               plotOutput("glm_gam_difference_plot", height = "335px"),
                               textOutput("glm_gam_difference_stats"),
                               p(class = "hint", "Signed differences in occurrence probability, ",
                                 "not residuals or uncertainty. For a fair comparison of ",
                                 "model structure, choose the same predictors in GLM and GAM. ",
                                 "Without an uploaded prediction grid, only surveyed locations ",
                                 "are shown."))),
                 column(6, div(class = "panel-card",
                               plotOutput("disagreement_plot", height = "335px"),
                               p(class = "hint", "Absolute spread across GLM, GAM and ",
                                 "Random Forest; it does not indicate which model is right.")))
               ),
               div(class = "panel-card", h4("Location ranking / restoration shortlist"),
                   textOutput("candidate_note"), tableOutput("candidate_table"),
                   p(class = "hint", "High predicted occurrence does not establish habitat ",
                     "restorability, feasibility, connectivity or colonisation potential.")),
               conditionalPanel(
                 condition = "(input.data_source == 'simulation' || input.data_source == 'custom') && input.reveal_truth == true",
                 div(class = "panel-card", plotOutput("truth_plot", height = "325px"),
                     p(class = "hint", "True occurrence probabilities are available ",
                       "only for simulated datasets, never for uploaded field data."))))
        )
      )
  ),
  div(class = "footer-note",
      tags$hr(),
      tags$p(style = "margin: 2px 0;", tags$strong("EcoSDM Lab")," | Developed by Dr Oscar Rodriguez de Rivera"),
      tags$p(style = "margin: 2px 0;","Senior Lecturer in Statistics | University of Exeter"),
      tags$p(style = "margin: 2px 0;","Contact: ",tags$a(href = "mailto:o.m.rodriguez-de-rivera-ortega@exeter.ac.uk","Email")),
      tags$p(style = "margin: 2px 0;","© 2026 | Species Distribution Modelling ","Interactive Teaching Resource"))
)

# -------------------------------- Server ------------------------------------
server <- function(input, output, session) {
  uploaded_surveys <- reactive({
    req(input$survey_upload)
    read_upload(input$survey_upload)
  })
  uploaded_grid <- reactive({
    if (is.null(input$grid_upload)) return(NULL)
    read_upload(input$grid_upload)
  })

  # Runs once initially and thereafter only when Generate is clicked. Slider
  # changes do not unexpectedly invalidate every model during an explanation.
  own_species <- eventReactive(input$generate_species, {
    validate(need(is.finite(input$sim_seed) && input$sim_seed >= 1 &&
                    input$sim_seed <= 999999,
                  "Please enter a random seed between 1 and 999999."))
    simulate_own_species(
      sim_landscape, input$sim_n, input$sim_grassland,
      input$sim_temperature, input$sim_temp_width,
      input$sim_rainfall, input$sim_rain_width,
      input$sim_prevalence, input$sim_sampling, input$sim_seed)
  }, ignoreNULL = FALSE)

  # Set up column selectors when a new file is loaded; students can override.
  observeEvent(uploaded_surveys(), {
    d <- uploaded_surveys()
    cols <- names(d)
    numeric_cols <- cols[vapply(d, is.numeric, logical(1))]
    guess <- function(options, possible, fallback) {
      x <- intersect(possible, options)
      if (length(x)) x[1] else fallback
    }
    binary_cols <- cols[vapply(d, function(x) {
      v <- unique(stats::na.omit(suppressWarnings(as.numeric(as.character(x)))))
      length(v) > 0 && all(v %in% c(0, 1))
    }, logical(1))]
    resp <- guess(cols, c("presence", "occurrence", "detected", "pa", "observed"),
                  if (length(binary_cols)) binary_cols[1] else cols[1])
    xc <- guess(numeric_cols, c("x_km", "x", "easting", "Easting", "lon", "longitude"),
                if (length(numeric_cols)) numeric_cols[1] else cols[1])
    yc <- guess(numeric_cols, c("y_km", "y", "northing", "Northing", "lat", "latitude"),
                if (length(numeric_cols) > 1) numeric_cols[2] else cols[1])
    updateSelectInput(session, "response_col", choices = cols, selected = resp)
    updateSelectInput(session, "x_col", choices = numeric_cols, selected = xc)
    updateSelectInput(session, "y_col", choices = numeric_cols, selected = yc)
  })

  dataset <- reactive({
    if (input$data_source == "simulation") {
      prepare_data(sim_observations, sim_landscape,
                   response_col = "presence", x_col = "x_km", y_col = "y_km",
                   simulated = TRUE)
    } else if (input$data_source == "custom") {
      created <- own_species()
      prepare_data(created$observations, created$landscape,
                   response_col = "presence", x_col = "x_km", y_col = "y_km",
                   simulated = TRUE)
    } else {
      validate(need(!is.null(input$survey_upload),
                    "Upload a survey CSV under Data & upload to begin modelling."))
      req(input$response_col, input$x_col, input$y_col)
      prepare_data(uploaded_surveys(), uploaded_grid(),
                   response_col = input$response_col,
                   x_col = input$x_col, y_col = input$y_col)
    }
  })

  # Replace dropdown values for the new dataset. All tabs use the same columns.
  observeEvent(dataset(), {
    ds <- dataset()
    choice <- stats::setNames(ds$predictors, unname(ds$labels[ds$predictors]))
    fallback <- head(ds$predictors, 3L)
    for (type in names(model_names)) {
      original <- isolate(input[[paste0("vars_", type)]])
      selected <- intersect(original, ds$predictors)
      if (!length(selected)) selected <- fallback
      updateCheckboxGroupInput(session, paste0("vars_", type),
                               choices = choice, selected = selected)
      current <- isolate(input[[paste0("effect_", type)]])
      updateSelectInput(session, paste0("effect_", type),
                        choices = choice,
                        selected = if (current %in% ds$predictors) current else ds$predictors[1])
    }
    current_data_var <- isolate(input$data_plot_var)
    updateSelectInput(session, "data_plot_var", choices = choice,
                      selected = if (current_data_var %in% ds$predictors)
                        current_data_var else ds$predictors[1])
    current_cor <- intersect(isolate(input$cor_vars), ds$predictors)
    if (length(current_cor) < 2) current_cor <- head(ds$predictors, 6L)
    updateSelectizeInput(session, "cor_vars", choices = choice,
                         selected = head(current_cor, 12L), server = TRUE)
  })

  split <- reactive(make_split(input$validation, dataset()$observations))
  train <- reactive(split()$train)
  test <- reactive(split()$test)

  selected_predictors <- function(type) {
    ds <- dataset()
    selected <- intersect(input[[paste0("vars_", type)]], ds$predictors)
    validate(need(length(selected) >= 1L,
                  "Please select at least one environmental predictor for this model."))
    selected
  }

  # Each model's settings are captured only when Fit is clicked. Storing the
  # active dataset/validation in the snapshot prevents a model fitted on an old
  # dataset or split from being displayed as though it used the new settings.
  model_spec <- function(type) {
    req(input$validation)
    spec <- list(data = dataset(), validation = input$validation,
                 vars = selected_predictors(type))
    if (type == "gam") {
      req(input$gam_k)
      spec$gam_k <- input$gam_k
      spec$smooth_vars <- as.character(intersect(input$gam_smooth_vars, spec$vars))
    } else if (type == "rf") {
      req(input$rf_trees, input$rf_leaf)
      spec$trees <- input$rf_trees
      spec$leaf <- input$rf_leaf
    }
    spec
  }

  fit_snapshot <- function(type, spec) {
    # Never refit in renderPlot / renderTable. A single explicit click creates
    # one model that is shared by all outputs for this Shiny session.
    split_now <- make_split(spec$validation, spec$data$observations)
    fit <- withProgress({
      if (type == "glm") {
        fit_sdm("glm", spec$vars, split_now$train)
      } else if (type == "gam") {
        fit_sdm("gam", spec$vars, split_now$train,
                gam_k = spec$gam_k, gam_smooth_vars = spec$smooth_vars)
      } else {
        fit_sdm("rf", spec$vars, split_now$train,
                trees = spec$trees, leaf = spec$leaf)
      }
    }, message = paste("Fitting", model_names[[type]]), value = 0)
    list(model = fit, spec = spec)
  }

  # An initial fit is performed once per model on its first use, so the
  # workshop opens with results. Subsequent data/predictor/validation changes
  # never trigger fitting: students must use Fit or Fit all models.
  glm_snapshot <- eventReactive(list(input$fit_glm, input$fit_all_models), {
    fit_snapshot("glm", model_spec("glm"))
  }, ignoreNULL = FALSE)
  gam_snapshot <- eventReactive(list(input$fit_gam, input$fit_all_models), {
    fit_snapshot("gam", model_spec("gam"))
  }, ignoreNULL = FALSE)
  rf_snapshot <- eventReactive(list(input$fit_rf, input$fit_all_models), {
    fit_snapshot("rf", model_spec("rf"))
  }, ignoreNULL = FALSE)

  get_snapshot <- function(type) {
    switch(type, glm = glm_snapshot(), gam = gam_snapshot(), rf = rf_snapshot())
  }
  is_current <- function(type) {
    identical(get_snapshot(type)$spec, model_spec(type))
  }
  get_model <- function(type) {
    validate(need(is_current(type),
                  paste("Model settings or data have changed. Please click Fit",
                        model_names[[type]], "(or Fit all models in Compare & diagnose).")))
    get_snapshot(type)$model
  }

  for (type in names(model_names)) {
    local({
      m <- type
      output[[paste0("fit_status_", m)]] <- renderUI({
        if (is_current(m)) {
          div(class = "hint fit-ready", "Model fitted with current settings.")
        } else {
          div(class = "hint fit-stale", "Settings changed. Click Fit to update results.")
        }
      })
    })
  }
  output$compare_fit_status <- renderUI({
    ready <- vapply(names(model_names), is_current, logical(1))
    if (all(ready)) {
      div(class = "hint fit-ready", "All three models use the current settings and validation split.")
    } else {
      div(class = "hint fit-stale", "One or more models need updating. Click Fit all three models.")
    }
  })

  test_predictions <- reactive({
    d <- test()
    list(glm = predict_sdm(get_model("glm"), d),
         gam = predict_sdm(get_model("gam"), d),
         rf = predict_sdm(get_model("rf"), d))
  })
  grid_predictions <- reactive({
    d <- dataset()$landscape
    list(glm = predict_sdm(get_model("glm"), d),
         gam = predict_sdm(get_model("gam"), d),
         rf = predict_sdm(get_model("rf"), d))
  })

  # Intro tab: feedback, preview and environmental variable vs occurrence.
  output$data_status <- renderUI({
    if (input$data_source == "upload" && is.null(input$survey_upload)) {
      return(div(class = "panel-card", h4("No upload yet"),
                 p("Upload a survey CSV, then choose the response and coordinate columns.")))
    }
    ds <- dataset()
    div(class = "stat-row",
        div(class = "stat-pill", div(class = "value", nrow(ds$observations)),
            div(class = "label", "Complete surveys")),
        div(class = "stat-pill", div(class = "value", sum(ds$observations$presence == 1)),
            div(class = "label", "Observed presences")),
        div(class = "stat-pill", div(class = "value", length(ds$predictors)),
            div(class = "label", "Numeric predictors")),
        div(class = "stat-pill", div(class = "value", nrow(ds$landscape)),
            div(class = "label", if (ds$has_grid) "Prediction locations" else "Survey locations only")),
        if (ds$n_survey_dropped > 0L || ds$n_grid_dropped > 0L)
          div(class = "stat-pill", div(class = "value", ds$n_survey_dropped),
              div(class = "label", paste0("Survey rows dropped; prediction rows dropped: ",
                                          ds$n_grid_dropped))))
  })

  output$data_generation_description <- renderUI({
    if (input$data_source == "upload") {
      return(tagList(
        p("These records were supplied by you. The app does not assume ",
          "how your observations were collected."),
        p("Presence/absence must be 0 or 1. Environmental layers and coordinates ",
          "should describe each survey. A second location CSV is optional."),
        p(class = "hint", "Uploaded data have no known true probabilities. ",
          "Only use approved, non-sensitive ecological data.")))
    }
    if (input$data_source == "custom") {
      sp <- own_species()
      return(tagList(
        p("You created a fictional species across 3,640 landscape cells, ",
          "using selected temperature/rainfall optima, niche widths and ",
          "a grassland association."),
        tags$ol(
          tags$li("An ecological response assigns an occurrence probability ",
                  "to each landscape cell."),
          tags$li("Survey locations are sampled uniformly or closer to roads."),
          tags$li("Each location gets a 0/1 presence draw with perfect detection.")),
        p(sprintf("Mean true occurrence: %.1f%% across the landscape. Observed survey presence: %.1f%%.",
                  100 * sp$expected_prevalence, 100 * sp$observed_prevalence)),
        p(class = "hint", "This is a generated species, not a real Cornwall species. ",
          "True suitability can be revealed in Compare & diagnose.")))
    }
    tagList(
      p("A fictional butterfly was surveyed at 750 sites on a 3,640-cell landscape."),
      tags$ol(
        tags$li("Generate spatially varying rainfall, temperature, ",
                "elevation, grassland and road-distance layers."),
        tags$li("Define curved temperature/rainfall responses, grassland effects ",
                "and an environmental interaction."),
        tags$li("Sample presence/absence at 750 locations from independent Bernoulli ",
                "distributions, assuming perfect detection.")),
      p(class = "hint", "The geographic landscape is fictional. True occurrence ",
        "probabilities can be revealed in Compare & diagnose."))
  })

  output$correlation_plot <- renderPlot({
    correlation_plot(dataset(), input$cor_vars, input$cor_method)
  }, res = 100)

  output$data_preview <- renderTable({
    ds <- dataset()
    columns <- unique(c("presence", "x_km", "y_km", ds$predictors))
    head(ds$observations[, columns, drop = FALSE], 8L)
  }, striped = TRUE, spacing = "xs", digits = 3)

  output$data_columns <- renderTable({
    ds <- dataset()
    d <- ds$observations
    data.frame(
      Predictor = unname(ds$labels[ds$predictors]),
      Minimum = round(vapply(ds$predictors, function(v) min(d[[v]]), numeric(1)), 2),
      Median = round(vapply(ds$predictors, function(v) stats::median(d[[v]]), numeric(1)), 2),
      Maximum = round(vapply(ds$predictors, function(v) max(d[[v]]), numeric(1)), 2),
      row.names = NULL
    )
  }, striped = TRUE, spacing = "xs")

  output$data_scatter <- renderPlot({
    ds <- dataset()
    v <- input$data_plot_var
    validate(need(v %in% ds$predictors, "Choose an available environmental variable."))
    d <- ds$observations
    ggplot(d, aes(x = .data[[v]], y = factor(presence))) +
      geom_jitter(height = 0.14, width = 0, alpha = 0.28,
                  colour = "#27684e", size = 1.8) +
      scale_y_discrete(labels = c("0" = "Absent (0)", "1" = "Present (1)")) +
      labs(title = "Environmental conditions at survey sites",
           subtitle = "Each dot is a survey; vertical jitter separates overlapping records",
           x = ds$labels[[v]], y = "Survey outcome") + plot_theme()
  }, res = 100)

  output$download_survey_template <- downloadHandler(
    filename = function() "ecosdm_example_surveys.csv",
    content = function(file) write.csv(sim_observations, file, row.names = FALSE)
  )
  output$download_grid_template <- downloadHandler(
    filename = function() "ecosdm_example_prediction_grid.csv",
    content = function(file) write.csv(
      sim_landscape[, setdiff(names(sim_landscape), "true_probability")],
      file, row.names = FALSE)
  )

  # P-value-based filtering is exploratory, never confirmatory variable selection.
  for (type in names(model_names)) {
    local({
      m <- type
      observeEvent(input[[paste0("restore_all_", m)]], {
        ds <- dataset()
        choice <- stats::setNames(ds$predictors, unname(ds$labels[ds$predictors]))
        updateCheckboxGroupInput(session, paste0("vars_", m),
                                 choices = choice, selected = ds$predictors)
      })
      if (m %in% c("glm", "gam")) {
        observeEvent(input[[paste0("keep_significant_", m)]], {
          fit <- get_model(m)
          keep <- significant_predictors(fit)
          if (!length(keep)) {
            showNotification("No predictors currently have p < 0.05; ",
                             type = "warning", duration = 5)
          } else {
            updateCheckboxGroupInput(session, paste0("vars_", m),
                                     selected = keep)
            showNotification("Selection updated. Click Fit to refit this model. ",
                             type = "message", duration = 4)
          }
        })
      }
    })
  }

  # Keep the smooth-term choices in sync with currently selected GAM predictors.
  # An empty selection is intentional: it represents a fully linear GAM.
  observeEvent(input$vars_gam, {
    ds <- dataset()
    chosen <- intersect(input$vars_gam, ds$predictors)
    choice <- stats::setNames(chosen, unname(ds$labels[chosen]))
    keep <- intersect(input$gam_smooth_vars, chosen)
    updateCheckboxGroupInput(session, "gam_smooth_vars", choices = choice,
                             selected = keep)
  }, ignoreInit = TRUE)

  output$gam_partial_plot <- renderPlot({
    ds <- dataset()
    variable <- input$effect_gam
    plot_gam_partial(get_model("gam"), variable, train(), ds)
  }, res = 100)

  # GLM, GAM and Random Forest tabs.
  for (type in names(model_names)) {
    local({
      m <- type
      output[[paste0("map_", m)]] <- renderPlot({
        plot_prediction(get_model(m), input[[paste0("show_obs_", m)]], train(), dataset())
      }, res = 100)
      output[[paste0("response_", m)]] <- renderPlot({
        ds <- dataset()
        variable <- input[[paste0("effect_", m)]]
        validate(need(variable %in% ds$predictors, "Choose an available variable."))
        fitted <- get_model(m)
        # Do not silently re-fit a GLM inside a graph. Overlay the already
        # fitted GLM only when it uses exactly the same predictors as the GAM.
        baseline <- NULL
        if (m == "gam") {
          baseline <- tryCatch(get_model("glm"), error = function(e) NULL)
          if (!is.null(baseline) && !setequal(baseline$vars, fitted$vars))
            baseline <- NULL
        }
        plot_response(fitted, variable, train(), ds, baseline = baseline)
      }, res = 100)
      output[[paste0("formula_", m)]] <- renderText(model_formulation(get_model(m)))
      output[[paste0("summary_", m)]] <- renderTable({
        model_explanation(get_model(m), dataset()$labels)
      }, striped = TRUE, spacing = "xs", digits = 3)
      output[[paste0("stats_", m)]] <- renderUI({
        mo <- get_model(m)
        sc <- scores(test()$presence, predict_sdm(mo, test()), input$threshold)
        fmt <- function(x) ifelse(is.na(x), "N/A", sprintf("%.3f", x))
        div(class = "stat-row",
            div(class = "stat-pill", div(class = "value", fmt(sc["AUC"])),
                div(class = "label", "Test AUC ↑")),
            div(class = "stat-pill", div(class = "value", fmt(sc["Brier"])),
                div(class = "label", "Test Brier ↓")),
            div(class = "stat-pill", div(class = "value", nrow(train())),
                div(class = "label", "Training surveys")),
            div(class = "stat-pill", div(class = "value", nrow(test())),
                div(class = "label", "Unseen surveys")))
      })
    })
  }

  output$validation_note <- renderUI({
    y <- test()$presence
    if (length(unique(y)) < 2L) {
      div(class = "panel-card learning-card",
          tags$b("Holdout warning: "),
          "The spatial test sample contains only one outcome class. ",
          "AUC/ROC and one of sensitivity/specificity are unavailable. ",
          "Try random validation or a more balanced spatial dataset.")
    }
  })

  output$comparison_table <- renderTable({
    preds <- test_predictions()
    fmt <- function(x) ifelse(is.na(x), "N/A", sprintf("%.3f", x))
    do.call(rbind, lapply(names(preds), function(m) {
      z <- scores(test()$presence, preds[[m]], input$threshold)
      data.frame(Model = model_names[[m]],
                 AUC = fmt(z["AUC"]), Brier = fmt(z["Brier"]),
                 `Log loss` = fmt(z["LogLoss"]),
                 Sensitivity = ifelse(is.na(z["Sensitivity"]), "N/A",
                                      sprintf("%.1f%%", 100 * z["Sensitivity"])),
                 Specificity = ifelse(is.na(z["Specificity"]), "N/A",
                                      sprintf("%.1f%%", 100 * z["Specificity"])),
                 check.names = FALSE)
    }))
  }, striped = TRUE, spacing = "s")

  output$roc_plot <- renderPlot({
    preds <- test_predictions()
    validate(need(length(unique(test()$presence)) == 2L,
                  "ROC requires both presences and absences in the test set."))
    d <- do.call(rbind, lapply(names(preds), function(m)
      roc_points(test()$presence, preds[[m]], model_names[[m]])))
    ggplot(d, aes(fpr, tpr, colour = model)) +
      geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "#b5c2b9") +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = stats::setNames(unname(model_colors), unname(model_names))) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = "ROC curves", x = "False-positive rate", y = "True-positive rate",
           colour = NULL) + plot_theme() + theme(legend.position = "bottom")
  }, res = 100)

  output$calibration_plot <- renderPlot({
    preds <- test_predictions()
    d <- do.call(rbind, lapply(names(preds), function(m)
      calibration_points(test()$presence, preds[[m]], model_names[[m]])))
    ggplot(d, aes(predicted, observed, colour = model)) +
      geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "#b5c2b9") +
      geom_line(linewidth = 0.8) + geom_point(aes(size = n), alpha = 0.85) +
      scale_size_continuous(range = c(1.5, 4.5)) +
      scale_colour_manual(values = stats::setNames(unname(model_colors), unname(model_names))) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = "Calibration by probability bin", x = "Mean predicted probability",
           y = "Observed presence fraction", colour = NULL, size = "Sites") +
      plot_theme() + theme(legend.position = "bottom")
  }, res = 100)

  output$split_plot <- renderPlot({
    ds <- dataset()
    d <- rbind(transform(train(), Partition = "Training"),
               transform(test(), Partition = "Test"))
    if (ds$simulated) {
      p <- ggplot(ds$landscape, aes(x_km, y_km)) +
        geom_raster(aes(fill = grassland)) +
        scale_fill_gradient(low = "#edf0df", high = "#71a580", name = "Grassland (%)")
    } else {
      p <- ggplot()
      if (ds$has_grid) {
        p <- p + geom_point(data = ds$landscape, aes(x_km, y_km),
                            colour = "#dde9df", shape = 15, size = 1, alpha = 0.45)
      }
    }
    p + geom_point(data = d, aes(x_km, y_km, colour = Partition),
                   inherit.aes = FALSE, size = 1.6, alpha = 0.85) +
      scale_colour_manual(values = c(Training = "#277b57", Test = "#bd4f26")) +
      coord_equal(expand = FALSE) +
      labs(title = "Survey locations & held-out sites",
           x = if (ds$simulated) "Fictional km east" else "X coordinate",
           y = if (ds$simulated) "Fictional km north" else "Y coordinate") +
      plot_theme() + theme(panel.grid = element_blank(), legend.position = "bottom")
  }, res = 100)

  output$residual_plot <- renderPlot({
    ds <- dataset()
    d <- test()
    d$residual <- d$presence - test_predictions()[[input$res_model]]
    ggplot(d, aes(x_km, y_km, colour = residual)) + geom_point(size = 2.3, alpha = 0.9) +
      scale_colour_gradient2(low = "#bd5c48", mid = "#f4eee1", high = "#277b57",
                             midpoint = 0, limits = c(-1, 1), name = "Observed\n− predicted") +
      coord_equal(xlim = range(ds$landscape$x_km), ylim = range(ds$landscape$y_km)) +
      labs(title = paste(model_names[[input$res_model]], "· test residuals"),
           x = if (ds$simulated) "Fictional km east" else "X coordinate",
           y = if (ds$simulated) "Fictional km north" else "Y coordinate") + plot_theme()
  }, res = 100)

  # This map uses all available prediction locations, while test residuals above
  # continue to use only held-out observations.
  output$glm_gam_difference_plot <- renderPlot({
    plot_glm_gam_difference(dataset(), grid_predictions())
  }, res = 100)

  output$glm_gam_difference_stats <- renderText({
    delta <- glm_gam_difference(grid_predictions())
    sprintf("Across %s prediction locations: mean signed difference %+.3f; mean absolute difference %.3f; maximum absolute difference %.3f.",
            format(length(delta), big.mark = ",", scientific = FALSE),
            mean(delta), mean(abs(delta)), max(abs(delta)))
  })

  output$disagreement_plot <- renderPlot({
    ds <- dataset()
    preds <- grid_predictions()
    d <- ds$landscape
    d$disagreement <- do.call(pmax, preds) - do.call(pmin, preds)
    p <- spatial_layer(d, "disagreement", ds$simulated, limits = c(0, 1),
                       legend_title = "Max − min\nprobability")
    if (ds$simulated) {
      candidate_sites <- get_candidate_sites(ds$landscape)
      p <- p + geom_point(data = candidate_sites, aes(x_km, y_km),
                          inherit.aes = FALSE, shape = 21, fill = "white",
                          colour = "#193f34", size = 5, stroke = 1) +
        geom_text(data = candidate_sites, aes(x_km, y_km, label = short_label),
                  inherit.aes = FALSE, colour = "#193f34", fontface = "bold", size = 3)
    }
    p + labs(title = "Where do the three models disagree?",
             subtitle = "Prediction differences are NOT confidence intervals or posterior uncertainty",
             x = if (ds$simulated) "Fictional km east" else "X coordinate",
             y = if (ds$simulated) "Fictional km north" else "Y coordinate")
  }, res = 100)

  output$candidate_note <- renderText({
    if (dataset()$simulated) {
      "Three fictional restoration candidates. Which would you investigate first?"
    } else {
      paste0("Three locations with the highest mean prediction among the uploaded locations. ",
             "These are NOT independently identified restoration sites.")
    }
  })

  output$candidate_table <- renderTable({
    ds <- dataset()
    preds <- grid_predictions()
    if (ds$simulated) {
      candidate_sites <- get_candidate_sites(ds$landscape)
      ids <- match(paste(candidate_sites$x_km, candidate_sites$y_km),
                   paste(ds$landscape$x_km, ds$landscape$y_km))
      d <- data.frame(Site = candidate_sites$site,
                      `Grassland (%)` = sprintf("%.0f", candidate_sites$grassland),
                      `GLM probability` = sprintf("%.2f", preds$glm[ids]),
                      `GAM probability` = sprintf("%.2f", preds$gam[ids]),
                      `RF probability` = sprintf("%.2f", preds$rf[ids]),
                      check.names = FALSE)
      if (isTRUE(input$reveal_truth)) {
        d[["True (simulated) probability"]] <-
          sprintf("%.2f", candidate_sites$true_probability)
      }
      d
    } else {
      av <- (preds$glm + preds$gam + preds$rf) / 3
      ids <- head(order(av, decreasing = TRUE), 3L)
      data.frame(`Location row` = ids,
                 X = round(ds$landscape$x_km[ids], 3),
                 Y = round(ds$landscape$y_km[ids], 3),
                 `Mean probability` = sprintf("%.2f", av[ids]),
                 GLM = sprintf("%.2f", preds$glm[ids]),
                 GAM = sprintf("%.2f", preds$gam[ids]),
                 RF = sprintf("%.2f", preds$rf[ids]), check.names = FALSE)
    }
  }, striped = TRUE, spacing = "s")

  output$truth_plot <- renderPlot({
    req(dataset()$simulated)
    spatial_layer(dataset()$landscape, "true_probability", TRUE, c(0, 1),
                  "True simulated\noccurrence") +
      labs(title = "Simulation truth (unobservable in real studies)",
           x = "Fictional km east", y = "Fictional km north")
  }, res = 100)

  output$download_observations <- downloadHandler(
    filename = function() "ecosdm_active_surveys.csv",
    content = function(file) write.csv(dataset()$observations, file, row.names = FALSE)
  )
  output$download_predictions <- downloadHandler(
    filename = function() "ecosdm_active_model_predictions.csv",
    content = function(file) {
      ds <- dataset()
      d <- ds$landscape[, c("x_km", "y_km", ds$predictors), drop = FALSE]
      p <- grid_predictions()
      for (m in names(p)) d[[paste0("prediction_", m)]] <- p[[m]]
      write.csv(d, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
#' @title Get the correct model
#' @description Get the correct model
#' @param x data
#' @param method The type of analysis. Possible analysis: c('time_series')
#' @keywords vfmodels
#' @author Emiel Veersma
#' @family get_model
#' @examples
#' model <- get_model(lh, 'time_series')
#' @export
get_model <- function(x, method, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc <- mc[!names(mc) %in% 'method']
  mc[[1]] <- as.name(paste0("get_model.", method))
  res <- eval(mc, parent.frame())
  res$call <- sys.calls()[[1]]
  return(res)
}

#' @title Get Time Series Model
#' @description An algorithm to find the appropriate time series model. The algorithm performs a forward search on the formula.
#' As a starting point, a completely empty formula is chosen and every time a sum or product is added/removed. The local optimum is
#' found by minimizing the Bayesian information criterion-value (BIC-value).
#'
#' Some improvements still need to be made:
#' \itemize{
#' \item Adjust for heterskedasticity
#' \item Implement smoothing (now only sums and products of constants, linear relations and lags are incorporated)
#' \item Check for overfitting
#' }
#' @param x a univariate time series
#' @param lag the number of lags to incorporate in the model
#' @param labels the dates of the series
#' @keywords vfmodels time series
#' @author Emiel Veersma
#' @seealso \url{http://www.automaticstatistician.com/research.php} for the inspiration of this function.
#' @family get_model
#' @examples
#' model <- get_model(lh, 'time_series')
#' @export
get_model.time_series <- function(x, lag = NULL, labels = NULL, family = gaussian, train = NULL) {
  if (is.null(lag)) {
    lag <- 3
    warning('Lag is not set and therefore defaulted to 3. See the help of get_model.time_series for more information.')
  }
  built <- build_x(x, lag, labels)
  data <- x <- built$x
  if (!is.null(train)) {
    x <- x[train, ]
  }
  model <- init_model(x, family)
  optimal_model <- model
  good_sol <- TRUE
  for (i in 1:200) {
    old_BIC <- model$BIC
    new_BIC <- model$BIC - 1
    while (old_BIC > new_BIC & good_sol) {
      old_BIC <- model$BIC
      model <- local_search(x, model, lag, built$periods, family)
      new_BIC <- model$BIC
    }
    if (model$BIC < optimal_model$BIC) {
      optimal_model <- model
    }
    model <- shake(x, optimal_model, lag, built$periods, family = family)
    if (model$BIC < (optimal_model$BIC * 1.1)) {
      good_sol <- FALSE
    }
  }
  res <- build_results(x, optimal_model, built$periods, data, lag, train)
  return(res)
}

#' @export
plot.model <- function(model) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(paste0('plot_model.', model$method))
  eval(mc, parent.frame())
}

#' @export
summary.model <- function(model) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(paste0('summary_model.', model$method))
  eval(mc, parent.frame())
}

#' @export
summary_model.time_series <- function(model) {
  res <- summary(model$lm)
  res$periods <- model$periods
  return(res)
}

#' @export
predict.model <- function(model, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(paste0('predict_model.', model$method))
  eval(mc, parent.frame())
}

predict_one_series <- function(model, data, prob) {
  prediction <- predict(model$lm, newdata = data, se.fit = TRUE)
  family <- family(model$lm)
  res <- family$linkinv(prediction$fit + prediction$se.fit %*% t(qnorm(prob)))
  return(res)
}

#' @export
predict_model.time_series <- function(model, x, n.ahead = 1, labels = NULL, train = NULL) {
  all_labels <- labels
  if (!is.null(train)) {
    x <- x[train]
    labels <- labels[train]
  } else {
    train <- 1:length(x)
  }
  for (i in 1:n.ahead) {
    labels <- c(labels, all_labels[max(train) + i])
    model$periods[length(model$periods)] <- model$periods[length(model$periods)] + 1
    if (i == 1) {
      x <- c(x, NA)
      built <- build_x(x, model$lag, labels, periods = model$periods)
      res <- list()
      res[[i]] <- tail(predict_one_series(model, data = built$x, prob = (1:19)/20), 1)
    } else if (i > 1) {
      res[[i]] <- list()
      x <- c(x, NA)
      for (j in res[[i - 1]]) {
        x[(length(x) - i + 1):(length(x) - 1)] <- j
        built <- build_x(x, model$lag, labels, periods = model$periods)
        res[[i]] <- c(res[[i]], as.data.frame(tail(predict_one_series(model, data = built$x, prob = (1:19) / 20), i)))
      }
    }
  }
  res <- res[[length(res)]]
  if (class(res) == 'list') {
    res <- do.call(rbind, res)
  }
  if (i == 1) {
    outcome <- res
  } else {
    outcome <- NULL
    for (i in 1:n.ahead) {
      col <- res[, i][order(res[, i])]
      col <- col[(1:19) ^ n.ahead]
      outcome <- rbind(outcome, col)
    }
  }
  outcome <- as.data.frame(outcome)
  row.names(outcome) <- tail(labels, n.ahead)
  colnames(outcome) <- (1:19) / 20
  return(outcome)
}

#' @import ggplot2
#' @export
plot_model.time_series <- function(model) {
  periods <- substr(names(model$x), 1, 6) == 'period'
  model$x$periods <- factor(as.matrix(model$x[, periods]) %*% (1:sum(periods)))
  par(mfrow = c(2, 1))
  ggplot() + geom_line(data = model$x, aes(x = labels, y = x, colour = periods)) + xlab('date') + ylab('value') + theme_minimal()
  ggplot() + geom_point(aes(x = model$x$labels, y = model$lm$residuals, colour = model$x$periods)) + xlab('date') + ylab('residuals') + theme_minimal()
}

shake_reduce <- function(x, model, lag, periods, family) {
  i <- 1
  while (i < length(model$form)) {
    i <- i + 1
    if (rnorm(1) > 0.2) {
      mode <- NULL
      if (model$form[i] == '+ 1') {
        mode <- 'const'
      } else if (substr(model$form[i], nchar(model$form[i]), nchar(model$form[i])) == '(') {
        mode <- 'prod'
      } else if (substr(model$form[i], 1, 1) == '+') {
        mode <- 'sum'
      }
      if (!is.null(mode)) {
        model <- remove_element(i = i, x = x, form = model$form, mode = mode, family = family)
      }
    }
  }
  return(model)
}

shake_expand <- function(x, model, lag, periods, family) {
  variants <- build_variants(lag, periods)
  prob <- rnorm(length(variants))
  for (i in 1:length(prob)) {
    if (prob[i] > 1) {
      mode = 'sum'
      model <- add_element(addition = variants[i], x = x, form = model$form, mode = mode, family = family)
    } else if (prob[i] < -2) {
      mode = 'prod'
      model <- add_element(addition = variants[i], x = x, form = model$form, mode = mode, family = family)
    }
  }
  return(model)
}

shake <- function(x, model, lag, periods, family) {
  model <- shake_reduce(x, model, lag, periods, family)
  model <- shake_expand(x, model, lag, periods, family)
  return(model)
}

build_results <- function(x, model, periods, data, lag, train) {
  res <- model
  res$method <- 'time_series'
  res$periods <- periods
  res$x <- x
  res$lag <- lag
  res$train <- train
  class(res) <- 'model'
  return(res)
}

#' @importFrom strucchange breakpoints
build_x <- function(x, lag, labels = NULL, periods = NULL) {
  if (is.null(periods)) {
    periods <- breakpoints(x ~ 1)$breakpoints
  }
  len_x <- length(x)
  lag_1 <- c(NA, head(x, n = len_x - 1))
  linear <- 1:len_x
  if (!is.null(labels)) {
    x <- x[order(labels)]
  }
  x <- as.data.frame(cbind(x, linear))
  if (!is.null(labels)) {
    x$labels <- labels[order(labels)]
  } else {
    x$labels <- 1:nrow(x)
  }
  if (length(periods) > 0) {
    periods <- unique(c(0, periods, nrow(x)))
    rows <- 1:nrow(x)
    periods_set <- 1:(length(periods) - 1)
    names <- paste0('period_', periods_set)
    for (i in periods_set) {
      x[, names[i]] <- rows > periods[i] & rows <= periods[i + 1]
    }
  }
  lag_set <- 1:lag
  names <- paste0('lag_', lag_set)
  for (i in lag_set) {
    x[, names[i]] <- c(rep(NA, i), head(x$x, n = len_x - i))
  }
  x <- x[(lag + 1):nrow(x), ]
  return(list(x = x, periods = periods))
}

build_variants <- function(lag, periods) {
  variants <- c('Constant', 'Linear', paste('Periodic', 1:lag), paste('Periods', 1:(length(periods) - 1))) #White Noise and Smooth still need to be done
  return(variants)
}

local_search <- function(x, init_model, lag, periods, family) {
  variants <- build_variants(lag, periods)
  increased_models <- increase_models(x, init_model$form, variants, family)
  reduced_models <- reduce_models(x, init_model$form, family)
  models <- c(increased_models, reduced_models)
  model_overview <- do.call(Map, c(c, models))
  best <- which.min(model_overview$BIC)
  model <- models[[best]]
  model$form <- simplify(model$form)
  if (model$BIC >= init_model$BIC) {
    model <- init_model
  }
  return(model)
}

simplify <- function(form) {
  const <- which(form %in% c('- 1', '+ 1'))
  if (length(const) > 0) {
    form <- c(form[1], form[const], form[-c(1, const)])
    end_brack <- which(form == ')')
    if (length(end_brack) > 0 && substr(form[end_brack - 1], nchar(form[end_brack - 1]), nchar(form[end_brack - 1])) == '(') {
      form[end_brack - 1] <- substr(form[end_brack - 1], nchar(form[end_brack - 1]) - 3, nchar(form[end_brack - 1]))
      form <- form[-end_brack]
    }
  }
  return(form)
}

reduce_models <- function(x, form, family) {
  if (any(!form %in% c("x ~", '- 1'))) {
    models_sum <- models_prod <- models_const <- NULL
    rem_const <- which(form == '+ 1')
    rem_sum <- which(substr(form, 1, 1) == '+')
    rem_prod <- which(substr(form, nchar(form), nchar(form)) == '(')
    rem_sum <- setdiff(rem_sum, rem_prod)
    rem_sum <- setdiff(rem_sum, rem_const)
    if (length(rem_const) > 0) {
      models_const <- lapply(rem_const, remove_element, x = x, form = form, mode = 'const', family = family)
    }
    if (length(rem_sum) > 0) {
      models_sum <- lapply(rem_sum, remove_element, x = x, form = form, mode = 'sum', family = family)
    }
    if (length(rem_prod) > 0) {
      models_prod <- lapply(rem_prod, remove_element, x = x, form = form, mode = 'prod', family = family)
    }
    models <- c(models_const, models_sum, models_prod)
  } else {
    models <- NULL
  }
  return(models)
}

find_closing_bracket <- function(i, form) {
  rest_form <- form[(i + 1):length(form)]
  extra_prods <- which(substr(rest_form, nchar(rest_form), nchar(rest_form)) == '(') + i
  brackets <- which(form[i + 1:length(form)] == ')') + i
  res_brack <- 1
  remove <- NULL
  for (k in i:length(form)) {
    if (k %in% brackets) {
      res_brack <- res_brack - 1
    } else if (k %in% extra_prods) {
      res_brack <- res_brack + 1
    }
    if (res_brack == 0) {
      remove <- c(remove, k)
    }
  }
  k <- remove[1]
  return(k)
}

remove_element <- function(i, x, form, mode, family) {
  if (mode == 'prod') {
    k <- find_closing_bracket(i, form)
    form <- form[-k]
  }
  if (mode == 'const') {
    form[i] <- '- 1'
  } else {
    form <- form[-i]
    while (substr(form[i - 1], nchar(form[i - 1]), nchar(form[i - 1])) == '(' && form[i] == ')') {
      form <- form[-i]
      form[i - 1] <- substr(form[i - 1], 1, nchar(form[i - 1]) - 4)
    }
  }
  res <- form_to_model(x, form, family)
  return(res)
}

increase_models <- function(x, form, variants, family) {
  sum_models <- lapply(variants, add_element, x = x, form = form, mode = 'sum', family = family)
  mult_models <- lapply(variants, add_element, x = x, form = form, mode = 'prod', family = family)
  models <- c(sum_models, mult_models)
  return(models)
}

init_model <- function(x, family) {
  form <- c('x ~', '- 1')
  res <- form_to_model(x, form, family)
  return(res)
}

form_to_model <- function(x, form, family) {
  newlm <- glm(paste(form, collapse = ' '), data = x, family = family)
  BIC <- BIC(newlm)
  return(list(BIC = BIC, form = form, lm = newlm))
}

add_element <- function(addition, x, form, mode, family) {
  init_form <- form
  run <- TRUE
  if (substr(addition, 1, 8) == 'Periodic') {
    i <- substr(addition, 10, nchar(addition))
  }
  if (substr(addition, 1, 7) == 'Periods') {
    i <- substr(addition, 9, nchar(addition))
  }
  if (addition == 'Constant' && mode != 'sum') {
    run <- FALSE
  } else if (addition == 'Constant' && mode == 'sum') {
    const_pos <- which(form == '- 1')
    if (length(const_pos) > 0) {
      form[[which(form == '- 1')]] <- '+ 1'
    } else {
      run <- FALSE
    }
  } else if (addition == 'Linear' && mode == 'prod') {
    form <- c(form[1], '+ linear : (', form[2:length(form)], ')')
    form <- create_quadratic(form, 2)
  } else if (addition == 'Linear' && mode == 'sum') {
    form <- append(form, '+ linear')
  } else if (substr(addition, 1, 8) == 'Periodic' && mode == 'prod') {
    form <- c(form[1], paste0('+ lag_', i, ' : ('), form[2:length(form)], ')')
    form <- create_quadratic(form, 2)
  } else if (substr(addition, 1, 8) == 'Periodic' && mode == 'sum') {
    form <- append(form, paste0('+ lag_', i))
  } else if (substr(addition, 1, 7) == 'Periods' && mode == 'prod') {
    form <- c(form[1], paste0('+ period_', i, ' : ('), form[2:length(form)], ')')
  } else if (substr(addition, 1, 7) == 'Periods' && mode == 'sum') {
    form <- append(form, paste0('+ period_', i))
  }
  if (run) {
    res <- form_to_model(x, form, family)
  } else{
    res <- form_to_model(x, init_form, family)
  }
  return(res)
}

create_quadratic <- function(form, i) {
  similar_prod <- which(form == form[2])
  similar_sum <- which(form == substr(form[2], 1, nchar(form[2]) - 4))
  similar_prod <- setdiff(similar_prod, 2)
  similar_sum <- setdiff(similar_sum, 2)
  if (length(similar_prod) > 0) {
    form[similar_prod] <- paste0('+ I(', substr(form[similar_prod], 3, nchar(form[similar_prod]) - 4), ') : (')
  }
  if (length(similar_sum) > 0) {
    form[similar_sum] <- paste0('+ I(', substr(form[similar_sum], 3, nchar(form[similar_sum])), ')')
  }
  return(form)
}

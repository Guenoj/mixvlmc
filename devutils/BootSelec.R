#' Predictor of the following state of a discrete time series for a vlmc
#'
#' This function predict the next state of a time series from the distribution estimated by the
#' given vlmc object.
#'
#'
#' @param object a fitted vlmc object.
#' @param ctx a time series adapted to the vlmc object.
#' @param ... additional arguments.
#'
#' @export
#' @examples
#' pc <- powerconsumption[powerconsumption$week == 5, ]
#' dts <- cut(pc$active_power, breaks = c(0, quantile(pc$active_power, probs = c(0.25, 0.5, 0.75, 1))))
#' model <- vlmc(dts, min_size = 5)
#' predict(model, dts[1:5])
predict.vlmc <- function(object, ts, only_last = F, ...) {
  max_depth <- depth(object)
  if (!is.null(ts)) {
    assertthat::assert_that((typeof(ts) == typeof(object$vals)) && (class(ts) == class(object$vals)),
      msg = "ts is not compatible with the model state space"
    )
    dts <- to_dts(ts, object$vals)
    ctx <- rev(dts$ix) + 1
  } else {
    # Prediction du 1ER ??
  }
  pred_vals <- object$vals ### ATTENTION
  if (only_last) {
    if (max_depth == 0) {
      return(pred_vals[which.max(match_context(object, ctx)$tree$f_by)]) # In simulate there is no $tree but here I got a numeric(0)
    } else {
      return(pred_vals[which.max(match_context(object, ctx)$tree$f_by)])
    }
  } else {
    if (max_depth == 0) {
      pred <- pred_vals[which.max(match_context(object, ctx)$f_by)]
      return(matrix(c(ts, rep(pred, length(ctx))), ncol = 2))
    } else {
      thing <- matrix(nrow = length(ctx), ncol = 2)
      for (ii in 1:length(ctx)) {
        ac_ctx <- ctx[(length(ctx) - ii + 1):length(ctx)]
        thing[ii, ] <- c(ts[ii], pred_vals[which.max(match_context(object, ac_ctx)$tree$f_by)])
      }
      return(thing)
    }
  }
}


#' Select the best model vlmc by a bootstrap method
#'
#' This function returns the best vlmc model fitting a time series by a bootstrap method.
#'
#'
#' @param dts a time series.
#' @param n the length of the bootstrap sequences.
#' @param B the number of bootstrapp.
#' @param ... additional arguments.
#'
#' @export
#' @examples
#' California_centre <- data.frame(longitude = -119.449444, latitude = 37.166111)
#' distances <- geodist(globalearthquake[, c("longitude", "latitude")], California_centre, measure = "geodesic")
#' California_earth_quakes <- globalearthquake[distances < 2e6, ] ## distances are in meters
#' California_weeks <- rep(0, max(globalearthquake$nbweeks))
#' California_weeks[California_earth_quakes$nbweeks] <- 1
#' selec_by_boot(California_weeks, 400, 60)$bestmodel
selec_by_boot <- function(dts, n, B) {
  thing <- list()
  obj <- vlmc(dts)
  alphas <- c(obj$alpha, cutoff(obj))
  thing$alphas <- alphas
  # print(alphas)
  score <- rep(0, length(alphas))
  # print(score)
  for (b in 1:B) {
    # print(b)
    new_dts <- simulate.vlmc(obj, nsim = length(dts) + n + 1, init = dts)[(length(dts) + 1):(length(dts) + n + 1)]
    # print(new_dts)
    for (ii in 1:length(alphas)) {
      v0 <- vlmc(new_dts[1:(length(new_dts) - 1)], alpha = alphas[ii])
      score[ii] <- score[ii] + (predict.vlmc2(v0, new_dts[1:(length(new_dts) - 1)], only_last = T) == new_dts[length(new_dts)]) * 1
    }
  }
  thing$score <- score
  thing$bestmodel <- prune(obj, alpha = alphas[which.max(score)])
  return(thing)
}

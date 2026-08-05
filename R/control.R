#' Control Numerical Optimization in `regcorr`
#'
#' Creates a list of numerical controls for the safeguarded Newton optimizer
#' used by \code{\link{regcorr}}. The defaults require both a small relative
#' coefficient
#' change and a small score norm before convergence is reported.
#'
#' @param maxit Maximum number of Newton iterations.
#' @param reltol Relative tolerance for the maximum absolute coefficient
#'   change.
#' @param gradtol Tolerance for the maximum absolute component of the score.
#' @param max_halving Maximum number of times a Newton step may be halved.
#' @param min_step Smallest step size considered during step-halving.
#' @param hessian_max_condition Largest permitted Hessian condition number.
#'   Larger values are treated as numerically ill-conditioned.
#' @param boundary_eps Numerical safety distance used when checking correlation
#'   and bivariate Bernoulli joint-probability boundaries.
#'
#' @return A list of class `"regcorr_control"`.
#'
#' @examples
#' ctrl <- regcorr_control(maxit = 150, gradtol = 1e-7)
#' ctrl
#'
#' @export
regcorr_control <- function(maxit = 100L,
                            reltol = 1e-8,
                            gradtol = 1e-6,
                            max_halving = 20L,
                            min_step = 1e-8,
                            hessian_max_condition = 1e12,
                            boundary_eps = 1e-10) {
  .check_control_integer(maxit, "maxit", lower = 1L)
  .check_control_scalar(reltol, "reltol", lower = 0, open_lower = TRUE)
  .check_control_scalar(gradtol, "gradtol", lower = 0, open_lower = TRUE)
  .check_control_integer(max_halving, "max_halving", lower = 0L)
  .check_control_scalar(min_step, "min_step", lower = 0, upper = 1,
                        open_lower = TRUE)
  .check_control_scalar(hessian_max_condition, "hessian_max_condition",
                        lower = 1, open_lower = TRUE)
  .check_control_scalar(boundary_eps, "boundary_eps", lower = 0, upper = 0.1,
                        open_lower = TRUE, open_upper = TRUE)

  structure(
    list(
      maxit = as.integer(maxit),
      reltol = as.numeric(reltol),
      gradtol = as.numeric(gradtol),
      max_halving = as.integer(max_halving),
      min_step = as.numeric(min_step),
      hessian_max_condition = as.numeric(hessian_max_condition),
      boundary_eps = as.numeric(boundary_eps)
    ),
    class = "regcorr_control"
  )
}

.check_control_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                  open_lower = FALSE, open_upper = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("`", name, "` must be one finite numeric value.", call. = FALSE)
  }
  below <- if (open_lower) x <= lower else x < lower
  above <- if (open_upper) x >= upper else x > upper
  if (below || above) {
    left <- if (open_lower) "(" else "["
    right <- if (open_upper) ")" else "]"
    stop("`", name, "` must lie in ", left, lower, ", ", upper, right, ".",
         call. = FALSE)
  }
  invisible(x)
}

.check_control_integer <- function(x, name, lower = 0L) {
  .check_control_scalar(x, name, lower = lower)
  if (x != floor(x)) {
    stop("`", name, "` must be a whole number.", call. = FALSE)
  }
  invisible(x)
}

.validate_regcorr_control <- function(control) {
  if (is.null(control)) {
    return(regcorr_control())
  }
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control)))) {
    stop("`control` must be a `regcorr_control()` object or a named list.",
         call. = FALSE)
  }

  defaults <- unclass(regcorr_control())
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown `control` component(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  for (name in names(control)) {
    defaults[[name]] <- control[[name]]
  }
  do.call(regcorr_control, defaults)
}

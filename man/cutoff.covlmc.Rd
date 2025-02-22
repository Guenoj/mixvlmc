% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covlmc.R
\name{cutoff.covlmc}
\alias{cutoff.covlmc}
\title{Cutoff values for pruning the context tree of a VLMC with covariates}
\usage{
\method{cutoff}{covlmc}(
  vlmc,
  mode = c("quantile", "native"),
  raw = FALSE,
  tolerance = .Machine$double.eps^0.5,
  ...
)
}
\arguments{
\item{vlmc}{a fitted covlmc model.}

\item{mode}{specify whether the results should be "native" likelihood ratio
values or expressed in a "quantile" scale of a chi-squared distribution.
For covlmc, only the quantile scale is supported.}

\item{raw}{specify whether the returned values should be limit values computed in the model or
modified values that guarantee pruning (see details)}

\item{tolerance}{specify the minimum separation between two consecutive values of
the cut off in native mode (before any transformation). See details.}

\item{...}{additional arguments for the cutoff function.}
}
\value{
a vector of cut off values, \code{NULL} is none can be computed
}
\description{
This function returns all the cutoff values that should induce a pruning of
the context tree of a VLMC with covariates.
}
\details{
Notice that the list of cutoff values returned by the function is not as
complete as the one computed for a VLMC without covariates. Indeed, pruning
the coVLMC tree creates new pruning opportunities that are not evaluated
during the construction of the initial model, while all pruning opportunities
are computed during the construction of a VLMC context tree. Nevertheless,
the largest value returned by the function is guaranteed to produce the least
pruned tree consistent with the reference one.

For large COVLMC, some cut off values can be almost identical, with a
difference of the order of the machine epsilon value. The \code{tolerance} parameter
is used to keep only values that are different enough. This is done in the
quantile scale, before transformations implemented when \code{raw} is \code{FALSE}.

Notice that the loglikelihood scale is not directly useful in coVLMC as
the differences in model sizes are not constant through the pruning process.
As a consequence, the "native" scale is not supported by this function.

Setting \code{raw} to \code{TRUE} removes the small perturbation that are subtracted from
the log-likelihood ratio values computed from the coVLMC (in quantile scale).
}
\examples{
pc <- powerconsumption[powerconsumption$week == 5, ]
dts <- cut(pc$active_power, breaks = c(0, quantile(pc$active_power, probs = c(0.5, 1))))
m_nocovariate <- vlmc(dts)
draw(m_nocovariate)
dts_cov <- data.frame(day_night = (pc$hour >= 7 & pc$hour <= 17))
m_cov <- covlmc(dts, dts_cov, min_size = 5)
draw(m_cov)
cutoff(m_cov)
}

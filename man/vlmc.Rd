% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vlmc.R
\name{vlmc}
\alias{vlmc}
\title{Fit a Variable Length Markov Chain (VLMC)}
\usage{
vlmc(
  x,
  alpha = 0.05,
  cutoff = NULL,
  min_size = 2L,
  max_depth = 100L,
  prune = TRUE,
  keep_match = FALSE
)
}
\arguments{
\item{x}{a discrete time series; can be numeric, character, factor or logical.}

\item{alpha}{number in (0,1] (default: 0.05) cut off value in quantile scale
in the pruning phase.}

\item{cutoff}{non negative number: cut off value in native (likelihood ratio)
scale in the pruning phase. Defaults to the value obtained from \code{alpha}.
Takes precedence over \code{alpha} is specified.}

\item{min_size}{integer >= 1 (default: 2). Minimum number of observations for
a context in the growing phase of the context tree.}

\item{max_depth}{integer >= 1 (default: 100). Longest context considered in
growing phase of the context tree.}

\item{prune}{logical: specify whether the context tree should be pruned
(default behaviour).}

\item{keep_match}{logical: specify whether to keep the context matches (default to FALSE)}
}
\value{
a fitted vlmc model.
}
\description{
This function fits a  Variable Length Markov Chain (VLMC) to a discrete time
series.
}
\details{
The VLMC is built using Bühlmann and Wyner's algorithm which consists in
fitting a context tree (see \code{\link[=ctx_tree]{ctx_tree()}}) to a time series and then pruning
it in such as way that the conditional distribution of the next state of the
time series given the context is significantly different from the
distribution given a truncated version of the context.

The construction of the context tree is controlled by \code{min_size} and
\code{max_depth}, exactly as in \code{\link[=ctx_tree]{ctx_tree()}}. Significativity is measured using a
likelihood ratio test (threshold can be specified in terms of the ratio
itself with \code{cutoff}) or in quantile scale with \code{alpha}.

Pruning can be postponed by setting \code{prune=FALSE}. Using a combination of
\code{\link[=cutoff]{cutoff()}} and \code{\link[=prune]{prune()}}, the complexity of the VLMC can then be adjusted.
Any VLMC model can be pruned after construction, \code{prune=FALSE} is a
convenience parameter to avoid setting \code{alpha=1} (which essentially prevents
any pruning). Automated model selection is provided by \code{\link[=tune_vlmc]{tune_vlmc()}}.
}
\examples{
pc <- powerconsumption[powerconsumption$week == 5, ]
dts <- cut(pc$active_power, breaks = c(0, quantile(pc$active_power, probs = c(0.25, 0.5, 0.75, 1))))
model <- vlmc(dts)
draw(model)
depth(model)
## reduce the detph of the model
shallow_model <- vlmc(dts, max_depth = 3)
draw(shallow_model, prob = FALSE)
## improve probability estimates
robust_model <- vlmc(dts, min_size = 25)
draw(robust_model, prob = FALSE) ## show the frequencies
draw(robust_model)
}
\references{
Bühlmann, P. and Wyner, A. J. (1999), "Variable length Markov
chains. Ann. Statist." 27 (2)
480-513 \doi{10.1214/aos/1018031204}
}
\seealso{
\code{\link[=cutoff]{cutoff()}}, \code{\link[=prune]{prune()}} and \code{\link[=tune_vlmc]{tune_vlmc()}}
}

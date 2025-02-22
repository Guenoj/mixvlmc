% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covlmc_tune.R
\name{tune_covlmc}
\alias{tune_covlmc}
\title{Fit an optimal Variable Length Markov Chain with Covariates (coVLMC)}
\usage{
tune_covlmc(
  x,
  covariate,
  criterion = c("BIC", "AIC"),
  initial = c("truncated", "specific", "extended"),
  min_size = 5,
  max_depth = 100,
  verbose = 0,
  save = c("best", "initial", "all"),
  trimming = c("full", "partial", "none"),
  best_trimming = c("none", "partial", "full")
)
}
\arguments{
\item{x}{a discrete time series; can be numeric, character, factor and logical.}

\item{covariate}{a data frame of covariates.}

\item{criterion}{criterion used to select the best model. Either \code{"BIC"} (default)
or \code{"AIC"} (see details).}

\item{initial}{specifies the likelihood function, more precisely the way the
first few observations for which contexts cannot be calculated are integrated
in the likelihood. See \code{\link[=loglikelihood]{loglikelihood()}} for details.}

\item{min_size}{integer >= 1 (default: 2). Minimum number of observations for
a context in the growing phase of the initial context tree.}

\item{max_depth}{integer >= 1 (default: 100). Longest context considered in
growing phase of the initial context tree (see details).}

\item{verbose}{integer >= 0 (default: 0). Verbosity level of the pruning process.}

\item{save}{specify which BIC models are saved during the pruning process. The default
value \code{"best"} asks the function to keep only the best model according to
the \code{criterion}. When \code{save="initial"} the function keeps \emph{in addition} the
initial (complex) model which is then pruned during the selection process.
When \code{save="all"}, the function returns all the models considered during the
selection process. See details for memory occupation.}

\item{trimming}{specify the type of trimming used when saving the intermediate models,
see details.}

\item{best_trimming}{specify the type of trimming used when saving the best model
and the initial one (see details).}
}
\value{
a list with the following components:
\itemize{
\item \code{best_model}: the optimal VLMC
\item \code{criterion}: the criterion used to select the optimal VLMC
\item \code{initial}: the likelihood function used to select the optimal VLMC
\item \code{results}: a data frame with details about the pruning process
\item \code{saved_models}: a list of intermediate coVLMCs if \code{save="initial"} or
\code{save="all"}. It contains an \code{initial} component with the large coVLMC obtained
first and an \code{all} component with a list of all the \emph{other} coVLMC obtained
by pruning the initial one.
}
}
\description{
This function fits a Variable Length Markov Chain with Covariates (coVLMC)
to a discrete time series coupled with a time series of covariates
by optimizing an information criterion (BIC or AIC).
}
\details{
This function automates the process of fitting a large coVLMC to a discrete time
series with \code{\link[=covlmc]{covlmc()}} and of pruning the tree (with \code{\link[=cutoff]{cutoff()}} and \code{\link[=prune]{prune()}})
to get an optimal with respect to an information criterion. To avoid missing
long term dependencies, the function uses the \code{max_depth} parameter as an initial
guess but then relies on an automatic increase of the value to make sure the
initial context tree is only limited by the \code{min_size} parameter. The initial
value of the \code{alpha} parameter of \code{\link[=covlmc]{covlmc()}} is also set to a conservative value
to avoid prior simplification of the context tree.

Once the initial coVLMC is obtained, the \code{\link[=cutoff]{cutoff()}} and \code{\link[=prune]{prune()}} functions are
used to build all the coVLMC models that could be generated using smaller values of
the alpha parameter. The best model is selected from this collection,
including the initial complex tree, as the one that minimizes the chosen
information criterion.
}
\section{Memory occupation}{


\code{covlmc} objects tend to be large and saving all the models during the search for
the optimal model can lead to an unreasonable use of memory. To avoid this problem,
models are kept in trimmed form only using \code{\link[=trim.covlmc]{trim.covlmc()}} with \code{keep_model=FALSE}.
Both the initial model and the best one are saved untrimmed. This default
behaviour corresponds to \code{trimming="full"}. Setting \code{trimming="partial"} asks the function
to use \code{keep_model=TRUE} in \code{\link[=trim.covlmc]{trim.covlmc()}} for intermediate models. Finally,
\code{trimming="none"} turns off trimming, which is discouraged expected for small data sets.

In parallel processing contexts (e.g. using \link[foreach:foreach]{foreach::\%dopar\%}), the memory
occupation of the results can become very large as models tend to keep
environments attached to the formulas. In this situation, it is highly recommended
to trim all saved models, including the best one and the initial one. This can
be done via the \code{best_trimming} parameter whose possible values are identical
to the ones of \code{trimming}.
}

\examples{
pc <- powerconsumption[powerconsumption$week \%in\% 5:7, ]
dts <- cut(pc$active_power, breaks = c(0, quantile(pc$active_power, probs = c(0.5, 1))))
dts_cov <- data.frame(day_night = (pc$hour >= 7 & pc$hour <= 17))
dts_best_model_tune <- tune_covlmc(dts, dts_cov)
draw(as_covlmc(dts_best_model_tune))
}
\seealso{
\code{\link[=covlmc]{covlmc()}}, \code{\link[=cutoff]{cutoff()}} and \code{\link[=prune]{prune()}}
}

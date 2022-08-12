
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Variable Length Markov Chains

<!-- badges: start -->

[![R-CMD-check](https://github.com/fabrice-rossi/mixvlmc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fabrice-rossi/mixvlmc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`mixvlmc` implements variable length Markov chains (VLMC) and variable
length Markov chains with covariates, as described in:

-   [Bühlmann, P. and Wyner, A. J. (1999), Variable length Markov
    chains. Ann. Statist. 27 (2)
    480-513](https://dx.doi.org/10.1214/aos/1018031204)
-   [Zanin Zambom, A., Kim, S. and Lopes Garcia, N. (2022), Variable
    length Markov chain with exogenous covariates. J. Time Ser. Anal.,
    43 (2) 312-328](https://doi.org/10.1111/jtsa.12615)

## Installation

The package is only available from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("fabrice-rossi/mixvlmc")
```

## Usage

### Variable length Marvkoc chains

Variable length Marvkov chains (VLMC) are sparse high order Markov
chains. They can be used to model time series with discrete values
(states) with a mix of small order dependencies for certain states and
higher order dependencies for other states. For instance, with a binary
time series, the probability of observing 1 a time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
could be constant whatever the past states if the last one (at time
![t-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t-1 "t-1"))
was 1, but could depend on states at time
![t-3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t-3 "t-3")
and
![t-2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t-2 "t-2")
if the state was 0 at time
![t-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t-1 "t-1").
A collection of past states that determines completely the transition
probabilities is a *context* of the VLMC. Read
`vignette("context-trees")` for details about contexts and context tree.

VLMC with covariates are extension of VLMC in which transition
probabilities (probabilities of the next state given the past) can be
influenced by the past values of some covariates (in addition to the
past values of the time series itself). Each context is associated to a
logistic model that maps the (past values of the) covariates to
transition probabilities.

### Fitting a VLMC

The package is loaded in a standard way.

``` r
library(mixvlmc)
## we use ggplot2 for graphics but this is not a dependency of the package
library(ggplot2)
```

The main function of VLMC is `vlmc()` which can be called on a time
series represented by a numerical vector or a factor.

``` r
set.seed(0)
x <- sample(c(0L, 1L), 200, replace = TRUE)
model <- vlmc(x)
model
#> VLMC context tree on 0, 1 
#>  cutoff: 1.921 (quantile: 0.05)
#>  Maximum context length: 6 
#>  Number of contexts: 11
```

The default parameters of `vlmc()` will tend to produce overly complex
VLMC in order to avoid missing potential structure in the time series.
In the example above, we expect the optimal VLMC to be a constant
distribution as the sample is independent and uniformly distributed (is
has not temporal structure). The default paramaters give here an overly
complex model, as illustrated by its text based representation

``` r
draw(model)
#> * (0.505, 0.495)
#> '-- 1 (0.4848, 0.5152)
#>     +-- 0 (0.5319, 0.4681)
#>     |   '-- 1 (0.5, 0.5)
#>     |       '-- 0 (0.4444, 0.5556)
#>     |           '-- 0 (0.4286, 0.5714)
#>     |               +-- 0 (1, 0)
#>     |               '-- 1 (0, 1)
#>     '-- 1 (0.4314, 0.5686)
#>         '-- 0 (0.2727, 0.7273)
#>             '-- 0 (0.3846, 0.6154)
#>                 '-- 0 (0.5, 0.5)
#>                     +-- 0 (0, 1)
#>                     '-- 1 (1, 0)
```

The representation uses crude ascii art to display the contexts of the
VLMC organized in to tree:

-   the root `*` corresponds to an empty context;
-   one can read contexts by following branches (represented by 2 dashes
    `--`) down to its end (the leaf): for instance
    ![(1, 0, 1, 0, 0, 0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2C%200%2C%201%2C%200%2C%200%2C%200%29 "(1, 0, 1, 0, 0, 0)")
    is one of the context of the tree.

Here the context
![(1, 0, 1, 0, 0, 0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2C%200%2C%201%2C%200%2C%200%2C%200%29 "(1, 0, 1, 0, 0, 0)")
is associated to the transition probabilities
![(1, 0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2C%200%29 "(1, 0)").
This means that when one observe this context in the time series, it is
*always* followed by a 0. Notice that contexts are traditionally written
from the most recent value to the oldest one. Thus, the context
![(1, 0, 1, 0, 0, 0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281%2C%200%2C%201%2C%200%2C%200%2C%200%29 "(1, 0, 1, 0, 0, 0)")
corresponds to the sub time series
![(0, 0, 0, 1, 0, 1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280%2C%200%2C%200%2C%201%2C%200%2C%201%29 "(0, 0, 0, 1, 0, 1)").

#### BIC based model selection

The VLMC above is obviously overfitting to the time series, as
illustrated by the 0/1 transition probabilities. A classical way to
select a good model is to maximize the
[BIC](https://en.wikipedia.org/wiki/Bayesian_information_criterion). In
`mixvlmc` this can be done easily by *pruning* a complex VLMC using a
combination of `cutoff` and `prune`, as follows:

``` r
alphas <- cutoff(model)
pruned_models <- c(list(model), lapply(alphas, function(a) prune(model, a)))
bics <- sapply(pruned_models, BIC)
best_model <- pruned_models[[which.min(bics)]]
draw(best_model)
#> * (0.505, 0.495)
```

As expected, we end up with a constant model.

In time series with actual temporal patterns, the optimal model will be
more complex. As a very basic illustrative example, let us consider the
`sunspot.year` time series and turn it into a binary one, with high
activity associated to a number of sun spots larger than the median
number.

``` r
sun_activity <- as.factor(ifelse(sunspot.year >= median(sunspot.year), "high", "low"))
```

We adjust a VLMC with a very conservative value of `alpha`:

``` r
sun_model <- vlmc(sun_activity, alpha=0.1)
sun_model
#> VLMC context tree on high, low 
#>  cutoff: 1.353 (quantile: 0.1)
#>  Maximum context length: 14 
#>  Number of contexts: 27
```

Then we proceed as above:

``` r
sun_alphas <- cutoff(sun_model)
pruned_sun_models <- c(list(sun_model), lapply(sun_alphas, function(a) prune(sun_model, a)))
sun_bics <- sapply(pruned_sun_models, BIC)
sun_bic_fig <- data.frame(alpha=c(0.1, sun_alphas), BIC = sun_bics)
ggplot(sun_bic_fig, aes(x=alpha, y=BIC)) + 
  geom_line() 
```

<img src="man/figures/README-sunspots_bic-1.png" width="100%" />

The plot shows that simpler models are too simple as the BIC increases
when pruning becomes strong enough. The best model remains rather
complex (as expected based on the periodicity of the Solar cycle):

``` r
best_sun_model <- pruned_sun_models[[which.min(sun_bics)]]
draw(best_sun_model)
#> * (0.5052, 0.4948)
#> +-- high (0.8207, 0.1793)
#> |   +-- high (0.7899, 0.2101)
#> |   |   +-- high (0.7447, 0.2553)
#> |   |   |   +-- high (0.6571, 0.3429)
#> |   |   |   |   '-- low (0.9167, 0.08333)
#> |   |   |   '-- low (1, 0)
#> |   |   '-- low (0.96, 0.04)
#> |   '-- low (0.9615, 0.03846)
#> '-- low (0.1888, 0.8112)
#>     +-- high (0, 1)
#>     '-- low (0.2328, 0.7672)
#>         +-- high (0, 1)
#>         '-- low (0.3034, 0.6966)
#>             '-- high (0.07692, 0.9231)
```

### Fitting a VLMC with covariates

To illustrate the use of covariates, we use the power consumption data
set included in the package. We consider a week of electricity usage as
follows:

``` r
pc_week_5 <- powerconsumption[powerconsumption$week == 5, ]
elec <- pc_week_5$active_power
ggplot(pc_week_5, aes(x=date_time, y=active_power)) + 
  geom_line() + 
  xlab("t") +
  ylab("Activer power (kW)")
```

<img src="man/figures/README-active_power_week_5-1.png" width="100%" />
The time series displays some typical patterns of electricity usage:

-   low active power at night (typically below 0.4 kW);
-   standard use between 0.4 and 2 kW;
-   peak use above 2 kW.

We build a discrete time series from those (somewhat arbitrary)
thresholds:

``` r
elec_dts <- cut(elec, breaks = c(0, 0.4, 2, 8), labels=c("low", "typical", "high"))
```

The best VLMC model is quite simple. It is almost a standard order one
Markov chain, up to the order 2 context used when the active power is
typical.

``` r
elec_vlmc <- vlmc(elec_dts)
e_alphas <- cutoff(elec_vlmc)
pruned_elec_vlmcs <- c(list(elec_vlmc), lapply(e_alphas, function(a) prune(elec_vlmc, a)))
e_bics <- sapply(pruned_elec_vlmcs, BIC)
best_elec_vlmc <- pruned_elec_vlmcs[[which.min(e_bics)]]
draw(best_elec_vlmc)
#> * (0.1667, 0.5496, 0.2837)
#> +-- low (0.7665, 0.2335, 0)
#> '-- typical (0.0704, 0.8466, 0.08303)
#> |   '-- low (0.3846, 0.5385, 0.07692)
#> '-- high (0.003497, 0.1573, 0.8392)
```

As pointed about above, low active power tend to correspond to night
phase. We can include this information by introducing a day covariate as
follows:

``` r
elec_cov <- data.frame(day = (pc_week_5$hour >= 7 & pc_week_5$hour <= 17))
```

A VLMC with covariate is estimated using the `covlmc` function:

``` r
elec_covlmc <- covlmc(elec_dts, elec_cov, min_size = 5, alpha=0.5)
draw(elec_covlmc, time_sep= "|")
#> *
#> +-- low (0.03031 [ -1.558|1.006 ])
#> '-- typical
#> |   +-- low (0.2181 [ 0.3567|-27.81|27.81
#> |   |                 -1.253|-14.39|13.69 ])
#> |   '-- typical (collapsing: 0.3619)
#> |   |   +-- low (0.07525 [ 1.099 |18.33
#> |   |   |                  -18.22|35.57 ])
#> |   |   '-- typical (0.2843 [ 2.791|0.4573
#> |   |   |                     0.47 |-0.06454 ])
#> |   |   '-- high (0.4491 [ 18.89|-16.25
#> |   |                      16.25|-15.56 ])
#> |   '-- high (0.1715 [ 2.015 |16.18
#> |                      0.6931|16.61 ])
#> '-- high (0.05754 [ 17.41|-14.23
#>                     19.38|-14.88 ])
```

Pruning a covlmc model is slightly more complicated than in the case of
vlmc models. The main difficulty is that cut off values cannot be
computed directly from a complex model as new values will appear during
the pruning process. A typical way to circumvent this difficulty is to
recompute cut off values after each pruning operation, as follows:

``` r
elec_covlmc_models <- list(elec_covlmc)
current_model <- elec_covlmc
alpha <- cutoff(current_model)[1]
while(TRUE) {
  new_model <- prune(current_model, alpha)
  elec_covlmc_models <- c(elec_covlmc_models, list(new_model))
  current_model <- new_model
  new_alphas <- cutoff(current_model)
  if(is.null(new_alphas) || all(new_alphas>=alpha)) {
    break
  } else {
    alpha <- max(new_alphas[new_alphas<alpha])
  }
}
elec_covlmc_bics <- sapply(elec_covlmc_models, BIC)
best_elec_covlmc <- elec_covlmc_models[[which.min(elec_covlmc_bics)]]
draw(best_elec_covlmc, model = "full", time_sep = " | ")
#> *
#> +-- low (0.03031 [ (I)    | day_1TRUE
#> |                  -1.558 | 1.006 ])
#> '-- typical (collapsing: 3.663e-09)
#> |   +-- low (0.8071 [ (I)
#> |   |                 0.3365
#> |   |                 -1.609 ])
#> |   '-- typical (0.3168 [ (I)
#> |   |                     2.937
#> |   |                     0.3747 ])
#> |   '-- high (0.1715 [ (I)
#> |                      2.773
#> |                      1.705 ])
#> '-- high (0.05754 [ (I)
#>                     3.807
#>                     5.481 ])
```

As in the VLMC case, the optimal model is simple. The two contexts
`typical` and `high` do not use the covariate. The `low` context does
not switch to a `high` context (hence the single row of parameters) but
uses the covariate. As expected, the probability of switching from `low`
to `typical` is larger during the day.

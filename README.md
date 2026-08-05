# regcorr: Regression Models of Pearson Correlation Coefficient

The `regcorr` R package provides statistical tools to evaluate how covariates of interest influence the strength of the Pearson correlation coefficient between two responses. It supports both **continuous (bivariate normal)** and **bivariate binary (Bernoulli)** responses without requiring repeated measures.

This package implements the likelihood-based estimators of Dufera, Liu, and Xu (2023), with bootstrap-based inference and explicit numerical safeguards. The statistical likelihoods are unchanged from the paper.

## 👥 Authors
* **Ze Lin**
* **Bo Li**
* **Jinyao Shen**
* **Jiangyan Zhao**
* **Jin Xu**

## ✨ Key Features
* **Bivariate Normal Responses:** Models correlation using either the logistic link for positive correlations or the hyperbolic tangent link for correlations in `(-1, 1)`.
* **Bivariate Binary Responses:** Uses the same two correlation links while retaining the Bernoulli likelihood and fitted marginal logistic models of Dufera et al. (2023).
* **Safeguarded Newton Optimization:** Uses deterministic step-halving, finite-domain checks, Hessian condition-number checks, and convergence criteria based jointly on relative coefficient change and the final score norm. Random jitter restarts are not part of the routine optimization path.
* **Explicit Bernoulli Feasibility:** Every candidate update must keep all four fitted joint probabilities strictly inside their numerical domain.
* **Bootstrap Robustness:** Bootstrap inference remains the supported inference method. A failed or non-converged replicate is excluded rather than terminating the fit, and the returned object records valid, non-converged, and errored replicate counts.

## 📥 Installation

`regcorr` is available on CRAN, and the development version can be installed from GitHub:

```r
# CRAN release
install.packages("regcorr")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("lonze-nb/regcorr")
```

## 🚀 Quick Start (Usage)
Here is a basic example of fitting the correlation regression model with the unified `regcorr()` function:
```r
library(regcorr)

# 1. Generate simulated bivariate normal data (n = 500)
set.seed(123)
n <- 500
x <- runif(n)
rho <- plogis(0.25 + 1 * x)              # logistic link for the correlation
z1 <- rnorm(n)
y1 <- z1
y2 <- z1 * rho + rnorm(n) * sqrt(1 - rho^2)
dat <- data.frame(y1 = y1, y2 = y2, x = x)

# 2. Fit the model (non-0/1 response is detected as "normal" automatically)
fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 100)

# 3. View results
print(fit)
summary(fit)
```

Numerical tolerances can be changed without altering existing calls:

```r
fit <- regcorr(
  cbind(y1, y2) ~ x,
  data = dat,
  nboot = 100,
  control = regcorr_control(maxit = 150, gradtol = 1e-7)
)
```

For bivariate binary responses (both columns 0/1), `regcorr()` selects the
binary model automatically; set `type = "binary"` to force it.

## 📖 References
Dufera, A. G., Liu, T., & Xu, J. (2023). Regression models of Pearson correlation coefficient. *Statistical Theory and Related Fields*, 7(2), 97-106. https://doi.org/10.1080/24754269.2023.2164970

# regcorr: Regression Models of Pearson Correlation Coefficient

The `regcorr` R package provides statistical tools to evaluate how covariates of interest influence the strength of the Pearson correlation coefficient between two responses. It supports both **continuous (bivariate normal)** and **bivariate binary (Bernoulli)** responses without requiring repeated measures.

This package replicates and robustly extends the methodologies for likelihood-based inference using Newton-Raphson estimation and bootstrap-based significance testing. 

## 👥 Authors
* **Ze Lin**
* **Bo Li**
* **Jinyao Shen**

## ✨ Key Features
* **Bivariate Normal Responses:** Models the Fisher z-transformed correlation (hyperbolic tangent link) against linear combinations of covariates.
* **Bivariate Binary Responses:** Models the correlation of dichotomous outcomes using a logistic link function.
* **Algorithmic Robustness:** Implements short-circuit evaluation and `try-catch` mechanisms to handle the notorious "perfect separation" and numerical instabilities (e.g., Hessian matrix explosion, `NA/NaN` generation) inherent in small-sample bivariate Bernoulli distributions.

## 📥 Installation

You can install the development version of `regcorr` directly from GitHub:

```{r}
# Install devtools if you haven't already
# install.packages("devtools")

devtools::install_github("lonze-nb/regcorr")
```

## 🚀 Quick Start (Usage)
Here is a basic example of fitting the correlation regression model with the unified `regcorr()` function:
```{r}
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

For bivariate binary responses (both columns 0/1), `regcorr()` selects the
binary model automatically; set `type = "binary"` to force it.

## 📖 References
This package is built based on the statistical framework proposed in related literature regarding regression models of Pearson correlation coefficients.

Dufera, A. G., Liu, T., & Xu, J. (2023). Regression models of Pearson correlation coefficient. Statistical Theory and Related Fields, 7(2), 97-106.
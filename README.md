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
Here is a basic example of generating bivariate normal data and fitting the correlation regression model:
```{r}
library(regcorr)

# 1. Set true parameters for simulation
true_beta <- c(0.25, 1, 0) # Intercept and two covariates
true_eta <- c(0, 0, 0)

# 2. Generate simulated Bivariate Normal data (n = 500)
set.seed(123)
my_data <- genDataBN(numSample = 500, p = 2, 
                     betaTrue = true_beta, 
                     eta1True = true_eta, 
                     eta2True = true_eta, 
                     link = "1") # 1: logistic link, 2: tanh link

# 3. Fit the model using Newton-Raphson iteration
fit_result <- NRfitBivNormal(Y = my_data$Y, X = my_data$X, 
                             betaIni = c(0.25, 0, 0), 
                             link = "1")

# 4. View the estimated parameters and number of iterations
print(fit_result$betaCurrent)
print(fit_result$numIter)
```

## 📖 References
This package is built based on the statistical framework proposed in related literature regarding regression models of Pearson correlation coefficients.

Dufera, A. G., Liu, T., & Xu, J. (2023). Regression models of Pearson correlation coefficient. Statistical Theory and Related Fields, 7(2), 97-106.

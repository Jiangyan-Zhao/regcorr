# Reproduce the simulation studies in Dufera, Liu & Xu (2023),
# "Regression models of Pearson correlation coefficient",
# Statistical Theory and Related Fields, 7(2), 97-106, Section 4.
#
# Generates Table 3 (type-I error, null model) and Table 4 (ERMSE, CR and
# power, nonnull model) under the design of Section 4.1:
# model x link x p x n, with numSimu simulation replications and numBoot
# bootstrap replications per fit.
library(regcorr)

# ---- settings (paper Section 4.1) ----
numSimu  <- 1000              # number of simulation replications
numBoot  <- 1000              # number of bootstrap replications (paper's B)
cSample  <- c(100, 200, 400, 1000)
nModel   <- 2                 # 1: bivariate normal; 2: bivariate Bernoulli
nLink    <- 2                 # 1: logistic; 2: tanh
nCova    <- 2:3               # number of covariates p

# nonnull model: only the first covariate is significant (paper Table 2)
beta_nonnull <- function(iLink, iCova)
  switch(iLink,
         "1" = c(0.25, 1,  rep(0, iCova - 1)),
         "2" = c(0.25, -1, rep(0, iCova - 1)))
# null model: no covariate effect, for type-I error (paper Table 3)
beta_null <- function(iLink, iCova) c(0.25, rep(0, iCova))

# run one beta configuration over model x link x p x n
simulate_table <- function(beta_fun) {
  tableTest <- NULL
  for (iModel in 1:nModel) {
    for (iLink in 1:nLink) {
      res_model_link <- NULL
      for (iCova in nCova) {
        # zero marginal mean models (paper sets gamma1 = gamma2 = 0)
        eta1True <- eta2True <- rep(0, iCova + 1)
        betaTrue <- beta_fun(iLink, iCova)
        betaIni  <- c(0.25, rep(0, iCova))

        res_model_link_p <- NULL
        for (iNumSample in cSample) {
          cat("model", iModel, "| link", iLink, "| p", iCova,
              "| n", iNumSample, "\n")
          res <- regcorr:::subRoutineTest(
            iNumSample, iCova, iLink, iModel,
            betaTrue, betaIni, eta1True, eta2True,
            numSimu, numBoot
          )
          # power reported in %, matching the paper's tables
          res_model_link_p <- rbind(
            res_model_link_p,
            c(res$RMSE, res$ConsistRate * 100, res$power[-2] * 100)
          )
        }
        res_model_link <- cbind(res_model_link, res_model_link_p)
      }
      tableTest <- rbind(tableTest, cbind(iModel, iLink, cSample, res_model_link))
    }
  }
  tableTest
}

# column names distinguish the two covariate numbers, e.g. p2.RMSE, p3.b1
make_cols <- function(p) {
  c(paste0("p", p, ".RMSE"), paste0("p", p, ".CR"),
    paste0("p", p, ".Global"), paste0("p", p, ".b", seq_len(p)))
}
write_table <- function(tab, file) {
  tab <- as.data.frame(tab)
  names(tab) <- c("model", "link", "n", unlist(lapply(nCova, make_cols)))
  print(tab)
  write.csv(tab, file = file, row.names = FALSE)
}

out_dir <- "inst/ext_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("===== Table 4: nonnull model (ERMSE, CR, power) =====\n")
write_table(simulate_table(beta_nonnull), file.path(out_dir, "tableTestPower.csv"))

cat("\n===== Table 3: null model (type-I error, %) =====\n")
write_table(simulate_table(beta_null), file.path(out_dir, "tableTestTypeI.csv"))

test_that("paper reproduction helper remains compatible with optimizer controls", {
  for (model in c("1", "2")) {
    set.seed(if (model == "1") 1601 else 1602)
    result <- subRoutineTest(
      numSample = 100, p = 1, link = "1", model = model,
      betaTrue = c(-0.5, 0.25), betaIni = c(0.1, 0.1),
      eta1True = c(0, 0), eta2True = c(0, 0),
      numSimu = 1, numBoot = 6
    )

    expect_equal(names(result), c("RMSE", "ConsistRate", "power"))
    expect_true(is.finite(result$RMSE))
    expect_true(is.finite(result$ConsistRate))
    expect_length(result$power, 3)
  }
})

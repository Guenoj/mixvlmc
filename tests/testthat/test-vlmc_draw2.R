test_that("draw of non-zero depth model", {
  withr::local_seed(0)
  data_set <- sample(c("A","B","C"), 500, replace = TRUE)
  d_vlmc <- vlmc(data_set)
  expect_snapshot_output(drawV(d_vlmc))
  expect_snapshot_output(drawV(d_vlmc,prob = F))
  expect_snapshot_output(drawV(d_vlmc,prob = NULL))
})

test_that("draw of zero depth model", {
  withr::local_seed(0)
  data_set <- sample(1:5, 50, replace = TRUE)
  d_vlmc <- vlmc(data_set)
  expect_snapshot_output(drawV(d_vlmc))
  expect_snapshot_output(drawV(d_vlmc,prob = F))
  expect_snapshot_output(drawV(d_vlmc,prob = NULL))
})

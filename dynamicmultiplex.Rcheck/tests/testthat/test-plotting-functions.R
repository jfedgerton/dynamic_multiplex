test_that("plotting functions are available", {
  expect_true(is.function(plot_multilayer_series))
  expect_true(is.function(animate_multilayer_gif))
  expect_true(is.function(plot_multilayer_alluvial))
})

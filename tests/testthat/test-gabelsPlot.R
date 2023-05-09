test_that("Gabel function moving average plot unit test", {
  b <- runif(1000, min=-2, max=2)
  c <- sample(2000:1000000, 1000, replace=TRUE)
  df <- data.frame(logFC.crude = b, gene.length = c)
  expected <- gabelsPlot(mat = df)
  #output list
  expect_type(expected, "list")
  #input list
  expect_type(df, "list")
  # logFC real numbers
  expect_true(typeof(b) == "double")
  # gene.length positive integers
  expect_true(typeof(c) == "integer")
})

test_that("Overlay Gabel function moving average plot unit test", {
  a <- runif(1000, min=-2, max=2)
  b <- runif(1000, min=-2, max=2)
  c <- sample(2000:1000000, 1000, replace=TRUE)
  df <- data.frame(comp.mat = a, logFC.crude = b, gene.length = c)
  expected <- overlayGabelsPlot(mat = df,comp.between1 = "(WT/WT)",
    comp.between2 = "(KO/WT)",bin.size = 200,
    shift.size = 40, confidenceinterval=0.50)
  #output list
  expect_type(expected, "list")
  #input list
  expect_type(df, "list")
  # comparison mat real numbers
  expect_true(typeof(a) == "double")
  # logFC real numbers
  expect_true(typeof(b) == "double")
  # gene.length positive integers
  expect_true(typeof(c) == "integer")
})


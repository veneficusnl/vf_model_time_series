context('get_models')

data('EIA')
train <- 1:(nrow(data) - 3)
labels <- data$Date
x <- data$Value

test_that("get_model works", {
  model <- get_model(x = x,method = 'time_series', lag = 2, labels = labels, family = gaussian('log') ,train = train )
  expect_that(get_model(x = x,method = 'time_series', lag = 2, labels = labels, family = gaussian('log') ,train = train ), takes_less_than(30))
  expect_is(model, 'model')
  expect_named(model, c("BIC", "form", "lm", "method", "periods", "x", "lag", "train",
                        "call"))
})

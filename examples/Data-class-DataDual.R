my_data <- DataDual(
  w = rnorm(8),
  x = c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
  y = c(0, 0, 0, 0, 0, 0, 1, 0),
  doseGrid = c(
    0.1, 0.5, 1.5, 3, 6,
    seq(from = 10, to = 80, by = 2)
  )
)
my_data

my_data <- DataDA(
  u = c(42, 30, 15, 5, 20, 25, 30, 60),
  t0 = c(0, 15, 30, 40, 55, 70, 75, 85),
  Tmax = 60,
  x = c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
  y = c(0, 0, 1, 1, 0, 0, 1, 0),
  doseGrid = c(0.1, 0.5, 1.5, 3, 6, seq(from = 10, to = 80, by = 2))
)

# Set up an empty data set.
empty_data <- DataDA(
  doseGrid = c(0.1, 0.5, 1, 1.5, 3, 6, seq(from = 10, to = 80, by = 2)),
  Tmax = 60
)
empty_data

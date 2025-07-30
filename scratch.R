dust2::dust_package(".")
devtools::load_all()

o <- run_model(
  event_time_on = 25
)
o$events

o$data$ipr
o$data$ipr_est

plot(o$data$ipr, type = "l")
lines(o$data$ipr_est, col = "red", type = "b")
abline(h = 0.1, col = "blue")

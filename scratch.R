dust2::dust_package(".")
devtools::load_all()

o <- run_model(
  event_time_on = 100,
  event_time_off = 100,
  time_end = 70
)
o$events

o$data$ipr
o$data$ipr_est

plot(o$data$ipr, type = "s")
lines(o$data$ipr_est, col = "red", type = "l")
abline(h = 0.1, col = "blue")

dust2::dust_package(".")
devtools::load_all()

o <- run_model(
  event_time_on = 10,
  event_time_off = 100,
  time_end = 150
)
o$events

o$data$rt_cumul
o$data$rt_est
o$data$flag

plot(o$data$rt_est, type = "s")
lines(o$data$rt_saved, col = "red", type = "l")
abline(h = 0.1, col = "blue")

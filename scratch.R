dust2::dust_package(".")
devtools::load_all()

# starting early: no state threshold crossed
o = run_model(
  event_time_on = 25
)
o
plot(o$data$I, type = "l")

# starting a bit later: weird, roots are completely off
o = run_model(
  event_time_on = 30
)
o
plot(o$data$I, type = "l")
abline(h = 100, col = 2)

t1_on = o$events[[1]]$time[1]
t2_on = o$events[[1]]$time[3]
t2_off = o$events[[1]]$time[4]

abline(h = 100, col = 2)
abline(v = t1_on)
abline(v = t2_on)
abline(v = t2_off)

plot(o$data$flag1, type = "l")

# event time is null, triggers only on state: seems correct
o = run_model(event_time_on = NULL, event_time_off = NULL)
o

plot(o$data$I, type = "l")
abline(h = 100, col = 2)

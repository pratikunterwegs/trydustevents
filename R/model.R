#' Run a simple SIR model
#'
#' @export
run_model <- function(
  I0 = 10,
  N = 1000,
  beta = 0.2,
  gamma = 0.1,
  nu = 0.2,
  infect_cap = 100.0,
  event_time_on = 30,
  event_time_off = 60,
  time_end = 100
) {
  # create a dust2 system
  sys <- dust2::dust_system_create(
    sirode,
    list(
      I0 = I0,
      N = N,
      beta = beta,
      gamma = gamma,
      nu = nu,
      infect_cap = infect_cap,
      event_time_on = event_time_on,
      event_time_off = event_time_off
    )
  )

  state <- c(N, I0, 0)
  flags <- c(0.0, 0.0) # flag is part of state
  dust2::dust_system_set_state(sys, c(state, flags))

  # simulate with default values
  state <- dust2::dust_system_simulate(sys, seq(0, time_end))

  # unpack the ODE state and return
  list(
    data = dust2::dust_unpack_state(sys, state),
    events = dust2::dust_system_internals(sys)[["events"]]
  )
}

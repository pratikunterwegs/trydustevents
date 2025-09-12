#' Run a simple SIR model
#'
#' @export
run_model <- function(
  I0 = 10,
  N = 1000,
  beta = 0.2,
  gamma = 0.1,
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
      infect_cap = infect_cap,
      event_time_on = event_time_on,
      event_time_off = event_time_off
    ),
    deterministic = TRUE
  )

  state <- c(N, I0, 0)
  flags <- c(0.0) # flag is part of state
  r0 <- beta / gamma
  rt <- c(r0, r0)
  dust2::dust_system_set_state(sys, c(state, rt, flags))

  # simulate with default values
  state <- dust2::dust_system_simulate(sys, seq(0, time_end, 1))

  # unpack the ODE state and return
  list(
    data = dust2::dust_unpack_state(sys, state),
    events = dust2::dust_system_internals(sys)[["events"]]
  )
}

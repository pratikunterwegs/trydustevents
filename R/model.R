#' Run a simple SIR model
#'
#' @export
run_model = function(
  I0 = 10,
  N = 1000,
  beta = 0.2,
  gamma = 0.1,
  time_end = 100
) {
  # create a dust2 system
  sys = dust2::dust_system_create(
    sirode
  )

  state = c(N, I0, 0)
  dust2::dust_system_set_state(sys, state)

  # simulate with default values
  state = dust2::dust_system_simulate(sys, seq(0, time_end))

  # unpack the ODE state and return
  dust2::dust_unpack_state(sys, state)
}

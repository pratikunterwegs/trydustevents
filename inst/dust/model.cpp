// class sirode copied and modified from
// https://mrc-ide.github.io/dust2/articles/writing.html#continuous-time-ode-models

#include <dust2/common.hpp>

// [[dust2::class(sirode)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(I0, constant = TRUE)]]
// [[dust2::parameter(N, constant = TRUE)]]
// [[dust2::parameter(beta, constant = TRUE)]]
// [[dust2::parameter(gamma, constant = TRUE)]]
// [[dust2::parameter(event_time_on, constant = TRUE)]]
// [[dust2::parameter(event_time_off, constant = TRUE)]]
class sirode {
public:
  sirode() = delete;

  using real_type = double;

  // shared model parameters
  struct shared_state {
    real_type N, I0, beta, gamma, event_time_on, event_time_off;
    size_t i_flag;
  };

  // scratch space, not used here
  struct internal_state {};

  // needed for compilation but not otherwise used
  using rng_state_type = monty::random::generator<real_type>;

  // specify how the data should be returned
  static dust2::packing packing_state(const shared_state &shared) {
    return dust2::packing{{"S", {}}, {"I", {}}, {"R", {}}, {"flag", {}}};
  }

  // pass model params from R to C++
  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = dust2::r::read_real(pars, "I0", 10.0);
    const real_type N = dust2::r::read_real(pars, "N", 1000.0);
    const real_type beta = dust2::r::read_real(pars, "beta", 0.2);
    const real_type gamma = dust2::r::read_real(pars, "gamma", 0.1);

    const real_type event_time_on = dust2::r::read_real(pars, "event_time_on");
    const real_type event_time_off =
        dust2::r::read_real(pars, "event_time_off");

    // index of the flag
    const size_t i_flag = 3; // hard-coded, comes after S, I, R, 0-indexed

    return shared_state{N,     I0, beta, gamma, event_time_on, event_time_off,
                        i_flag};
  }

  static void update_shared(cpp11::list pars, const shared_state &shared) {
    // NOTE: we are setting these constant
  }

  // pass the initial state from R to C++/dust2
  static void initial(real_type time, const shared_state &shared,
                      internal_state &internal, rng_state_type &rng_state,
                      real_type *state_next) {
    state_next[0] = shared.N - shared.I0;
    state_next[1] = shared.I0;
    state_next[2] = 0.0;
    state_next[3] = 0.0; // events initialised as off
  }

  // RHS of ODE
  static void rhs(real_type time, const real_type *state,
                  const shared_state &shared, internal_state &internal,
                  real_type *state_deriv) {
    const auto S = state[0];
    const auto I = state[1];

    // when the event (NPI) is on, reduce beta to 50%
    const double beta_now = shared.beta * (1.0 - 0.5 * state[shared.i_flag]);

    const auto rate_SI = beta_now * S * I / shared.N;
    const auto rate_IR = shared.gamma * I;

    state_deriv[0] = -rate_SI;
    state_deriv[1] = rate_SI - rate_IR;
    state_deriv[2] = rate_IR;
  }

  // model events; this uses manual dust2 events
  static auto events(const shared_state &shared,
                     const internal_state &internal) {
    /* create a time-limited event that tests for a time-on and time-off */
    // make a function that tests for time
    auto fn_time_test_on = [&](const double t, const double *y) {
      return t - shared.event_time_on;
    };

    auto fn_state_test = [&](const double t, const double *y) {
      const double infect_cap = 100.0; // an imaginary infectious capacity
      return y[0] - infect_cap;
    };

    auto fn_time_test_off = [&](const double t, const double *y) {
      return t - shared.event_time_off;
    };

    // make a function that changes a flag value to on
    auto fn_flag_on = [&](const double t, const double sign, double *y) {
      y[shared.i_flag] = 1.0;
    };

    // function to switch flag off
    auto fn_flag_off = [&](const double t, const double sign, double *y) {
      y[shared.i_flag] = 0.0;
    };

    // a dummy function
    auto fn_dummy_action = [&](const double t, const double sign, double *y) {};

    // make on and off events
    std::string name_event_on = "event_on";
    dust2::ode::event<real_type> event_time_on(name_event_on, {},
                                               fn_time_test_on, fn_flag_on,
                                               dust2::ode::root_type::increase);

    std::string name_event_off = "event_off";
    dust2::ode::event<real_type> event_time_off(
        name_event_off, {}, fn_time_test_off, fn_flag_off,
        dust2::ode::root_type::decrease);

    // an event that launches when 100 individuals are in I
    std::string name_state_on = "event_state_on";
    dust2::ode::event<real_type> event_state_on(
        name_state_on, {1}, fn_state_test, fn_dummy_action,
        dust2::ode::root_type::increase);

    std::string name_state_off = "event_state_off";
    dust2::ode::event<real_type> event_state_off(
        name_state_off, {1}, fn_state_test, fn_dummy_action,
        dust2::ode::root_type::decrease);

    // return events vector
    return dust2::ode::events_type<real_type>(
        {event_time_on, event_time_off, event_state_on, event_state_off});
  }
};

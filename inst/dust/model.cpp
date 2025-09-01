// class sirode copied and modified from
// https://mrc-ide.github.io/dust2/articles/writing.html#continuous-time-ode-models

#include <algorithm>
#include <dust2/common.hpp>

const int state_size = 6L;

// [[dust2::class(sirode)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(I0, constant = TRUE)]]
// [[dust2::parameter(N, constant = TRUE)]]
// [[dust2::parameter(beta, constant = TRUE)]]
// [[dust2::parameter(gamma, constant = TRUE)]]
// [[dust2::parameter(infect_cap, constant = TRUE)]]
// [[dust2::parameter(event_time_on, constant = TRUE)]]
// [[dust2::parameter(event_time_off, constant = TRUE)]]
class sirode {
public:
  sirode() = delete;

  using real_type = double;

  // shared model parameters
  struct shared_state {
    real_type N, I0, beta, gamma, nu, infect_cap, event_time_on, event_time_off;
    size_t i_flag;
  };

  // scratch space, not used here
  struct internal_state {};

  // needed for compilation but not otherwise used
  using rng_state_type = monty::random::generator<real_type>;

  // specify how the data should be returned
  static dust2::packing packing_state(const shared_state &shared) {
    return dust2::packing{{"S", {}},        {"I", {}},      {"R", {}},
                          {"rt_cumul", {}}, {"rt_est", {}}, {"rt_saved", {}},
                          {"flag", {}}};
  }

  static size_t size_special() {
    return 2; // rt_saved and flags
  }

  // pass model params from R to C++
  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = dust2::r::read_real(pars, "I0", 10.0);
    const real_type N = dust2::r::read_real(pars, "N", 1000.0);
    const real_type beta = dust2::r::read_real(pars, "beta", 0.2);
    const real_type gamma = dust2::r::read_real(pars, "gamma", 0.1);

    // trigger NPI on 100 infections by default
    const real_type infect_cap = dust2::r::read_real(pars, "infect_cap", 100.0);

    const real_type event_time_on =
        dust2::r::read_real(pars, "event_time_on", NA_REAL);
    const real_type event_time_off =
        dust2::r::read_real(pars, "event_time_off", NA_REAL);

    // index of the flag
    const size_t i_flag =
        state_size; // hard-coded, comes after S, I, R, Rts, 0-indexed

    // clang-format off
    return shared_state{
      N, I0, beta, gamma, infect_cap, event_time_on, event_time_off,
      i_flag
    };
    // clang-format on
  }

  static void update_shared(cpp11::list pars, const shared_state &shared) {
    // NOTE: we are setting these constant
  }

  // pass the initial state from R to C++/dust2
  static void initial(real_type time, const shared_state &shared,
                      internal_state &internal, rng_state_type &rng_state,
                      real_type *state_next) {
    const double p_susc = (shared.N - shared.I0) / shared.N;
    const double rt = shared.beta * p_susc / shared.gamma;

    state_next[0] = shared.N - shared.I0;
    state_next[1] = shared.I0;
    state_next[2] = 0.0;

    state_next[3] = rt; // rt_cumul
    state_next[4] = rt; // rt_est
    // state_next[5] = rt; // rt_saved

    // flags init as 0
    state_next[6] = 0.0;
  }

  // RHS of ODE
  static void rhs(real_type time, const real_type *state,
                  const shared_state &shared, internal_state &internal,
                  real_type *state_deriv) {
    const real_type S = state[0];
    const real_type I = state[1];

    const real_type rate_SI = shared.beta * S * I / shared.N;
    const real_type rate_IR = shared.gamma * I;

    const double p_susc = S / shared.N;
    const double rt = shared.beta * p_susc / shared.gamma;

    state_deriv[0] = -rate_SI;
    state_deriv[1] = rate_SI - rate_IR;
    state_deriv[2] = rate_IR;
    state_deriv[3] = rt; // cumulative rt
  }

  // model events; this uses manual dust2 events
  static auto events(const shared_state &shared,
                     const internal_state &internal) {
    // make a function that tests for time as integerish
    auto fn_time_test_integer = [&](const double t, const double *y) {
      const double time_mod = std::fmod(t, 1.0);
      const bool is_t_integerish = time_mod < 1e-3;

      if (is_t_integerish) {
        return 0.0;
      } else {
        return 1.0;
      }
    };

    // check the first element of the array pointer y against infect_cap
    // NOTE: the way dust2 events work is that the condition function works on
    // a pointer to an array comprising of user-specified indices (how does it
    // handle non-contiguous indices - unclear), so the starting index i0 is
    // accessed in the condition function as 0.
    // NOTE: the action function does not work this way, it accesses the whole
    // state array.
    auto fn_state_test_on = [&](const double t, const double *y) {
      const double sum_infected = std::accumulate(y, y + 1, 0);
      Rprintf("sum_infected = %f \n", sum_infected);
      return sum_infected - shared.infect_cap;
    };

    auto fn_state_test_off = [&](const double t, const double *y) {
      const double diff = y[0] - 1.0; // rt < 1.0
      return diff < 0.0 ? 0.0 : 1.0;
    };

    // make a function that changes a flag value to on
    auto fn_flag_on = [&](const double t, const double sign, double *y) {
      y[shared.i_flag] = 1.0;
    };

    // function to switch flag off
    auto fn_flag_off = [&](const double t, const double sign, double *y) {
      y[shared.i_flag] = 0.0;
    };

    // a dummy function, not used
    auto fn_dummy_action = [&](const double t, const double sign, double *y) {
      Rprintf("Action at time = %f \n", t);
    };

    // save rt estimated to r_saved, refer to absolute positions in y
    auto fn_save_internal = [&](const double t, const double sign, double *y) {
      y[5] = y[4];
    };

    // an event that launches when 100 individuals are in I
    std::string name_state_on = "event_state_on";
    dust2::ode::event<real_type> event_state_on(name_state_on, {1},
                                                fn_state_test_on, fn_flag_on);

    // even that ends when Rt_saved hits 1 while decreasing
    std::string name_state_off = "event_state_off";
    dust2::ode::event<real_type> event_state_off(
        name_state_off, {5}, fn_state_test_off, fn_flag_off);

    // event that triggers when time is integerish
    std::string name_time_integer = "event_time_integer";
    dust2::ode::event<real_type> event_time_integer(
        name_time_integer, {}, fn_time_test_integer, fn_save_internal);

    // return events vector
    return dust2::ode::events_type<real_type>(
        {event_state_on, event_state_off});
  }

  static size_t size_output() { return 1; }

  static auto delays(const shared_state &shared) {
    return dust2::ode::delays<real_type>({{1, {{3, 1}}}});
  }

  static void output(real_type time, real_type *state,
                     const shared_state &shared, internal_state &internal,
                     const dust2::ode::delay_result_type<real_type> &delays) {
    const auto &delay_rt = delays[0].data[0];
    state[4] = state[3] - delay_rt;
  }

  static auto zero_every(const shared_state &shared) {
    return dust2::zero_every_type<real_type>{{1, {4}}}; // zero IPR value
  }
};

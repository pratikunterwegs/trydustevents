// class sirode copied and modified from
// https://mrc-ide.github.io/dust2/articles/writing.html#continuous-time-ode-models

#include <dust2/common.hpp>

// [[dust2::class(sirode)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(I0, constant = TRUE)]]
// [[dust2::parameter(N, constant = TRUE)]]
// [[dust2::parameter(beta, constant = TRUE)]]
// [[dust2::parameter(gamma, constant = TRUE)]]
class sirode {
 public:
  sirode() = delete;

  using real_type = double;

  // shared model parameters
  struct shared_state {
    real_type N;
    real_type I0;
    real_type beta;
    real_type gamma;
  };

  // scratch space, not used here
  struct internal_state {};

  // needed for compilation but not otherwise used
  using rng_state_type = monty::random::generator<real_type>;

  // specify how the data should be returned
  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"S", {}}, {"I", {}}, {"R", {}}};
  }

  // pass model params from R to C++
  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = dust2::r::read_real(pars, "I0", 10);
    const real_type N = dust2::r::read_real(pars, "N", 1000);
    const real_type beta = dust2::r::read_real(pars, "beta", 0.2);
    const real_type gamma = dust2::r::read_real(pars, "gamma", 0.1);
    return shared_state{N, I0, beta, gamma};
  }

  static void update_shared(cpp11::list pars, const shared_state& shared) {
    // NOTE: we are setting these constant
  }

  // pass the initial state from R to C++/dust2
  static void initial(real_type time, const shared_state& shared,
                      internal_state& internal, rng_state_type& rng_state,
                      real_type* state_next) {
    state_next[0] = shared.N - shared.I0;
    state_next[1] = shared.I0;
    state_next[2] = 0;
  }

  // RHS of ODE
  static void rhs(real_type time, const real_type* state,
                  const shared_state& shared, internal_state& internal,
                  real_type* state_deriv) {
    const auto S = state[0];
    const auto I = state[1];

    const auto rate_SI = shared.beta * S * I / shared.N;
    const auto rate_IR = shared.gamma * I;

    state_deriv[0] = -rate_SI;
    state_deriv[1] = rate_SI - rate_IR;
    state_deriv[2] = rate_IR;
  }
};

// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// model.cpp
SEXP dust2_system_sirode_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::list r_time_control, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic, cpp11::sexp r_n_threads);
extern "C" SEXP _trydustevents_dust2_system_sirode_alloc(SEXP r_pars, SEXP r_time, SEXP r_time_control, SEXP r_n_particles, SEXP r_n_groups, SEXP r_seed, SEXP r_deterministic, SEXP r_n_threads) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_time), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_time_control), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_particles), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_groups), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_seed), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_deterministic), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_threads)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_internals(cpp11::sexp ptr, bool include_coefficients, bool include_history);
extern "C" SEXP _trydustevents_dust2_system_sirode_internals(SEXP ptr, SEXP include_coefficients, SEXP include_history) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_internals(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(include_coefficients), cpp11::as_cpp<cpp11::decay_t<bool>>(include_history)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time);
extern "C" SEXP _trydustevents_dust2_system_sirode_run_to_time(SEXP ptr, SEXP r_time) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_run_to_time(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_time)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_state(cpp11::sexp ptr, cpp11::sexp r_index_state, cpp11::sexp r_index_particle, cpp11::sexp r_index_group, bool preserve_particle_dimension, bool preserve_group_dimension);
extern "C" SEXP _trydustevents_dust2_system_sirode_state(SEXP ptr, SEXP r_index_state, SEXP r_index_particle, SEXP r_index_group, SEXP preserve_particle_dimension, SEXP preserve_group_dimension) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_state(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index_state), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index_particle), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index_group), cpp11::as_cpp<cpp11::decay_t<bool>>(preserve_particle_dimension), cpp11::as_cpp<cpp11::decay_t<bool>>(preserve_group_dimension)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_time(cpp11::sexp ptr);
extern "C" SEXP _trydustevents_dust2_system_sirode_time(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_time(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_set_state_initial(cpp11::sexp ptr);
extern "C" SEXP _trydustevents_dust2_system_sirode_set_state_initial(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_set_state_initial(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_set_state(cpp11::sexp ptr, cpp11::list r_state);
extern "C" SEXP _trydustevents_dust2_system_sirode_set_state(SEXP ptr, SEXP r_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_set_state(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_state)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_reorder(cpp11::sexp ptr, cpp11::integers r_index);
extern "C" SEXP _trydustevents_dust2_system_sirode_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_reorder(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::integers>>(r_index)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_rng_state(cpp11::sexp ptr);
extern "C" SEXP _trydustevents_dust2_system_sirode_rng_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_rng_state(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state);
extern "C" SEXP _trydustevents_dust2_system_sirode_set_rng_state(SEXP ptr, SEXP r_rng_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_set_rng_state(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_rng_state)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_set_time(cpp11::sexp ptr, cpp11::sexp r_time);
extern "C" SEXP _trydustevents_dust2_system_sirode_set_time(SEXP ptr, SEXP r_time) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_set_time(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_time)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_update_pars(cpp11::sexp ptr, cpp11::list pars);
extern "C" SEXP _trydustevents_dust2_system_sirode_update_pars(SEXP ptr, SEXP pars) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_update_pars(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(pars)));
  END_CPP11
}
// model.cpp
SEXP dust2_system_sirode_simulate(cpp11::sexp ptr, cpp11::sexp r_times, cpp11::sexp r_index_state, bool preserve_particle_dimension, bool preserve_group_dimension);
extern "C" SEXP _trydustevents_dust2_system_sirode_simulate(SEXP ptr, SEXP r_times, SEXP r_index_state, SEXP preserve_particle_dimension, SEXP preserve_group_dimension) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust2_system_sirode_simulate(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_times), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index_state), cpp11::as_cpp<cpp11::decay_t<bool>>(preserve_particle_dimension), cpp11::as_cpp<cpp11::decay_t<bool>>(preserve_group_dimension)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_trydustevents_dust2_system_sirode_alloc",             (DL_FUNC) &_trydustevents_dust2_system_sirode_alloc,             8},
    {"_trydustevents_dust2_system_sirode_internals",         (DL_FUNC) &_trydustevents_dust2_system_sirode_internals,         3},
    {"_trydustevents_dust2_system_sirode_reorder",           (DL_FUNC) &_trydustevents_dust2_system_sirode_reorder,           2},
    {"_trydustevents_dust2_system_sirode_rng_state",         (DL_FUNC) &_trydustevents_dust2_system_sirode_rng_state,         1},
    {"_trydustevents_dust2_system_sirode_run_to_time",       (DL_FUNC) &_trydustevents_dust2_system_sirode_run_to_time,       2},
    {"_trydustevents_dust2_system_sirode_set_rng_state",     (DL_FUNC) &_trydustevents_dust2_system_sirode_set_rng_state,     2},
    {"_trydustevents_dust2_system_sirode_set_state",         (DL_FUNC) &_trydustevents_dust2_system_sirode_set_state,         2},
    {"_trydustevents_dust2_system_sirode_set_state_initial", (DL_FUNC) &_trydustevents_dust2_system_sirode_set_state_initial, 1},
    {"_trydustevents_dust2_system_sirode_set_time",          (DL_FUNC) &_trydustevents_dust2_system_sirode_set_time,          2},
    {"_trydustevents_dust2_system_sirode_simulate",          (DL_FUNC) &_trydustevents_dust2_system_sirode_simulate,          5},
    {"_trydustevents_dust2_system_sirode_state",             (DL_FUNC) &_trydustevents_dust2_system_sirode_state,             6},
    {"_trydustevents_dust2_system_sirode_time",              (DL_FUNC) &_trydustevents_dust2_system_sirode_time,              1},
    {"_trydustevents_dust2_system_sirode_update_pars",       (DL_FUNC) &_trydustevents_dust2_system_sirode_update_pars,       2},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_trydustevents(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

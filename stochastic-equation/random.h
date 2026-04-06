//GENERADOR DE NUMEROS ALEATORIOS
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define RANDOM gsl_rng_uniform(gsl_rng_r)
#define RANDOM_INT(A) gsl_rng_uniform_int(gsl_rng_r, A)
#define RANDOM_GAUSS(S) gsl_ran_gaussian(gsl_rng_r, S)
#define RANDOM_EXP(MU) gsl_ran_exponential(gsl_rng_r, MU)
#define RANDOM_PARETO(A,B) gsl_ran_pareto(gsl_rng_r, A, B)
// semilla del reloj
#define INIT_RANDOM {gsl_rng_env_setup(); if(!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0); gsl_rng_T=gsl_rng_default;  gsl_rng_r=gsl_rng_alloc(gsl_rng_T);}
// semilla fija
//#define INITIALIZE_RANDOM {gsl_rng_env_setup(); gsl_rng_T=gsl_rng_default;  gsl_rng_r=gsl_rng_alloc(gsl_rng_T);}
#define FREE_RANDOM gsl_rng_free(gsl_rng_r);
const gsl_rng_type * gsl_rng_T;
gsl_rng * gsl_rng_r;


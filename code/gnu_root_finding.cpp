#include "math.h"
#include <chrono>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

struct functional_params {
  vector<double> theta;
};

void print_state(size_t iter, gsl_multiroot_fdfsolver *s) {
  cout << "iter = " << iter << " x = " << gsl_vector_get(s->x, 0) << " "
       << gsl_vector_get(s->x, 1)
       << " "
          "f(x) = "
       << gsl_vector_get(s->f, 0) << " " << gsl_vector_get(s->f, 1) << endl;
}

int functional_f(const gsl_vector *x, void *p, gsl_vector *f) {
  struct functional_params *params = (struct functional_params *)p;
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  double f0 = 0, f1 = 0;

  for (int i = 0; i < 4; i++) {
    f0 += sin(x0 - params->theta[i * 3]) +
          sin(x0 + x1 + params->theta[i * 3 + 2]) +
          sin(x0 + params->theta[(i + 4) * 3]) +
          sin(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    f1 += sin(x1 - params->theta[i * 3 + 1]) +
          sin(x0 + x1 + params->theta[i * 3 + 2]) +
          sin(x1 + params->theta[(i + 4) * 3 + 1]) +
          sin(x0 + x1 - params->theta[(i + 4) * 3 + 2]);
  }

  gsl_vector_set(f, 0, f0);
  gsl_vector_set(f, 1, f1);

  return GSL_SUCCESS;
}

int functional_df(const gsl_vector *x, void *p, gsl_matrix *J) {
  struct functional_params *params = (struct functional_params *)p;
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  double J00 = 0, J01 = 0, J10 = 0, J11 = 0;

  for (int i = 0; i < 4; i++) {
    J00 += cos(x0 - params->theta[i * 3]) +
           cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x0 + params->theta[(i + 4) * 3]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    J01 += cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    J10 += cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    J11 += cos(x1 - params->theta[i * 3 + 1]) +
           cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x1 + params->theta[(i + 4) * 3 + 1]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);
  }

  gsl_matrix_set(J, 0, 0, J00);
  gsl_matrix_set(J, 0, 1, J01);
  gsl_matrix_set(J, 1, 0, J10);
  gsl_matrix_set(J, 1, 1, J11);

  return GSL_SUCCESS;
}

int functional_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J) {
  struct functional_params *params = (struct functional_params *)p;
  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  double f0 = 0, f1 = 0;
  double J00 = 0, J01 = 0, J10 = 0, J11 = 0;

  for (int i = 0; i < 4; i++) {
    f0 += sin(x0 - params->theta[i * 3]) +
          sin(x0 + x1 + params->theta[i * 3 + 2]) +
          sin(x0 + params->theta[(i + 4) * 3]) +
          sin(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    f1 += sin(x1 - params->theta[i * 3 + 1]) +
          sin(x0 + x1 + params->theta[i * 3 + 2]) +
          sin(x1 + params->theta[(i + 4) * 3 + 1]) +
          sin(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    J00 += cos(x0 - params->theta[i * 3]) +
           cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x0 + params->theta[(i + 4) * 3]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    J01 += cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    J10 += cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);

    J11 += cos(x1 - params->theta[i * 3 + 1]) +
           cos(x0 + x1 + params->theta[i * 3 + 2]) +
           cos(x1 + params->theta[(i + 4) * 3 + 1]) +
           cos(x0 + x1 - params->theta[(i + 4) * 3 + 2]);
  }

  gsl_vector_set(f, 0, f0);
  gsl_vector_set(f, 1, f1);

  gsl_matrix_set(J, 0, 0, J00);
  gsl_matrix_set(J, 0, 1, J01);
  gsl_matrix_set(J, 1, 0, J10);
  gsl_matrix_set(J, 1, 1, J11);

  return GSL_SUCCESS;
}

int main() {

  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 2;
  struct functional_params p = {
      {-0.003430567378423177, 2.5187546645278616,   -0.18790892227498368,
       1.8580870595325507,    0.5462893265144331,   -0.11351799849439281,
       0.5752081963953266,    2.161798709021106,    -0.10056541433062449,
       2.4770558204233915,    0.014931896904666786, 2.1720296730661524,
       -1.2834853852652746,   0.7885560009794599,   1.1986857960623745,
       2.7952004417190537,    2.930628506228693,    0.48003331213481815,
       0.527335578519212,     2.297331991870097,    2.3595674142926963,
       -3.109486430507647,    -2.4471192629700935,  -1.1906340553863846}};

  gsl_multiroot_function_fdf f = {&functional_f, &functional_df,
                                  &functional_fdf, n, &p};

  double x_init[2] = {0.123, 0.4363};
  gsl_vector *x = gsl_vector_alloc(n);

  gsl_vector_set(x, 0, x_init[0]);
  gsl_vector_set(x, 1, x_init[1]);

  T = gsl_multiroot_fdfsolver_hybridsj;
  s = gsl_multiroot_fdfsolver_alloc(T, n);
  gsl_multiroot_fdfsolver_set(s, &f, x);

  print_state(iter, s);

  gsl_vector *x_tmp = gsl_vector_alloc(n);

  start_time = clock();

  for (int i = 0; i < 1000000; i++) {
    do {
      iter++;

      status = gsl_multiroot_fdfsolver_iterate(s);

      // x_tmp = gsl_multiroot_fdfsolver_root(s);
      // cout << gsl_vector_get(x_tmp, 0) << " " << gsl_vector_get(x_tmp, 1)
      //      << endl;
      // print_state(iter, s);

      if (status)
        break;

      status = gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);
  }

  end_time = clock();
  search_time = end_time - start_time;
  cout << "roots time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  printf("status = %s\n", gsl_strerror(status));

  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
}
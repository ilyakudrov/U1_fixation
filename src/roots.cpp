#include "../include/roots.h"

using namespace std;

double function(vector<double> &x, const functional_params *params) {
  double result = 0;
  for (int i = 0; i < 4; i++) {
    result += cos(x[0] - params->theta[i * 3]) +
              cos(x[1] - params->theta[i * 3 + 1]) +
              cos(x[0] + x[1] + params->theta[i * 3 + 2]) +
              cos(x[0] + params->theta[(i + 4) * 3]) +
              cos(x[1] + params->theta[(i + 4) * 3 + 1]) +
              cos(x[0] + x[1] - params->theta[(i + 4) * 3 + 2]);
  }

  return result;
}

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

vector<double> find_root_max(gsl_multiroot_fdfsolver *s,
                             gsl_multiroot_function_fdf *f,
                             vector<vector<double>> &x_init,
                             const functional_params *params,
                             vector<double> &x_former) {
  int status;
  size_t iter = 0;

  vector<double> root_max = x_former;
  double function_tmp;
  double function_max = function(x_former, params);
  // cout << "functional former in find_root_max "
  //      << " " << function_max << endl;

  gsl_vector *x = gsl_vector_alloc(2);
  gsl_vector *x_tmp;

  vector<double> x_tmp1;

  // cout << "params in find_root_max:" << endl;
  // for (int i = 0; i < 24; i++) {
  //   cout << params->theta[i] << " ";
  // }
  // cout << endl;

  for (int i = 0; i < x_init.size(); i++) {
    gsl_vector_set(x, 0, x_init[i][0]);
    gsl_vector_set(x, 1, x_init[i][1]);
    gsl_multiroot_fdfsolver_set(s, f, x);

    do {
      iter++;

      status = gsl_multiroot_fdfsolver_iterate(s);

      // x_tmp = gsl_multiroot_fdfsolver_root(s);
      // cout << gsl_vector_get(x_tmp, 0) << " " << gsl_vector_get(x_tmp, 1) <<
      // endl; print_state(iter, s);

      if (status)
        break;

      status = gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);

    x_tmp = gsl_multiroot_fdfsolver_root(s);
    x_tmp1 = {gsl_vector_get(x_tmp, 0), gsl_vector_get(x_tmp, 1)};

    function_tmp = function(x_tmp1, params);

    // cout << "roots: " << x_tmp1[0] << " " << x_tmp1[1] << " functional "
    //      << function_tmp << endl;

    if (function_tmp > function_max) {
      function_max = function_tmp;
      root_max = x_tmp1;
    }
  }

  return root_max;
}
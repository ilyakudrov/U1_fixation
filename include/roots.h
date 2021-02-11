#pragma once

#include "math.h"
#include <algorithm>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <vector>

using namespace std;

struct functional_params {
  vector<double> theta;
};

double function(vector<double> &x, const functional_params *params);

void print_state(size_t iter, gsl_multiroot_fdfsolver *s);

int functional_f(const gsl_vector *x, void *p, gsl_vector *f);

int functional_df(const gsl_vector *x, void *p, gsl_matrix *J);

int functional_fdf(const gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J);

vector<double> find_root_max(gsl_multiroot_fdfsolver *s,
                             gsl_multiroot_function_fdf *f,
                             vector<vector<double>> &x_init,
                             const functional_params *params,
                             vector<double> &x_former);
#include "../include/maximization.h"

// iteration over lattice
#define SPACE_ITER_START                                                       \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {                                     \
          link.go(x, y, z, t);                                                 \
          link.update(0);                                                      \
          link.update(1);                                                      \
          link.update(2);                                                      \
          link.update(3);

#define SPACE_ITER_END                                                         \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

using namespace std;

void claculate_params(vector<angles> &angles, vector<gauge> &gauge, link1 &link,
                      functional_params &params) {

  for (int mu = 0; mu < 4; mu++) {

    params.theta[3 * mu] = angles[link.place + mu].angle[0];
    params.theta[3 * mu + 1] = angles[link.place + mu].angle[1];
    params.theta[3 * mu + 2] = angles[link.place + mu].angle[2];

    link.move(mu, 1);

    params.theta[3 * mu] =
        params.theta[3 * mu] + gauge[link.place / 4].alpha[0];
    params.theta[3 * mu + 1] =
        params.theta[3 * mu + 1] + gauge[link.place / 4].alpha[1];
    params.theta[3 * mu + 2] = params.theta[3 * mu + 2] -
                               gauge[link.place / 4].alpha[0] -
                               gauge[link.place / 4].alpha[1];

    link.move(mu, -2);

    params.theta[3 * mu + 12] =
        angles[link.place + mu].angle[0] - gauge[link.place / 4].alpha[0];
    params.theta[3 * mu + 13] =
        angles[link.place + mu].angle[1] - gauge[link.place / 4].alpha[1];
    params.theta[3 * mu + 14] = angles[link.place + mu].angle[2] +
                                gauge[link.place / 4].alpha[0] +
                                gauge[link.place / 4].alpha[1];

    link.move(mu, 1);
  }
}

double functional_test(vector<angles> &angles, vector<gauge> &gauge) {
  double functional = 0;

  functional_params params;
  params.theta.reserve(24);
  params.theta.resize(24);

  vector<double> x_tmp;

  link1 link(x_size, y_size, z_size, t_size);

  SPACE_ITER_START

  // if (x == 0 && y == 0 && z == 0 && t == 0)
  //   cout << "gauge at 0 " << gauge[link.place / 4].alpha[0] << " "
  //        << gauge[link.place / 4].alpha[1] << endl;

  claculate_params(angles, gauge, link, params);

  // if (x == 0 && y == 0 && z == 0 && t == 0)
  //   cout << "link in functional_test " << link;

  // if (x == 0 && y == 0 && z == 0 && t == 0) {
  //   cout << "params in functional_test:" << endl;
  //   for (int i = 0; i < 24; i++) {
  //     cout << params.theta[i] << " ";
  //   }
  //   cout << endl;
  // }

  // if (x == 0 && y == 0 && z == 0 && t == 0)
  //   cout << "functional in functional_test before " << functional << endl;

  x_tmp = {gauge[link.place / 4].alpha[0], gauge[link.place / 4].alpha[1]};

  functional += function(x_tmp, &params);

  // if (x == 0 && y == 0 && z == 0 && t == 0)
  //   cout << "roots in functional_test: " << x_tmp[0] << " " << x_tmp[1]
  //        << " functional in functiona_test " << functional << endl;

  // if (x == 0 && y == 0 && z == 0 && t == 0)
  //   cout << "functional at 0 " << functional << endl;

  SPACE_ITER_END

  return functional / (x_size * y_size * z_size * t_size * 8);
}

double functional(vector<angles> &angles, vector<gauge> &gauge) {
  double functional = 0;
  link1 link(x_size, y_size, z_size, t_size);

  vector<double> theta(3);

  SPACE_ITER_START

  // if (x == 0 && y == 0 && z == 0 && t == 0)
  //   cout << "gauge at 0 " << gauge[link.place / 4].alpha[0] << " "
  //        << gauge[link.place / 4].alpha[1] << endl;

  for (int mu = 0; mu < 4; mu++) {
    theta[0] =
        angles[link.place + mu].angle[0] - gauge[link.place / 4].alpha[0];
    theta[1] =
        angles[link.place + mu].angle[1] - gauge[link.place / 4].alpha[1];
    theta[2] = angles[link.place + mu].angle[2] -
               gauge[link.place / 4].alpha[0] - gauge[link.place / 4].alpha[1];

    link.move(mu, 1);

    theta[0] = theta[0] + gauge[link.place / 4].alpha[0];
    theta[1] = theta[1] + gauge[link.place / 4].alpha[1];
    theta[2] = theta[2] + gauge[link.place / 4].alpha[0] +
               gauge[link.place / 4].alpha[1];

    link.move(mu, -1);

    functional += cos(theta[0]) + cos(theta[1]) + cos(theta[2]);
  }

  // if (x == 0 && y == 0 && z == 0 && t == 0)
  //   cout << "functional at 0 " << functional << endl;

  SPACE_ITER_END

  return functional / (x_size * y_size * z_size * t_size * 4);
}

vector<maximize_data> maximize(vector<angles> &angles, vector<gauge> &gauge) {

  const size_t n = 2;

  functional_params params;
  params.theta.reserve(24);
  params.theta.resize(24);

  vector<vector<double>> x_init_default;
  int discretization = 4;
  for (int i = 0; i <= discretization; i++) {
    for (int j = 0; j <= discretization; j++) {
      x_init_default.push_back(
          {2 * M_PI * i / discretization, 2 * M_PI * j / discretization});
    }
  }
  vector<double> x_former;

  int iteration_number = 100;
  vector<maximize_data> maximization_data;
  maximization_data.reserve(iteration_number);
  maximize_data data_tmp;

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;

  T = gsl_multiroot_fdfsolver_hybridsj;
  s = gsl_multiroot_fdfsolver_alloc(T, n);

  gsl_multiroot_function_fdf f = {&functional_f, &functional_df,
                                  &functional_fdf, n, &params};

  link1 link(x_size, y_size, z_size, t_size);

  vector<double> root;

  data_tmp.iteration = 0;
  data_tmp.functional = functional_test(angles, gauge);
  maximization_data.push_back(data_tmp);

  for (int i = 0; i < iteration_number; i++) {

    SPACE_ITER_START

    claculate_params(angles, gauge, link, params);

    // TODO: fix root returning later
    x_former = {gauge[link.place / 4].alpha[0], gauge[link.place / 4].alpha[1]};
    root = find_root_max(s, &f, x_init_default, &params, x_former);

    gauge[link.place / 4].alpha[0] = root[0];
    gauge[link.place / 4].alpha[1] = root[1];

    SPACE_ITER_END

    data_tmp.iteration = i + 1;
    data_tmp.functional = functional_test(angles, gauge);
    maximization_data.push_back(data_tmp);

    // cout << "functional " << functional(angles, gauge) << endl;
    // cout << "functional_test " << functional_test(angles, gauge) << endl;
  }

  return maximization_data;
}

void maximize_test(vector<angles> &angles, vector<gauge> &gauge) {

  const size_t n = 2;

  functional_params params;
  params.theta.reserve(24);
  params.theta.resize(24);

  vector<vector<double>> x_init_default;
  int discretization = 4;
  for (int i = 0; i <= discretization; i++) {
    for (int j = 0; j <= discretization; j++) {
      x_init_default.push_back(
          {2 * M_PI * i / discretization, 2 * M_PI * j / discretization});
    }
  }
  vector<double> x_former;

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;

  T = gsl_multiroot_fdfsolver_hybridsj;
  s = gsl_multiroot_fdfsolver_alloc(T, n);

  gsl_multiroot_function_fdf f = {&functional_f, &functional_df,
                                  &functional_fdf, n, &params};

  link1 link(x_size, y_size, z_size, t_size);

  vector<double> root;

  // for (int i = 0; i < 100; i++) {

  // SPACE_ITER_START

  claculate_params(angles, gauge, link, params);

  cout << "functional before " << functional(angles, gauge) << endl;
  cout << "functional_test before " << functional_test(angles, gauge) << endl;

  // TODO: fix root returning later
  x_former = {gauge[link.place / 4].alpha[0], gauge[link.place / 4].alpha[1]};
  root = find_root_max(s, &f, x_init_default, &params, x_former);

  cout << "final roots " << root[0] << " " << root[1] << endl;

  cout << link << endl;

  gauge[link.place / 4].alpha[0] = root[0];
  gauge[link.place / 4].alpha[1] = root[1];

  // SPACE_ITER_END

  cout << "functional after " << functional(angles, gauge) << endl;
  cout << "functional_test after " << functional_test(angles, gauge) << endl;
  // }
}
#include "../include/link.h"
#include "../include/matrix.h"
#include "../include/maximization.h"
#include <fstream>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <vector>

using namespace std;

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

vector<angles> read_double_fortran(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  vector<double> v(data_size * 4);
  vector<angles> array(data_size);

  ifstream stream(file_name);

  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size * 4 - 1) * sizeof(double)))
    cout << "read_double_fortran error: " << file_name << endl;

  link1 link(x_size, y_size, z_size, t_size);
  angles angles_tmp;

  for (int mu = 0; mu < 4; mu++) {

    SPACE_ITER_START

    angles_tmp.angle[0] = v[mu * data_size + link.place / 2];
    angles_tmp.angle[1] = v[mu * data_size + link.place / 2 + data_size / 2];
    angles_tmp.angle[2] = 0 - angles_tmp.angle[0] - angles_tmp.angle[1];
    array[link.place + mu] = angles_tmp;

    SPACE_ITER_END
  }
  stream.close();
  return array;
}

int x_size;
int y_size;
int z_size;
int t_size;

int main() {

  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 16;
  y_size = 16;
  z_size = 16;
  t_size = 16;

  cout.precision(10);

  string file_path = "../conf/bqcd_ab.b5p25.k1360501002.lat";

  vector<angles> data = read_double_fortran(file_path);

  // link1 link(x_size, y_size, z_size, t_size);
  // for (int x = 0; x < 16; x++) {
  //   link.go(x, 0, 0, 0);
  //   link.update(0);
  //   for (int mu = 0; mu < 4; mu++) {
  //     // cout << link;
  //     cout << x << " " << mu << " " << data[link.place + mu].angle[0] << " "
  //          << data[link.place + mu].angle[1] << " "
  //          << data[link.place + mu].angle[2] << endl;
  //   }
  // }

  vector<gauge> gauge_conf(x_size * y_size * z_size * t_size);

  for (int i = 0; i < gauge_conf.size(); i++) {
    gauge_conf[i].alpha[0] = 0;
    gauge_conf[i].alpha[1] = 0;
  }

  vector<maximize_data> maximization_data;

  maximization_data = maximize(data, gauge_conf);

  ofstream output;
  output.precision(17);
  output.open("../tuning/data/"
              "maximization_test.csv");

  for (int i = 0; i < maximization_data.size(); i++) {
    output << maximization_data[i].iteration << ", "
           << maximization_data[i].functional << endl;
  }

  output.close();

  // cout << "functional " << functional(data, gauge_conf) << endl;
  // cout << "functional_test " << functional_test(data, gauge_conf) << endl;

  // maximize_test(data, gauge_conf);
}
#include "../include/link.h"
// #include "../include/data.h"
#include <cmath>

#define data_size                                                              \
  4 * lattice_size[0] * lattice_size[1] * lattice_size[2] * lattice_size[3]
#define PLACE3_LINK_DIR                                                        \
  (coordinate[3]) * 3 * lattice_size[0] * lattice_size[1] * lattice_size[2] +  \
      (coordinate[2]) * 3 * lattice_size[0] * lattice_size[1] +                \
      (coordinate[1]) * 3 * lattice_size[0] + (coordinate[0]) * 3 + direction

#define PLACE1_LINK_NODIR                                                      \
  (coordinate[3]) * lattice_size[0] * lattice_size[1] * lattice_size[2] +      \
      (coordinate[2]) * lattice_size[0] * lattice_size[1] +                    \
      (coordinate[1]) * lattice_size[0] + (coordinate[0])

link1::link1(int lattice_size_x, int lattice_size_y, int lattice_size_z,
             int lattice_size_t) {
  lattice_size[0] = lattice_size_x;
  lattice_size[1] = lattice_size_y;
  lattice_size[2] = lattice_size_z;
  lattice_size[3] = lattice_size_t;
  direction = 0;
  for (int i = 0; i < 4; i++) {
    coordinate[i] = 0;
    coordinate_old[i] = 0;
  }
  place = 0;
}

ostream &operator<<(ostream &os, const link1 &link) {
  os << "x: " << link.coordinate[0] << " y: " << link.coordinate[1]
     << " z: " << link.coordinate[2] << " t: " << link.coordinate[3]
     << " dir: " << link.direction << endl;
  return os;
}

void link1::move(int dir, int step) {
  coordinate[dir] += (lattice_size[dir] + step);
  coordinate[dir] = coordinate[dir] % lattice_size[dir];
  place += (coordinate[dir] - coordinate_old[dir]) * multiplier[dir];
  coordinate_old[dir] = coordinate[dir];
}

void link1::update(int dir) {
  place += (coordinate[dir] - coordinate_old[dir]) * multiplier[dir];
  coordinate_old[dir] = coordinate[dir];
}

void link1::go(int x, int y, int z, int t) {
  coordinate[0] = x;
  coordinate[1] = y;
  coordinate[2] = z;
  coordinate[3] = t;
}
void link1::move_dir(int dir) { direction = dir; }

template <class T> const T *link1::get_matrix(const vector<T> &array) {
  return &array[place + direction];
}

// PLACE1_LINK_NODIR is defined on the top of the file
// it's a place of site in lexicographical order
// though link has a direction, it's not important here since spin sits on a
// site
const spin *link1::get_spin(const vector<spin> &vec) { return &vec[place / 4]; }

const spin *link1::get_consecutive_spin(const vector<spin> &vec, int mu) {
  int coordinate_new = coordinate[mu] + 1;
  coordinate_new = coordinate_new % lattice_size[mu];
  return &vec[(place + (coordinate_new - coordinate[mu]) * multiplier[mu]) / 4];
}

template <class T> T link1::plaket_mu(const vector<T> &array, int mu) {
  int dir = direction;
  T A = *get_matrix(array);
  move(dir, 1);
  move_dir(mu);
  A = A * get_matrix(array);
  move_dir(dir);
  move(dir, -1);
  move(mu, 1);
  A = A ^ get_matrix(array);
  move_dir(mu);
  move(mu, -1);
  A = A ^ get_matrix(array);
  move_dir(dir);
  return A;
}

template <class T> T link1::polyakov_loop(const vector<T> &array) {
  T A;
  for (int i = 0; i < lattice_size[3]; i++) {
    A = A * get_matrix(array);
    move(direction, 1);
  }
  return A;
}

template <class T> T link1::wilson_loop(const vector<T> &array, int r, int t) {
  int dir = direction;
  T A;
  for (int i = 0; i < r; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  move_dir(3);
  for (int i = 0; i < t; i++) {
    A = A * get_matrix(array);
    move(3, 1);
  }
  move_dir(dir);
  for (int i = 0; i < r; i++) {
    move(dir, -1);
    A = A ^ get_matrix(array);
  }
  move_dir(3);
  for (int i = 0; i < t; i++) {
    move(3, -1);
    A = A ^ get_matrix(array);
  }
  move_dir(dir);
  return A;
}

template <class T> T link1::wilson_line(const vector<T> &array, int length) {
  int dir = direction;
  T A;
  for (int i = 0; i < length; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  return A;
}

// specializations

// su2
template const su2 *link1::get_matrix(const vector<su2> &vec);
template su2 link1::plaket_mu(const vector<su2> &array, int mu);
template su2 link1::polyakov_loop(const vector<su2> &array);
template su2 link1::wilson_loop(const vector<su2> &array, int r, int t);
template su2 link1::wilson_line(const vector<su2> &array, int length);

// abelian
template const abelian *link1::get_matrix(const vector<abelian> &vec);
template abelian link1::plaket_mu(const vector<abelian> &array, int mu);
template abelian link1::polyakov_loop(const vector<abelian> &array);
template abelian link1::wilson_loop(const vector<abelian> &array, int r, int t);
template abelian link1::wilson_line(const vector<abelian> &array, int length);
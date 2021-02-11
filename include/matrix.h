#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include <iostream>

using namespace std;

// su2 matrix in sigma matrices representation (a0 + I * ai * sigma[i])
// where I = sqrt(-1)
class su2 {
public:
  FLOAT a0, a1, a2, a3;
  su2(FLOAT b0, FLOAT b1, FLOAT b2, FLOAT b3);
  su2();

  // calculate trace of the matrix
  FLOAT tr();

  // calculate inverse of the matrix
  su2 inverse();

  // calculate conjugate of the matrix
  su2 conj() const;

  // gets projection onto su2 group
  su2 proj();

  // calculates module of vector in sigma matrices representation
  FLOAT module();

  // left and right multiplication of the matrix on sigma[3] matrix
  su2 sigma3_mult() const;
};

su2 operator+(const su2 &A, const su2 &B);
su2 operator-(const su2 &A, const su2 &B);
su2 operator*(const FLOAT &x, const su2 &A);
su2 operator*(const su2 &A, const FLOAT &x);

// matrix multiplication A * B
su2 operator*(const su2 &A, const su2 &B);

// matrix multiplication A * B
su2 operator*(const su2 &A, const su2 *B);

// matrix multiplication A * B.conj()
su2 operator^(const su2 &A, const su2 *B);

ostream &operator<<(ostream &os, const su2 &A);

// abelian variable (module * exp(i * phi))
class abelian {
public:
  FLOAT r, phi;
  abelian();
  abelian(FLOAT r1, FLOAT phi1);

  // trace
  FLOAT tr();

  // inverse
  abelian inverse();
  // conjugated
  abelian conj() const;
  abelian proj();
  FLOAT module();
};

abelian operator+(const abelian &A, const abelian &B);
abelian operator-(const abelian &A, const abelian &B);
abelian operator*(const FLOAT &x, const abelian &A);
abelian operator*(const abelian &A, const FLOAT &x);

abelian operator*(const abelian &A, const abelian &B);
abelian operator*(const abelian &A, const abelian *B);
abelian operator^(const abelian &A, const abelian *B);

ostream &operator<<(ostream &os, const abelian &A);

// 3D vector realisation for spin model.
class spin {
public:
  FLOAT a1, a2, a3;

  spin(FLOAT a1, FLOAT a2, FLOAT a3);
  spin(su2 U);
  spin();

  FLOAT norm();
  // reflect spin-vector through the V-vector axis
  FLOAT reflect(spin &V);
  // rorate spin-vector to the same direction as V-vector
  FLOAT parallel(spin &V);
  // check identity matrix
  bool IsUnit();
  // contribution from neighbour site
  spin contribution(const su2 *A) const;
  // contribution from neighbour site in negative direction (conjugated matrix)
  spin contribution_conj(const su2 *A) const;
  // calculation of gauge matrix G from spin vector.
  su2 GetGaugeMatrix();
};

spin operator+(const spin &A, const spin &B);
spin operator-(const spin &A, const spin &B);
spin operator*(const FLOAT &x, const spin &A);
spin operator*(const spin &A, const FLOAT &x);
FLOAT operator*(const spin &A, const spin &B);

ostream &operator<<(ostream &os, const spin &A);

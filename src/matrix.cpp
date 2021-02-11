#include "../include/matrix.h"
#include <cmath>

// su2 methods
su2::su2() {
  a0 = 1;
  a1 = 0;
  a2 = 0;
  a3 = 0;
}

su2::su2(FLOAT b0, FLOAT b1, FLOAT b2, FLOAT b3) {
  a0 = b0;
  a1 = b1;
  a2 = b2;
  a3 = b3;
}

FLOAT su2::tr() { return 2 * a0; }

su2 su2::inverse() {
  FLOAT rho = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
  return su2(a0 / rho, -a1 / rho, -a2 / rho, -a3 / rho);
}
FLOAT su2::module() { return a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3; }
su2 su2::conj() const { return su2(a0, -a1, -a2, -a3); }
su2 su2::proj() {
  FLOAT rho = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
  return su2(a0 / powf(rho, 0.5), a1 / powf(rho, 0.5), a2 / powf(rho, 0.5),
             a3 / powf(rho, 0.5));
}
su2 su2::sigma3_mult() const { return su2(a0, -a1, -a2, a3); }

su2 operator+(const su2 &A, const su2 &B) {
  return su2(A.a0 + B.a0, A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
};
su2 operator-(const su2 &A, const su2 &B) {
  return su2(A.a0 - B.a0, A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
};
su2 operator*(const FLOAT &x, const su2 &A) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
};
su2 operator*(const su2 &A, const FLOAT &x) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
};

su2 operator*(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 - A.a1 * B.a1 - A.a2 * B.a2 - A.a3 * B.a3,
             A.a0 * B.a1 + B.a0 * A.a1 + A.a3 * B.a2 - A.a2 * B.a3,
             A.a0 * B.a2 + B.a0 * A.a2 + A.a1 * B.a3 - A.a3 * B.a1,
             A.a0 * B.a3 + B.a0 * A.a3 + A.a2 * B.a1 - A.a1 * B.a2);
};
su2 operator*(const su2 &A, const su2 *B) {
  return su2(A.a0 * B->a0 - A.a1 * B->a1 - A.a2 * B->a2 - A.a3 * B->a3,
             A.a0 * B->a1 + B->a0 * A.a1 + A.a3 * B->a2 - A.a2 * B->a3,
             A.a0 * B->a2 + B->a0 * A.a2 + A.a1 * B->a3 - A.a3 * B->a1,
             A.a0 * B->a3 + B->a0 * A.a3 + A.a2 * B->a1 - A.a1 * B->a2);
};
su2 operator^(const su2 &A, const su2 *B) {
  return su2(A.a0 * B->a0 + A.a1 * B->a1 + A.a2 * B->a2 + A.a3 * B->a3,
             -A.a0 * B->a1 + B->a0 * A.a1 - A.a3 * B->a2 + A.a2 * B->a3,
             -A.a0 * B->a2 + B->a0 * A.a2 - A.a1 * B->a3 + A.a3 * B->a1,
             -A.a0 * B->a3 + B->a0 * A.a3 - A.a2 * B->a1 + A.a1 * B->a2);
};

ostream &operator<<(ostream &os, const su2 &A) {
  os << "a0 = " << A.a0 << " "
     << "a1 = " << A.a1 << " "
     << "a2 = " << A.a2 << " "
     << "a3 = " << A.a3;
  return os;
}

// abelian methods
abelian::abelian() {
  r = 1;
  phi = 0;
}

abelian::abelian(FLOAT r1, FLOAT phi1) {
  r = r1;
  phi = phi1;
}

FLOAT abelian::tr() { return cos(phi); }

abelian abelian::inverse() { return abelian(1 / r, -phi); }
FLOAT abelian::module() { return r; }
abelian abelian::conj() const { return abelian(r, -phi); }
abelian abelian::proj() { return abelian(1, phi); }

abelian operator+(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) + B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) + B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) + B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) + B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) + B.r * sin(B.phi),
                       A.r * cos(A.phi) + B.r * cos(B.phi)));
};
abelian operator-(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) - B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) - B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) - B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) - B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) - B.r * sin(B.phi),
                       A.r * cos(A.phi) - B.r * cos(B.phi)));
};
abelian operator*(const FLOAT &x, const abelian &A) {
  return abelian(A.r * x, A.phi);
};
abelian operator*(const abelian &A, const FLOAT &x) {
  return abelian(A.r * x, A.phi);
};

abelian operator*(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, A.phi + B.phi);
};
abelian operator*(const abelian &A, const abelian *B) {
  return abelian(A.r * B->r, A.phi + B->phi);
};
abelian operator^(const abelian &A, const abelian *B) {
  return abelian(A.r * B->r, A.phi - B->phi);
};

ostream &operator<<(ostream &os, const abelian &A) {
  os << "r = " << A.r << " "
     << "phi = " << A.phi << " " << endl;
  return os;
}

spin::spin() {
  a1 = (FLOAT)0.0;
  a2 = (FLOAT)0.0;
  a3 = (FLOAT)1.0;
}
spin::spin(su2 U) {
  a1 = 2 * (U.a0 * U.a2 + U.a3 * U.a1);
  a2 = 2 * (U.a3 * U.a2 + U.a0 * U.a1);
  a3 = 1. - 2 * (U.a2 * U.a2 + U.a1 * U.a1);
}
spin::spin(FLOAT a1, FLOAT a2, FLOAT a3) : a1(a1), a2(a2), a3(a3) {}

FLOAT spin::norm() { return sqrt(a1 * a1 + a2 * a2 + a3 * a3); }

// Function reflect spin vector through arbitrary spin vector
// it's more convenient to just change the spin variable
FLOAT spin::reflect(spin &V) {
  FLOAT tmp1 = (V.a1 * V.a1 - V.a2 * V.a2 - V.a3 * V.a3) * a1 +
               2 * V.a1 * (V.a2 * a2 + V.a3 * a3);

  FLOAT tmp2 = (-V.a1 * V.a1 + V.a2 * V.a2 - V.a3 * V.a3) * a2 +
               2 * V.a2 * (V.a1 * a1 + V.a3 * a3);

  FLOAT tmp3 = (-V.a1 * V.a1 - V.a2 * V.a2 + V.a3 * V.a3) * a3 +
               2 * V.a3 * (V.a1 * a1 + V.a2 * a2);

  FLOAT norm2 = V.norm() * V.norm();

  FLOAT b1 = tmp1 / norm2;
  FLOAT b2 = tmp2 / norm2;
  FLOAT b3 = tmp3 / norm2;

  FLOAT c = b1 * a1 + b2 * a2 + b3 * a3;

  a1 = tmp1 / norm2;
  a2 = tmp2 / norm2;
  a3 = tmp3 / norm2;

  return c;
}

FLOAT spin::parallel(spin &V) {
  FLOAT vNorm = V.norm();

  FLOAT b1 = V.a1 / vNorm;
  FLOAT b2 = V.a2 / vNorm;
  FLOAT b3 = V.a3 / vNorm;

  FLOAT c = b1 * a1 + b2 * a2 + b3 * a3;

  a1 = b1;
  a2 = b2;
  a3 = b3;

  return c;
}

bool spin::IsUnit() { return (this->norm() - 1.) > 1e-6 ? false : true; }

// vector which spin variable is multiplied on
// it's contribution from single neighbour site in one direction
spin spin::contribution(const su2 *A) const {
  FLOAT q1 = A->a1 * A->a1;
  FLOAT q2 = A->a2 * A->a2;
  FLOAT q3 = A->a3 * A->a3;
  FLOAT q12 = 2. * A->a1 * A->a2;
  FLOAT q13 = 2. * A->a1 * A->a3;
  FLOAT q23 = 2. * A->a2 * A->a3;
  FLOAT q01 = 2. * A->a0 * A->a1;
  FLOAT q02 = 2. * A->a0 * A->a2;
  FLOAT q03 = 2. * A->a0 * A->a3;
  FLOAT A5 = A->a0 * A->a0 - q1 - q2 - q3;

  // matrix A(i, j) = tr(sigma(i) * U * sigma(j) * U.conj) is multiplied by
  // spin b variable from the right (A * b)
  return spin(a1 * (A5 + 2 * q1) + a2 * (q12 + q03) + a3 * (q13 - q02),
              a1 * (q12 - q03) + a2 * (A5 + 2 * q2) + a3 * (q23 + q01),
              a1 * (q13 + q02) + a2 * (q23 - q01) + a3 * (A5 + 2 * q3));
}

// supposed to be used for spin variable on current lattice cite
spin spin::contribution_conj(const su2 *A) const {
  FLOAT q1 = A->a1 * A->a1;
  FLOAT q2 = A->a2 * A->a2;
  FLOAT q3 = A->a3 * A->a3;
  FLOAT q12 = 2. * A->a1 * A->a2;
  FLOAT q13 = 2. * A->a1 * A->a3;
  FLOAT q23 = 2. * A->a2 * A->a3;
  FLOAT q01 = 2. * A->a0 * A->a1;
  FLOAT q02 = 2. * A->a0 * A->a2;
  FLOAT q03 = 2. * A->a0 * A->a3;
  FLOAT A5 = A->a0 * A->a0 - q1 - q2 - q3;

  // matrix A(i, j) = tr(sigma(i) * U * sigma(j) * U.conj) is multiplied by
  //  spin b variable from the left (b * A)
  return spin(a1 * (A5 + 2 * q1) + a2 * (q12 - q03) + a3 * (q13 + q02),
              a1 * (q12 + q03) + a2 * (A5 + 2 * q2) + a3 * (q23 - q01),
              a1 * (q13 - q02) + a2 * (q23 + q01) + a3 * (A5 + 2 * q3));
}

// Calculation have been done for specific element of gauge matrix G with g_3 =
// 0 Here:
//      | g_0 + I g_1    g_2 + I g_3 |
// G =  |                            |
//      |-g_2 + I g_3    g_0 - I g_1 |
//
// spinVector = (a, b, c) = (a1, a2, a3)
//
// g_0 = - a/sqrt(2-2*c)
// g_1 = - b/sqrt(2-2*c)
// g_2 = - sqrt(2-2*c) / 2
// g_3 = 0
su2 spin::GetGaugeMatrix() {
  FLOAT sq = sqrt(2. - 2 * a3);
  return su2(-a1 / sq, 0., -sq / 2., -a2 / sq);
}
spin operator*(const FLOAT &x, const spin &A) {
  return spin(A.a1 * x, A.a2 * x, A.a3 * x);
}
spin operator*(const spin &A, const FLOAT &x) {
  return spin(A.a1 * x, A.a2 * x, A.a3 * x);
}
spin operator+(const spin &A, const spin &B) {
  return spin(A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
}
spin operator-(const spin &A, const spin &B) {
  return spin(A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
}
FLOAT operator*(const spin &A, const spin &B) {
  return A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3;
}

ostream &operator<<(ostream &os, const spin &A) {
  os << "a1 = " << A.a1 << " "
     << "a2 = " << A.a2 << " "
     << "a3 = " << A.a3 << " ";
  return os;
}
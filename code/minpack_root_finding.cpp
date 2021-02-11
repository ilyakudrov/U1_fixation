#include "cminpack.h"
#include "math.h"
#include <chrono>
#include <iostream>
#include <vector>

using namespace std;

vector<double> theta(24);

int function_test(void *p, int n, const double *x, double *fvec, int iflag) {
  // positive directions
  fvec[0] = 0;
  fvec[1] = 0;

  for (int i = 0; i < 4; i++) {
    fvec[0] += sin(x[0] - theta[i * 3]) + sin(x[0] + x[1] + theta[i * 3 + 2]) +
               sin(x[0] + theta[(i + 4) * 3]) +
               sin(x[0] + x[1] - theta[(i + 4) * 3 + 2]);

    fvec[1] += sin(x[1] - theta[i * 3 + 1]) +
               sin(x[0] + x[1] + theta[i * 3 + 2]) +
               sin(x[1] + theta[(i + 4) * 3 + 1]) +
               sin(x[0] + x[1] - theta[(i + 4) * 3 + 2]);
  }

  return 0;
}

int main() {

  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  theta = {-0.003430567378423177, 2.5187546645278616,   -0.18790892227498368,
           1.8580870595325507,    0.5462893265144331,   -0.11351799849439281,
           0.5752081963953266,    2.161798709021106,    -0.10056541433062449,
           2.4770558204233915,    0.014931896904666786, 2.1720296730661524,
           -1.2834853852652746,   0.7885560009794599,   1.1986857960623745,
           2.7952004417190537,    2.930628506228693,    0.48003331213481815,
           0.527335578519212,     2.297331991870097,    2.3595674142926963,
           -3.109486430507647,    -2.4471192629700935,  -1.1906340553863846};

  int j, n, info, lwa;
  double tol, fnorm;
  double x[2], fvec[2], wa[19];

  n = 2;

  /*      the following starting values provide a rough solution. */

  x[0] = 0.123;
  x[1] = 0.4363;

  lwa = 19;

  tol = sqrt(__cminpack_func__(dpmpar)(1));

  start_time = clock();
  for (int i = 0; i < 1000000; i++) {
    info =
        __cminpack_func__(hybrd1)(function_test, 0, n, x, fvec, tol, wa, lwa);
  }

  end_time = clock();
  search_time = end_time - start_time;
  cout << "roots time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
  fnorm = __cminpack_func__(enorm)(n, fvec);

  cout << "     final L2 norm of the residuals " << (double)fnorm << endl;
  cout << "     exit parameter " << info << endl;
  cout << "     final approximates solution " << endl;

  for (j = 1; j <= n; j++) {
    printf("%s%15.7g", j % 3 == 1 ? "\n     " : "", (double)x[j - 1]);
  }
  printf("\n");
}
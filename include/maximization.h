#pragma once

#include "../include/link.h"
#include "../include/roots.h"
#include <math.h>

struct angles {
  double angle[3];
};

struct gauge {
  double alpha[2];
};

void claculate_params(vector<angles> &angles, vector<gauge> &gauge, link1 &link,
                      functional_params &params);

double functional_test(vector<angles> &angles, vector<gauge> &gauge);

double functional(vector<angles> &angles, vector<gauge> &gauge);

struct maximize_data {
  int iteration;
  double functional;
};

vector<maximize_data> maximize(vector<angles> &angles, vector<gauge> &gauge);

void maximize_test(vector<angles> &angles, vector<gauge> &gauge);
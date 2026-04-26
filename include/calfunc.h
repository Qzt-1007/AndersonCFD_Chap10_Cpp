#pragma once
#include <constant.h>

#include <cmath>

double Vis(double T) {
    double pTa = T / T0;
    double pTb = (T0 + 110.0) / (T + 110.0);
    return Vis0 * pow(pTa, 1.50) * pTb;
}

double Therm_k(double Vis) { return Vis * Cp / Pr; }
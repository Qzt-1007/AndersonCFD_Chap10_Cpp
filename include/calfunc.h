#pragma once
#include <constant.h>

#include <cmath>

double SutVis(double T) {
    double pTa = T / T0;
    double pTb = (T0 + 110.0) / (T + 110.0);
    return Vis0 * pow(pTa, 1.50) * pTb;
}

double Therm_k(double SutVis) { return SutVis * Cp / Pr; }

double Rho(double p, double T) { return p / (R * T); }

double Sonic(double T) { return sqrt(gamma * R * T); }
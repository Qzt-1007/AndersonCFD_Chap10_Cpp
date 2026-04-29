#pragma once
#include <FArray.h>
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

// 根据U反解出原始量(Q:u,v,p,t)
Array3D Calc_Prim(Array3D U) {
    Array3D Q = createArray3D(IMAX, JMAX, 4, 0.0);  // Q:u,v,p,t
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            Q[i][j][0] = U[i][j][1] / U[i][j][0];
            Q[i][j][1] = U[i][j][2] / U[i][j][0];
            Q[i][j][3] = (U[i][j][3] / U[i][j][0] - 0.50 * (Q[i][j][0] * Q[i][j][0] + Q[i][j][1] * Q[i][j][1])) / Cv;
            Q[i][j][2] = U[i][j][0] * R * Q[i][j][3];
        }
    }
    return Q;
}

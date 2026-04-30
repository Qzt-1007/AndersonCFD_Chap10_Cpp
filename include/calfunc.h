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

double Sonic(double T) { return sqrt(gamma * R * T); }

// 根据U反解出原始量(Q:u,v,p,t)
Array3D Calc_Prim(Array3D U) {
    Array3D Q = createArray3D(IMAX, JMAX, 4, 0.0);  // Q:u,v,p,t
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            Q[i][j][0] = U[i][j][1] / U[i][j][0];
            Q[i][j][1] = U[i][j][2] / U[i][j][0];
            Q[i][j][3] = ((U[i][j][3] / U[i][j][0]) - 0.50 * (Q[i][j][0] * Q[i][j][0] + Q[i][j][1] * Q[i][j][1])) / Cv;
            Q[i][j][2] = U[i][j][0] * R * Q[i][j][3];
        }
    }
    return Q;
}

// 计算时间步长
double Tstep(Array2D u, Array2D v, Array2D p, Array2D T) {
    double vd = 0.0;
    double deltaT;
    double deltaTmin = 1e9;
    for (int i = 1; i < IMAX - 1; i++) {
        for (int j = 1; j < JMAX - 1; j++) {
            vd = (4.0 * gamma * SutVis(T[i][j])) / (3 * Pr * p[i][j] / R / T[i][j]);
            deltaT = 1.0 / (abs(u[i][j]) / Dx + abs(v[i][j]) / Dy + Sonic(T[i][j]) * sqrt(1.0 / Dx / Dx + 1.0 / Dy / Dy) +
                          2.0 * vd * (1.0 / Dx / Dx + 1 / Dy / Dy));
            if (deltaT < deltaTmin) deltaTmin = deltaT;
        }
    }
    return Cfl * deltaT;
}

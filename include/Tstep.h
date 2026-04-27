#pragma once
#include <FArray.h>
#include <calfunc.h>
#include <constant.h>

#include <cmath>
#include <vector>

// 计算时间步长
double Tstep(Array2D u,Array2D v,Array2D p,Array2D T) {  // Q:u,v,p,T
    double vdmax = (4 * gamma * SutVis(T[0][0])) / (3 * Pr * Rho(T[0][0],p[0][0]));
    double vd;
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            vd = (4 * gamma * SutVis(T[i][j])) / (3 * Pr * Rho(T[i][j],p[i][j]));
            if (vd > vdmax) vdmax = vd;
        }
    }

    double deltaT = 1 / (abs(u[1][1]) / Dx + abs(v[1][1]) / Dy +
                         Sonic(T[1][1]) * sqrt(1 / Dx / Dx + 1 / Dy / Dy) + 2 * vdmax * (1 / Dx / Dx + 1 / Dy / Dy));
    double deltaTmin;
    for (int i = 1; i < IMAX; i++) {
        for (int j = 1; j < JMAX; j++) {
            deltaT =
                1 / (abs(u[i][j]) / Dx + abs(v[i][j]) / Dy + Sonic(T[i][j]) * sqrt(1 / Dx / Dx + 1 / Dy / Dy) +
                     2 * vdmax * (1 / Dx / Dx + 1 / Dy / Dy));
            if (deltaT < deltaTmin) deltaTmin = deltaT;
        }
    }
    return Cfl * deltaT;
}
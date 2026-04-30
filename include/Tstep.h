#pragma once
#include <FArray.h>
#include <calfunc.h>
#include <constant.h>

#include <cmath>
#include <vector>

// 计算时间步长
double Tstep(Array2D u, Array2D v, Array2D p, Array2D T) {
    double vd = 0.0;
    double deltaT;
    double deltaTmin = 1e9;
    for (int i = 1; i < IMAX - 1; i++) {
        for (int j = 1; j < JMAX - 1; j++) {
            vd = (4 * gamma * SutVis(T[i][j])) / (3 * Pr * p[i][j] / R / T[i][j]);
            deltaT = 1 / (abs(u[i][j]) / Dx + abs(v[i][j]) / Dy + Sonic(T[i][j]) * sqrt(1 / Dx / Dx + 1 / Dy / Dy) +
                          2 * vd * (1 / Dx / Dx + 1 / Dy / Dy));
            if (deltaT < deltaTmin) deltaTmin = deltaT;
        }
    }
    return Cfl * deltaT;
}
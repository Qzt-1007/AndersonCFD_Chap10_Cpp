#pragma once
#include <calfunc.h>
#include <constant.h>

#include <cmath>
#include <vector>

// 计算时间步长
double Tstep(std::vector<std::vector<std::vector<double>>> Q) {
    double vdmax = (4 * gamma * Vis(Q[0][0][3])) / (3 * Pr * Rho(Q[0][0][2], Q[0][0][3]));
    double vd;
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            vd = (4 * gamma * Vis(Q[i][j][3])) / (3 * Pr * Rho(Q[i][j][2], Q[i][j][3]));
            if (vd > vdmax) vdmax = vd;
        }
    }
    double deltaT = 1 / (abs(Q[1][1][0]) / Dx + abs(Q[1][1][1]) / Dy +
                         Sonic(Q[1][1][3]) * sqrt(1 / Dx / Dx + 1 / Dy / Dy) + 2 * vdmax * (1 / Dx / Dx + 1 / Dy / Dy));
    double deltaTmin;
    for (int i = 1; i < IMAX; i++) {
        for (int j = 1; j < JMAX; j++) {
            deltaT =
                1 / (abs(Q[i][j][0]) / Dx + abs(Q[i][j][1]) / Dy + Sonic(Q[i][j][3]) * sqrt(1 / Dx / Dx + 1 / Dy / Dy) +
                     2 * vdmax * (1 / Dx / Dx + 1 / Dy / Dy));
            if (deltaT < deltaTmin) deltaTmin = deltaT;
        }
    }
    return Cfl * deltaT;
}
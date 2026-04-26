#pragma once
#include <calfunc.h>
#include <constant.h>

#include <vector>
#include <cmath>

// 计算时间步长
double Tstep(std::vector<std::vector<std::vector<double>>> P) {
    double vdmax = (4 * gamma * Vis(P[0][0][3])) / (3 * Pr * Rho(P[0][0][2], P[0][0][3]));
    
}
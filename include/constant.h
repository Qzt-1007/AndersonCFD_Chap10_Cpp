#pragma once
#include <cmath>
// 用于定义常量

const int IMAX = 70;
const int JMAX = 70;
const int MaxIter = 8000;

const double Ma0 = 4.0;
const double T0 = 288.16;  // 海平面温度，也是来流温度
const double a0 = 340.28;
const double p0 = 101325.0;
const double gamma = 1.4;      // 绝热比
const double R = 287.0;        // 气体常数
const double LHORI = 0.00001;  // 平板长度

const double Tw = T0;
const double Pr = 0.71;
const double Vis0 = 1.789e-5;  // 海平面处空气粘度

const double Cv = R / (gamma - 1);
const double Cp = gamma * Cv;
const double Rho0 = p0 / (R * T0);
const double ReL = Rho0 * (Ma0 * a0) * LHORI / Vis0;
const double e0 = Cv * T0;
const double k0 = Vis0 * Cp / Pr;
const double Dd = 5 * LHORI / sqrt(ReL);
const double LVERT = 5 * Dd;
const double Dx = LHORI / (IMAX - 1);
const double Dy = LVERT / (JMAX - 1);
const double Cfl = 0.60;  // 柯朗数
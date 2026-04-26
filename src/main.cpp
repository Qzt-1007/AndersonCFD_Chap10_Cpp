#include <cmath>
#include <iostream>

constexpr int IMAX = 70;
constexpr int JMAX = 70;
constexpr int MaxIter = 5000;

constexpr double Ma0 = 4.0;
constexpr double T0 = 288.16;  // 海平面温度，也是来流温度
constexpr double a0 = 340.28;
constexpr double p0 = 101325.0;
constexpr double gamma = 1.4;      // 绝热比
constexpr double R = 287.0;        // 气体常数
constexpr double LHORI = 0.00001;  // 平板长度

constexpr double Tw = T0;
constexpr double Pr = 0.71;
constexpr double DynVis0 = 1.789e-5;  // 海平面处空气粘度

constexpr double Cv = R / (gamma - 1);
constexpr double Cp = gamma * Cv;
constexpr double Rho0 = p0 / (R * T0);
constexpr double ReL = Rho0 * (Ma0 * a0) * LHORI / DynVis0;
constexpr double e0 = Cv * T0;
constexpr double k0 = DynVis0 * Cp / Pr;
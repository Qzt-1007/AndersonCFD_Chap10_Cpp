#include <cmath>
#include <iostream>

constexpr int IMAX = 70;
constexpr int JMAX = 70;
constexpr int MaxIter = 5000;

constexpr double Ma0 = 4.0;
constexpr double T0 = 288.16;
constexpr double a0 = 340.28;
constexpr double p0 = 101325.0;
constexpr double gamma = 1.4;      // 绝热比
constexpr double R = 287.0;        // 气体常数
constexpr double LHORI = 0.00001;  // 平板长度

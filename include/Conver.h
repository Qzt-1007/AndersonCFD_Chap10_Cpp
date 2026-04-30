#include <FArray.h>
#include <constant.h>

#include <cmath>

constexpr double eps_p = 1e-8;

bool Conver(Array2D rho1, Array2D rho2) {
    double drho = 1e-9;
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            if (abs(rho1[i][j] - rho2[i][j]) > drho) drho = abs(rho1[i][j] - rho2[i][j]);
        }
    }
    if (drho < eps_p)
        return true;
    else
        return false;
}
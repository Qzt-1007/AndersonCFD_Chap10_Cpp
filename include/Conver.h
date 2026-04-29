#include <FArray.h>
#include <constant.h>

#include <cmath>

constexpr double eps_p = 1e-6;

bool Conver(Array2D p1, Array2D p2) {
    double dp = 0.0;
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            if (abs(p1[i][j] - p2[i][j]) > dp) dp = abs(p1[i][j] - p2[i][j]);
        }
    }
    if (dp < eps_p)
        return true;
    else
        return false;
}
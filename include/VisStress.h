#pragma once
#include <FArray.h>
#include <constant.h>

#include <stdexcept>
#include <vector>

Array2D Tau_xx(Array2D u, Array2D v, Array2D Vis, int indi)  // indi=0--前向差分,indi=1--后向差分
{
    Array2D Tau, dudx, dvdy = createArray2D(IMAX, JMAX, 0.0);
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            if (j == 0) {
                dvdy[i][j] = (-3 * v[i][j] + 4 * v[i][j + 1] - v[i][j + 2]) / (2 * Dy);
            } else if (j == JMAX - 1) {
                dvdy[i][j] = (-3 * v[i][j] + 4 * v[i][j - 1] - v[i][j - 2]) / (2 * Dy);
            } else {
                dvdy[i][j] = (v[i][j + 1] - v[i][j - 1]) / (2 * Dy);
            }
        }
    }

    if (indi == 0) {  // 前向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == IMAX - 1) {
                    dudx[i][j] = (-3 * u[i][j] + 4 * u[i - 1][j] - u[i - 2][j]) / (2 * Dx);
                } else {
                    dudx[i][j] = (u[i + 1][j] - u[i][j]) / Dx;
                }
            }
        }
    } else if (indi == 1) {  // 后向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == 0) {
                    dudx[i][j] = (-3 * u[i][j] + 4 * u[i + 1][j] - u[i + 2][j]) / (2 * Dx);
                } else {
                    dudx[i][j] = (u[i][j] - u[i - 1][j]) / Dx;
                }
            }
        }
    } else {
        throw std::runtime_error("indi should only be 0 or 1.");
    }

    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            Tau[i][j] = (4.0 / 3.0) * Vis[i][j] * dudx[i][j] - (2.0 / 3.0) * Vis[i][j] * dvdy[i][j];
        }
    }

    return Tau;
}

Array2D Tau_yy(Array2D u, Array2D v, Array2D Vis, int indi)  // indi=0--前向差分,indi=1--后向差分
{
    Array2D Tau, dudx, dvdy = createArray2D(IMAX, JMAX, 0.0);
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            if (j == 0) {
                dudx[i][j] = (-3 * u[i][j] + 4 * u[i][j + 1] - u[i][j + 2]) / (2 * Dx);
            } else if (j == JMAX - 1) {
                dudx[i][j] = (-3 * u[i][j] + 4 * u[i][j - 1] - u[i][j - 2]) / (2 * Dx);
            } else {
                dudx[i][j] = (u[i][j + 1] - u[i][j - 1]) / (2 * Dx);
            }
        }
    }

    if (indi == 0) {  // 前向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == IMAX - 1) {
                    dvdy[i][j] = (-3 * v[i][j] + 4 * v[i - 1][j] - v[i - 2][j]) / (2 * Dy);
                } else {
                    dvdy[i][j] = (v[i + 1][j] - v[i][j]) / Dx;
                }
            }
        }
    } else if (indi == 1) {  // 后向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == 0) {
                    dvdy[i][j] = (-3 * v[i][j] + 4 * v[i + 1][j] - v[i + 2][j]) / (2 * Dy);
                } else {
                    dvdy[i][j] = (v[i][j] - v[i - 1][j]) / Dy;
                }
            }
        }
    } else {
        throw std::runtime_error("indi should only be 0 or 1.");
    }

    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            Tau[i][j] = (4.0 / 3.0) * Vis[i][j] * dvdy[i][j] - (2.0 / 3.0) * Vis[i][j] * dudx[i][j];
        }
    }

    return Tau;
}

Array2D Tau_xy_E(Array2D u, Array2D v, Array2D Vis, int indi) {  // 计算E项中的tau_xy,E在方程中以pe/px形式存在
    Array2D Tau, dudy, dvdx = createArray2D(IMAX, JMAX, 0.0);
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            if (j == 0) {
                dudy[i][j] = (-3 * u[i][j] + 4 * u[i][j + 1] - u[i][j + 2]) / (2 * Dy);
            } else if (j == JMAX - 1) {
                dudy[i][j] = (-3 * u[i][j] + 4 * u[i][j - 1] - u[i][j - 2]) / (2 * Dy);
            } else {
                dudy[i][j] = (u[i][j + 1] - u[i][j - 1]) / (2 * Dy);
            }
        }
    }
    if (indi == 0) {  // 前向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == IMAX - 1) {
                    dvdx[i][j] = (-3 * v[i][j] + 4 * v[i - 1][j] - v[i - 2][j]) / (2 * Dx);
                } else {
                    dvdx[i][j] = (v[i + 1][j] - v[i][j]) / Dx;
                }
            }
        }
    } else if (indi == 1) {  // 后向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == 0) {
                    dvdx[i][j] = (-3 * v[i][j] + 4 * v[i + 1][j] - v[i + 2][j]) / (2 * Dx);
                } else {
                    dvdx[i][j] = (v[i][j] - v[i - 1][j]) / Dx;
                }
            }
        }
    } else {
        throw std::runtime_error("indi should only be 0 or 1.");
    }

    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            Tau[i][j] = Vis[i][j] * (dudy[i][j] + dvdx[i][j]);
        }
    }
    return Tau;
}

Array2D Tau_xy_F(Array2D u, Array2D v, Array2D Vis, int indi) {  // 计算F项中的tau_xy,E在方程中以pF/py形式存在
    Array2D Tau, dudy, dvdx = createArray2D(IMAX, JMAX, 0.0);
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            if (j == 0) {
                dvdx[i][j] = (-3 * v[i][j] + 4 * v[i][j + 1] - v[i][j + 2]) / (2 * Dx);
            } else if (j == JMAX - 1) {
                dvdx[i][j] = (-3 * v[i][j] + 4 * v[i][j - 1] - v[i][j - 2]) / (2 * Dx);
            } else {
                dvdx[i][j] = (v[i][j + 1] - v[i][j - 1]) / (2 * Dx);
            }
        }
    }
    if (indi == 0) {  // 前向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == IMAX - 1) {
                    dudy[i][j] = (-3 * u[i][j] + 4 * u[i - 1][j] - u[i - 2][j]) / (2 * Dy);
                } else {
                    dudy[i][j] = (u[i + 1][j] - u[i][j]) / Dy;
                }
            }
        }
    } else if (indi == 1) {  // 后向差分
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == 0) {
                    dudy[i][j] = (-3 * u[i][j] + 4 * u[i + 1][j] - u[i + 2][j]) / (2 * Dy);
                } else {
                    dudy[i][j] = (u[i][j] - u[i - 1][j]) / Dy;
                }
            }
        }
    } else {
        throw std::runtime_error("indi should only be 0 or 1.");
    }

    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            Tau[i][j] = Vis[i][j] * (dudy[i][j] + dvdx[i][j]);
        }
    }
    return Tau;
}
/*麦考马克方法*/
#pragma once
#include <FArray.h>
#include <Tstep.h>
#include <VisStress.h>
#include <calfunc.h>
#include <constant.h>

#include <vector>

// 计算热流量的函数
Array2D q_x(Array2D, Array2D, int);
Array2D q_y(Array2D, Array2D, int);

Array3D Maccormack(Array3D Q)  // Q:u,v,p,T
{
    Array2D u, v, p, T, rho, e, Vis, k = createArray2D(IMAX, JMAX, 0.0);
    Array3D U, E, F = createArray3D(IMAX, JMAX, 4);  // U,E,F有4个分量
    Array2D Tauxx, Tauyy, Tauxy = createArray2D(IMAX, JMAX, 0.0);
    Array2D qx, qy = createArray2D(IMAX, JMAX, 0.0);
    Array3D pU, pE, pF = createArray3D(IMAX, JMAX, 4, 0.0);
    double Dt = 0.0;

    // 初始化
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            u[i][j] = Q[i][j][0];
            v[i][j] = Q[i][j][1];
            p[i][j] = Q[i][j][2];
            T[i][j] = Q[i][j][3];
            rho[i][j] = p[i][j] / (R * T[i][j]);
            e[i][j] = Cv * T[i][j];
            Vis[i][j] = SutVis(T[i][j]);
            k[i][j] = (Vis[i][j] * Cp) / Pr;
        }
    }

    // 预估步准备(前向差分)
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {                                                         // 计算每个网格点上的U
            U[i][j][0] = rho[i][j];                                                              // U1
            U[i][j][1] = rho[i][j] * u[i][j];                                                    // U2
            U[i][j][2] = rho[i][j] * v[i][j];                                                    // U3
            U[i][j][3] = rho[i][j] * (e[i][j] + (u[i][j] * u[i][j] + v[i][j] * v[i][j]) / 2.0);  // U5
        }
    }
    // 为计算E,计算每个点上的切应力和热流量
    Tauxx = Tau_xx(u, v, Vis, 1);  // 预估步前向差分,对应的切应力后向差分
    Tauxy = Tau_xy_E(u, v, Vis, 1);
    qx = q_x(T, k, 1);
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {                                         // 计算每个网格点上的E
            E[i][j][0] = rho[i][j] * u[i][j];                                    // E1
            E[i][j][1] = rho[i][j] * u[i][j] * u[i][j] + p[i][j] - Tauxx[i][j];  // E2
            E[i][j][2] = rho[i][j] * u[i][j] * v[i][j] - Tauxy[i][j];            // E3
            E[i][j][3] = (U[i][j][3] + p[i][j]) - u[i][j] * Tauxx[i][j] - v[i][j] * Tauxy[i][j] + qx[i][j];  // E5
        }
    }
    // 为计算F,计算每个点上的切应力和热流量
    Tauyy = Tau_yy(u, v, Vis, 1);
    Tauxy = Tau_xy_F(u, v, Vis, 1);
    qy = q_y(T, k, 1);
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {  // 计算每个网格点上的F
            F[i][j][0] = rho[i][j] * v[i][j];
            F[i][j][1] = rho[i][j] * u[i][j] * v[i][j] - Tauxy[i][j];
            F[i][j][2] = rho[i][j] * v[i][j] * v[i][j] + p[i][j] - Tauyy[i][j];
            F[i][j][3] = (U[i][j][3] + p[i][j]) - v[i][j] * Tauyy[i][j] - u[i][j] * Tauxy[i][j] + qy[i][j];
        }
    }
    // 预估步-前向差分(只计算内部网格点)
    Dt = Tstep(u, v, p, T);
    for (int i = 0; i < IMAX - 1; i++) {
        for (int j = 0; j < JMAX - 1; j++) {
            for (int k = 0; k < 4; k++) {
                pU[i][j][k] =
                    U[i][j][k] - (Dt / Dx) * (E[i + 1][j][k] - E[i][j][k]) - (Dt / Dy) * (F[i][j + 1][k] - F[i][j][k]);
            }
        }
    }
}

Array2D q_x(Array2D T, Array2D k, int indi = 0) {
    Array2D q = createArray2D(IMAX, JMAX, 0.0);
    if (indi == 0) {
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == IMAX - 1) {
                    q[i][j] = k[i][j] * (-3 * T[i][j] + 4 * T[i - 1][j] - T[i - 2][j]) / (2 * Dx);
                } else {
                    q[i][j] = k[i][j] * (T[i + 1][j] - T[i][j]) / Dx;
                }
            }
        }
    } else if (indi == 1) {
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (i == 0) {
                    q[i][j] = k[i][j] * (-3 * T[i][j] + 4 * T[i + 1][j] - T[i + 2][j]) / (2 * Dx);
                } else {
                    q[i][j] = k[i][j] * (T[i][j] - T[i - 1][j]) / Dx;
                }
            }
        }
    } else {
        throw std::runtime_error("indi should only be 0 or 1.");
    }
    return q;
}

Array2D q_y(Array2D T, Array2D k, int indi = 0) {
    Array2D q = createArray2D(IMAX, JMAX, 0.0);
    if (indi == 0) {
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (j == JMAX - 1) {
                    q[i][j] = k[i][j] * (-3 * T[i][j] + 4 * T[i][j - 1] - T[i][j - 2]) / (2 * Dy);
                } else {
                    q[i][j] = k[i][j] * (T[i][j + 1] - T[i][j]) / Dy;
                }
            }
        }
    } else if (indi == 1) {
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                if (j == 0) {
                    q[i][j] = k[i][j] * (-3 * T[i][j] + 4 * T[i][j + 1] - T[i][j + 2]) / (2 * Dy);
                } else {
                    q[i][j] = k[i][j] * (T[i][j] - T[i][j - 1]) / Dy;
                }
            }
        }
    } else {
        throw std::runtime_error("indi should only be 0 or 1.");
    }
    return q;
}
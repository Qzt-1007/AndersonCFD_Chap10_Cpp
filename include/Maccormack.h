/*麦考马克方法*/
#pragma once
#include <FArray.h>
#include <calfunc.h>
#include <constant.h>

#include <vector>

Array3D Maccormack(Array3D Q)  // Q:u,v,p,T
{
    Array2D u, v, p, T, rho, e, Vis, k = createArray2D(IMAX, JMAX, 0.0);
    Array3D U, E, F = createArray3D(IMAX, JMAX, 4);  // U,E,F有4个分量

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

    // 预估步(前向差分)
    for (int i = 0; i < IMAX; i++) {
        for (int j = 0; j < JMAX; j++) {
            U[i][j][0] = rho[i][j];                                                              // U1
            U[i][j][1] = rho[i][j] * u[i][j];                                                    // U2
            U[i][j][2] = rho[i][j] * v[i][j];                                                    // U3
            U[i][j][3] = rho[i][j] * (e[i][j] + (u[i][j] * u[i][j] + v[i][j] * v[i][j]) / 2.0);  // U5
        }
    }
}
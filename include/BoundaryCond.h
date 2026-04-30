#pragma once
#include <FArray.h>
#include <constant.h>

Array3D ApplyBoundaryCond(Array3D Q) {  // Q:u,v,p,T
    // 1.前缘与入流条件
    for (int j = 0; j < JMAX; j++) {
        if (j == 0)
            Q[0][j][0] = 0.0;
        else
            Q[0][j][0] = a0 * Ma0;
        Q[0][j][1] = 0.0;
        Q[0][j][2] = p0;
        Q[0][j][3] = T0;
    }
    // 2.上边界
    for (int i = 1; i < IMAX; i++) {
        Q[i][JMAX - 1][0] = a0 * Ma0;
        Q[i][JMAX - 1][1] = 0.0;
        Q[i][JMAX - 1][2] = p0;
        Q[i][JMAX - 1][3] = T0;
    }
    // 3.出流
    for (int j = 1; j < JMAX - 1; j++) {
        Q[IMAX - 1][j][0] = 2 * Q[IMAX - 2][j][0] - Q[IMAX - 3][j][0];
        Q[IMAX - 1][j][1] = 2 * Q[IMAX - 2][j][1] - Q[IMAX - 3][j][1];
        Q[IMAX - 1][j][2] = 2 * Q[IMAX - 2][j][2] - Q[IMAX - 3][j][2];
        Q[IMAX - 1][j][3] = 2 * Q[IMAX - 2][j][3] - Q[IMAX - 3][j][3];
    }

    // 4.物面
    for (int i = 1; i < IMAX; i++) {
        Q[i][0][0] = 0.0;
        Q[i][0][1] = 0.0;
        Q[i][0][2] = 2 * Q[i][1][2] - Q[i][2][2];
        //恒温壁条件
        Q[i][0][3] = Tw;
        //绝热壁条件
        //Q[i][0][3] = (4.0*Q[i][1][3]-Q[i][2][3])/3.0;
    }
    return Q;
}
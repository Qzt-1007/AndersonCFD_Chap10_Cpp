#include <FArray.h>
#include <calfunc.h>
#include <constant.h>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    Array3D Q = createArray3D(IMAX, JMAX, 4, 0.0);  // u,v,p,T

    // 初值条件
    for (int ii = 0; ii < IMAX; ii++) {
        for (int jj = 0; jj < JMAX; jj++) {
            Q[ii][jj][0] = Ma0 * a0;                           // u
            Q[ii][jj][1] = 0.0;                                // v
            Q[ii][jj][2] = p0;                                 // p
            Q[ii][jj][3] = T0;                                 // T
        }
    }

    // 壁面初值条件
    for (int ii = 0; ii < IMAX; ii++) {
        Q[ii][0][0] = 0;
        Q[ii][0][1] = 0;
        Q[ii][0][2] = p0;
        Q[ii][0][3] = Tw;
    }
}
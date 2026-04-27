#include <FArray.h>
#include <Maccormack.h>
#include <calfunc.h>
#include <constant.h>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    Array3D Q = createArray3D(IMAX, JMAX, 4, 0.0);  // u,v,p,T

    // 初值
    for (int ii = 1; ii < IMAX - 1; ii++) {
        for (int jj = 1; jj < JMAX - 1; jj++) {
            Q[ii][jj][0] = Ma0 * a0;  // u
            Q[ii][jj][1] = 0.0;       // v
            Q[ii][jj][2] = p0;        // p
            Q[ii][jj][3] = T0;        // T
        }
    }
    // 边界条件
    Q = ApplyBoundaryCond(Q);
}
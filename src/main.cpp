#include <Conver.h>
#include <FArray.h>
#include <Maccormack.h>
#include <calfunc.h>
#include <constant.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// 输出物面压力分布
void outputSurfacePressure(const Array3D& Q, int Iter) {
    std::cout << "\n========== 迭代步 " << Iter << " ==========\n";
    std::cout << "物面 (j=0) 压力分布 p[i][0]:\n";
    std::cout << std::setw(6) << "i" << std::setw(15) << "p" << std::endl;
    for (int i = 0; i < IMAX; ++i) {
        double p_surface = Q[i][0][2];  // Q[i][j][2] 存储压力
        std::cout << std::setw(6) << i << std::setw(15) << p_surface << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    Array3D Q = createArray3D(IMAX, JMAX, 4);  // u,v,p,T
    Array3D Q1 = createArray3D(IMAX, JMAX, 4);
    Array2D rho = createArray2D(IMAX, JMAX);
    Array2D rho1 = createArray2D(IMAX, JMAX);

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

    for (int Iter = 1; Iter < MaxIter + 1; Iter++) {
        Q1 = Maccormack(Q);
        for (int i = 0; i < IMAX; i++) {
            for (int j = 0; j < JMAX; j++) {
                rho[i][j] = Q[i][j][2] / (R * Q[i][j][3]);
                rho1[i][j] = Q1[i][j][2] / (R * Q1[i][j][3]);
            }
        }
        if (Conver(rho, rho1)) {
            Q = Q1;
            break;
        } else {
            Q = Q1;
        }
        if (Iter == 1 || Iter % 100 == 0) std::cout << "Iter = " << Iter << std::endl;

        if (Iter == 1 || Iter % 1000 == 0) {
            outputSurfacePressure(Q1, Iter);  // 输出迭代后的结果
        }
    }
    return 0;
}
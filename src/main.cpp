#include <Conver.h>
#include <FArray.h>
#include <Maccormack.h>
#include <calfunc.h>
#include <constant.h>

#include <fstream>
#include <iomanip>
#include <iostream>

void saveSurfacePressure(const Array3D& Q, int Iter, const std::string& filename) {
    std::ofstream outFile("../output/" + filename);
    if (!outFile) {
        std::cerr << "警告：无法创建，请确保目录存在！\n";
    }
    outFile << "# Converged surface pressure at iteration " << Iter << "\n";
    outFile << "i" << "\t" << "p" << "\n";
    for (int i = 0; i < IMAX; ++i) outFile << i << "\t" << Q[i][0][2] << "\n";
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
    // 边界条件(默认为恒温壁边界条件)
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
            //saveSurfacePressure(Q, Iter, "pressure_isothermal.csv");
            saveSurfacePressure(Q,Iter,"pressure_adiabatic.csv");
            break;
        } else {
            Q = Q1;
        }
        if (Iter % 100 == 0) std::cout << "Iter = " << Iter << std::endl;
    }
    return 0;
}
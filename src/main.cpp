#include <Conver.h>
#include <FArray.h>
#include <Maccormack.h>
#include <calfunc.h>
#include <constant.h>

#include <fstream>
#include <iomanip>
#include <iostream>

// 输出全场 u,v,p,T 到四个CSV文件（每个文件名带后缀）
void saveAllFields(const Array3D& Q, const std::string& suffix) {
    std::string outputDir = "../output/";
    std::ofstream uFile(outputDir + "u_dist_" + suffix + ".csv");
    std::ofstream vFile(outputDir + "v_dist_" + suffix + ".csv");
    std::ofstream pFile(outputDir + "p_dist_" + suffix + ".csv");
    std::ofstream TFile(outputDir + "T_dist_" + suffix + ".csv");
    uFile<<std::fixed<<std::setprecision(3);
    vFile<<std::fixed<<std::setprecision(3);
    pFile<<std::fixed<<std::setprecision(3);
    TFile<<std::fixed<<std::setprecision(3);

    if (!uFile || !vFile || !pFile || !TFile) {
        std::cerr << "错误：无法创建输出文件，请确保目录 " << outputDir << " 存在！\n";
        return;
    }

    // 表头：第一行显示列索引 i
    auto writeHeader = [](std::ofstream& file) {
        file << "j";
        for (int i = 0; i < IMAX; ++i) file << "," << i;
        file << "\n";
    };
    writeHeader(uFile);
    writeHeader(vFile);
    writeHeader(pFile);
    writeHeader(TFile);
    

    // 写入数据：外层循环 j（行），内层循环 i（列）
    for (int j = 0; j < JMAX; ++j) {
        uFile << j;
        vFile << j;
        pFile << j;
        TFile << j;

        for (int i = 0; i < IMAX; ++i) {
            uFile << "," << Q[i][j][0]/a0;
            vFile << "," << Q[i][j][1]/a0;
            pFile << "," << Q[i][j][2]/p0;
            TFile << "," << Q[i][j][3]/T0;
        }

        uFile << "\n";
        vFile << "\n";
        pFile << "\n";
        TFile << "\n";
    }

    std::cout << "全场结果已保存:" << outputDir
              << "u_dist_" << suffix << ".csv 等四个文件\n";
}

// 封装单个边界条件类型的完整计算流程
void runSimulation(int wallType, const std::string& suffix) {
    std::cout << "\n===== 开始计算："
              << (wallType == 0 ? "等温壁 (isothermal)" : "绝热壁 (adiabatic)")
              << " =====\n";

    // 初始化数组
    Array3D Q = createArray3D(IMAX, JMAX, 4);
    Array3D Q1 = createArray3D(IMAX, JMAX, 4);
    Array2D rho = createArray2D(IMAX, JMAX);
    Array2D rho1 = createArray2D(IMAX, JMAX);

    // 初值（内部点）
    for (int i = 1; i < IMAX - 1; ++i) {
        for (int j = 1; j < JMAX - 1; ++j) {
            Q[i][j][0] = Ma0 * a0;   // u
            Q[i][j][1] = 0.0;        // v
            Q[i][j][2] = p0;         // p
            Q[i][j][3] = T0;         // T
        }
    }

    // 应用边界条件（根据 wallType）
    Q = ApplyBoundaryCond(Q, wallType);

    // 迭代求解
    for (int Iter = 1; Iter <= MaxIter; ++Iter) {
        Q1 = Maccormack(Q,wallType);
        // 计算密度用于收敛判断
        for (int i = 0; i < IMAX; ++i) {
            for (int j = 0; j < JMAX; ++j) {
                rho[i][j] = Q[i][j][2] / (R * Q[i][j][3]);
                rho1[i][j] = Q1[i][j][2] / (R * Q1[i][j][3]);
            }
        }
        if (Conver(rho, rho1)) {
            Q = Q1;
            std::cout << "在第 " << Iter << " 步收敛\n";
            // 输出全场变量
            saveAllFields(Q, suffix);
            break;
        } else {
            Q = Q1;
        }
        if (Iter % 200 == 0)
            std::cout << "Iter = " << Iter << std::endl;
    }
}

int main() {
    // 依次计算等温壁和绝热壁
    runSimulation(0, "iso");   // 等温壁 → 输出文件名带 _iso
    runSimulation(1, "adi");   // 绝热壁 → 输出文件名带 _adi

    std::cout << "\n全部计算完成！\n";
    return 0;
}
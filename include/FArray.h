#pragma once
#include <vector>

// 定义2D/3D数组及其初始化方法

// 2D 数组类型别名
using Array2D = std::vector<std::vector<double>>;
// 3D 数组类型别名 (x, y, component)
using Array3D = std::vector<std::vector<std::vector<double>>>;

Array2D createArray2D(int nx, int ny, double initVal = 0.0) { return Array2D(nx, std::vector<double>(ny, initVal)); }

Array3D createArray3D(int nx, int ny, int ncomp, double initVal = 0.0) {
    return Array3D(nx, std::vector<std::vector<double>>(ny, std::vector<double>(ncomp, initVal)));
}

// 复制二维数组
void copyArray2D(const Array2D& src, Array2D& dst) {
    int nx = src.size();
    int ny = src[0].size();
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}

// 复制三维数组
void copyArray3D(const Array3D& src, Array3D& dst) {
    int nx = src.size();
    int ny = src[0].size();
    int nz = src[0][0].size();
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                dst[i][j][k] = src[i][j][k];
            }
        }
    }
}
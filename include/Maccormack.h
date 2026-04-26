/*麦考马克方法*/
#pragma once
#include <Array.h>
#include <constant.h>

#include <vector>

Array3D Maccormack(Array3D Q)  // Q:u,v,p,T,rho,e,Vis,k
{
    Array2D u, v, p, T, rho, e, Vis, k = createArray2D(IMAX, JMAX, 0.0);
    Array3D U, E, F = createArray3D(IMAX, JMAX, 4);  // U,E,F有4个分量
}
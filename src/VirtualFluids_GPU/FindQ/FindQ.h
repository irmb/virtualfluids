#ifndef FIND_Q_H
#define FIND_Q_H

#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include "Parameter/Parameter.h"

extern "C" void findQ(Parameter* para, int lev);

extern "C" void findKforQ(Parameter* para, int lev);

extern "C" void findQ_MG(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, unsigned int* kk, unsigned int sizeQ, real* QQ, QforBoundaryConditions &QIN);

extern "C" void findKforQ_MG(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, QforBoundaryConditions &QIN);

extern "C" void findQInflow(Parameter* para);

extern "C" void findKforQInflow(Parameter* para);

extern "C" void findQPressInflow(Parameter* para);

extern "C" void findKforQPressInflow(Parameter* para);

extern "C" void findQOutflow(Parameter* para);

extern "C" void findKforQOutflow(Parameter* para);

//extern "C" void findQSchlaff( int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, unsigned int* kk, 
//                              unsigned int sizeQN, real* vxN, real* vyN, real* vzN, real*deltaVN, real* QQN, QforBoundaryConditions &QNin,
//                              unsigned int sizeQS, real* vxS, real* vyS, real* vzS, real*deltaVS, real* QQS, QforBoundaryConditions &QSin,
//                              unsigned int sizeQE, real* vxE, real* vyE, real* vzE, real*deltaVE, real* QQE, QforBoundaryConditions &QEin,
//                              unsigned int sizeQW, real* vxW, real* vyW, real* vzW, real*deltaVW, real* QQW, QforBoundaryConditions &QWin);
//
//extern "C" void findKforQSchlaff(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, QforBoundaryConditions &QN, QforBoundaryConditions &QS, QforBoundaryConditions &QE, QforBoundaryConditions &QW);


extern "C" void findKforQPressX1(Parameter* para, int lev);

extern "C" void findQPressX1(Parameter* para, int lev);

extern "C" void findKforQPressX0(Parameter* para, int lev);

extern "C" void findQPressX0(Parameter* para, int lev);

#endif

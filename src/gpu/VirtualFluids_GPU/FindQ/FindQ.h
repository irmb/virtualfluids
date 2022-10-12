#ifndef FIND_Q_H
#define FIND_Q_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "Parameter/Parameter.h"
#include "BoundaryConditions/BoundaryConditionStructs.h"


void findQ(Parameter* para, int lev);

void findKforQ(Parameter* para, int lev);

void findQ_MG(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, unsigned int* kk, unsigned int sizeQ, real* QQ, QforBoundaryConditions &QIN);

void findKforQ_MG(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, QforBoundaryConditions &QIN);

void findQInflow(Parameter* para);

void findKforQInflow(Parameter* para);

void findQPressInflow(Parameter* para);

void findKforQPressInflow(Parameter* para);

void findQOutflow(Parameter* para);

void findKforQOutflow(Parameter* para);

// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
//void findQSchlaff( int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, unsigned int* kk,
//                              unsigned int sizeQN, real* vxN, real* vyN, real* vzN, real*deltaVN, real* QQN, QforBoundaryConditions &QNin,
//                              unsigned int sizeQS, real* vxS, real* vyS, real* vzS, real*deltaVS, real* QQS, QforBoundaryConditions &QSin,
//                              unsigned int sizeQE, real* vxE, real* vyE, real* vzE, real*deltaVE, real* QQE, QforBoundaryConditions &QEin,
//                              unsigned int sizeQW, real* vxW, real* vyW, real* vzW, real*deltaVW, real* QQW, QforBoundaryConditions &QWin);
//
//void findKforQSchlaff(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, QforBoundaryConditions &QN, QforBoundaryConditions &QS, QforBoundaryConditions &QE, QforBoundaryConditions &QW);


void findKforQPressX1(Parameter* para, int lev);

void findQPressX1(Parameter* para, int lev);

void findKforQPressX0(Parameter* para, int lev);

void findQPressX0(Parameter* para, int lev);

#endif

#ifndef FIND_Q_H
#define FIND_Q_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "Parameter/Parameter.h"

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

void findKforQPressX1(Parameter* para, int lev);

void findQPressX1(Parameter* para, int lev);

void findKforQPressX0(Parameter* para, int lev);

void findQPressX0(Parameter* para, int lev);

#endif

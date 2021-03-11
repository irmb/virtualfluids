#include "CalcMac.h"


static const int E    = 0;
static const int W    = 1;
static const int N    = 2;
static const int S    = 3;
static const int T    = 4;
static const int B    = 5;
static const int NE   = 6;
static const int SW   = 7;
static const int SE   = 8;
static const int NW   = 9;
static const int TE   = 10;
static const int BW   = 11;
static const int BE   = 12;
static const int TW   = 13;
static const int TN   = 14;
static const int BS   = 15;
static const int BN   = 16;
static const int TS   = 17;
static const int TNE  = 18;
static const int TNW  = 19;
static const int TSE  = 20;
static const int TSW  = 21;
static const int BNE  = 22;
static const int BNW  = 23;
static const int BSE  = 24;
static const int BSW  = 25;
static const int REST = 26;

__host__ __device__ real LBM::getDensity(const real *const &f /*[27]*/)
{
    return ((f[TNE] + f[BSW]) + (f[TSE] + f[BNW])) + ((f[BSE] + f[TNW]) + (f[TSW] + f[BNE])) +
           (((f[NE] + f[SW]) + (f[SE] + f[NW])) + ((f[TE] + f[BW]) + (f[BE] + f[TW])) +
            ((f[BN] + f[TS]) + (f[TN] + f[BS]))) +
           ((f[E] + f[W]) + (f[N] + f[S]) + (f[T] + f[B])) + f[REST];
}


__host__ __device__ real LBM::getIncompVelocityX1(const real *const &f /*[27]*/)
{
    return ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[BSE] - f[TNW]) + (f[BNE] - f[TSW]))) +
            (((f[BE] - f[TW]) + (f[TE] - f[BW])) + ((f[SE] - f[NW]) + (f[NE] - f[SW]))) + (f[E] - f[W]));
}


__host__ __device__ real LBM::getIncompVelocityX2(const real *const &f /*[27]*/)
{
    return ((((f[TNE] - f[BSW]) + (f[BNW] - f[TSE])) + ((f[TNW] - f[BSE]) + (f[BNE] - f[TSW]))) +
            (((f[BN] - f[TS]) + (f[TN] - f[BS])) + ((f[NW] - f[SE]) + (f[NE] - f[SW]))) + (f[N] - f[S]));
}


__host__ __device__ real LBM::getIncompVelocityX3(const real *const &f /*[27]*/)
{
    return ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[TNW] - f[BSE]) + (f[TSW] - f[BNE]))) +
            (((f[TS] - f[BN]) + (f[TN] - f[BS])) + ((f[TW] - f[BE]) + (f[TE] - f[BW]))) + (f[T] - f[B]));
}


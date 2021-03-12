#include "MacroscopicQuantities.h"

#include "D3Q27.h"


__host__ __device__ real VF::LBM::getDensity(const real *const &f /*[27]*/)
{
    return ((f[DIR::TNE] + f[DIR::BSW]) + (f[DIR::TSE] + f[DIR::BNW])) + ((f[DIR::BSE] + f[DIR::TNW]) + (f[DIR::TSW] + f[DIR::BNE])) +
           (((f[DIR::NE] + f[DIR::SW]) + (f[DIR::SE] + f[DIR::NW])) + ((f[DIR::TE] + f[DIR::BW]) + (f[DIR::BE] + f[DIR::TW])) +
            ((f[DIR::BN] + f[DIR::TS]) + (f[DIR::TN] + f[DIR::BS]))) +
           ((f[DIR::E] + f[DIR::W]) + (f[DIR::N] + f[DIR::S]) + (f[DIR::T] + f[DIR::B])) + f[DIR::REST];
}


real VF::LBM::getIncompressibleVelocityX1(const real *const &f /*[27]*/)
{
    return ((((f[DIR::TNE] - f[DIR::BSW]) + (f[DIR::TSE] - f[DIR::BNW])) + ((f[DIR::BSE] - f[DIR::TNW]) + (f[DIR::BNE] - f[DIR::TSW]))) +
            (((f[DIR::BE] - f[DIR::TW]) + (f[DIR::TE] - f[DIR::BW])) + ((f[DIR::SE] - f[DIR::NW]) + (f[DIR::NE] - f[DIR::SW]))) + (f[DIR::E] - f[DIR::W]));
}


real VF::LBM::getIncompressibleVelocityX2(const real *const &f /*[27]*/)
{
    return ((((f[DIR::TNE] - f[DIR::BSW]) + (f[DIR::BNW] - f[DIR::TSE])) + ((f[DIR::TNW] - f[DIR::BSE]) + (f[DIR::BNE] - f[DIR::TSW]))) +
            (((f[DIR::BN] - f[DIR::TS]) + (f[DIR::TN] - f[DIR::BS])) + ((f[DIR::NW] - f[DIR::SE]) + (f[DIR::NE] - f[DIR::SW]))) + (f[DIR::N] - f[DIR::S]));
}


real VF::LBM::getIncompressibleVelocityX3(const real *const &f /*[27]*/)
{
    return ((((f[DIR::TNE] - f[DIR::BSW]) + (f[DIR::TSE] - f[DIR::BNW])) + ((f[DIR::TNW] - f[DIR::BSE]) + (f[DIR::TSW] - f[DIR::BNE]))) +
            (((f[DIR::TS] - f[DIR::BN]) + (f[DIR::TN] - f[DIR::BS])) + ((f[DIR::TW] - f[DIR::BE]) + (f[DIR::TE] - f[DIR::BW]))) + (f[DIR::T] - f[DIR::B]));
}


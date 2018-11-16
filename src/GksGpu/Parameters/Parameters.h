#ifndef Parameters_H
#define Parameters_H

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include <VirtualFluidsDefinitions.h>

enum class VF_PUBLIC ViscosityModel{
    constant,
    sutherlandsLaw
};

struct  VF_PUBLIC Parameters
{

    real mu = 0.01;
    real K  = 1;
    real Pr = 1.0;

    real D  = 0.01;

    real dt = 0.01;

    real dx = 0.01;

    Vec3 force;

    real lambdaRef = 1.0;

    ViscosityModel viscosityModel = ViscosityModel::constant;

    real boussinesqT0   = 1.0;
    real boussinesqBeta = 1.0;
};

#endif

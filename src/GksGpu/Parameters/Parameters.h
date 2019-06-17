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

    real mu = real(0.01);
    real K  = real(2.0);
    real Pr = real(1.0);
    real D  = real(0.01);

    real dt = real(0.01);
    real dx = real(0.01);

    Vec3 force;

    real lambdaRef = real(1.0);

    real rhoRef = real(1.0);

    ViscosityModel viscosityModel = ViscosityModel::constant;

    real boussinesqT0   = real(1.0);
    real boussinesqBeta = real(1.0);

    //////////////////////////////////////////////////////////////////////////

    bool useSmagorinsky = false;
    real smagorinskyConstant = real(0.2);

    //////////////////////////////////////////////////////////////////////////

    bool enableReaction = false;

    real heatOfReaction = real(8000.0); // kJ / kmol  

    bool useReactionLimiter      = false;
    bool useTemperatureLimiter   = false;
    bool usePassiveScalarLimiter = false;

    real reactionLimiter      = real(1.005);
    real temperatureLimiter   = real(1.0e-3);
    real passiveScalarLimiter = real(0.1);
};

#endif

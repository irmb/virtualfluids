#ifndef LBM_CUMULANT_CHIMERA_PARAMETER_H
#define LBM_CUMULANT_CHIMERA_PARAMETER_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/Core/DataTypes.h>

#include "Distribution27.h"

namespace vf
{
namespace lbm
{


struct CumulantChimeraParameter
{
    Distribution27& distribution;
    real omega;
    real* forces;
};



}
}

#endif

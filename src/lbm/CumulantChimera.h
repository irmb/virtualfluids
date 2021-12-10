#ifndef LBM_CUMULANT_CHIMERA_H
#define LBM_CUMULANT_CHIMERA_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/Core/DataTypes.h>

#include "KernelParameter.h"

namespace vf
{
namespace lbm
{

__host__ __device__ void setRelaxationRatesK17(real omega, real &OxxPyyPzz, real &OxyyPxzz, real &OxyyMxzz, real &Oxyz,
                                               real &O4, real &O5, real &O6);

__host__ __device__ void setRelaxationRatesK15(real omega, real &OxxPyyPzz, real &OxyyPxzz, real &OxyyMxzz, real &Oxyz,
                                               real &O4, real &O5, real &O6);

using RelaxationRatesFunctor = void(*)(real omega, real &OxxPyyPzz, real &OxyyPxzz, real &OxyyMxzz, real &Oxyz,
                                       real &O4, real &O5, real &O6);


__host__ __device__ void cumulantChimera(KernelParameter parameter, RelaxationRatesFunctor setRelaxationRates);

}
}
#endif

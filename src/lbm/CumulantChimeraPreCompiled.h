#ifndef LBM_CUMULANT_CHIMERA_PRE_H
#define LBM_CUMULANT_CHIMERA_PRE_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <cmath>

#include <basics/Core/DataTypes.h>
#include <basics/Core/RealConstants.h>

#include "Chimera.h"
#include "MacroscopicQuantities.h"

namespace vf
{
namespace lbm
{


struct Distribution27
{
    real f[27];

    inline __host__ __device__ real getDensity_() const
    {
        return getDensity(f);
    }
};




inline __host__ __device__ real abs_internal(real value)
{
#ifdef __CUDA_ARCH__
    return ::abs(value);
#else
    return std::abs(value);
#endif
}



//////////////////////////////////////////////////////////////////////////
//! Cumulant K17 Kernel is based on \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
//! and \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
//////////////////////////////////////////////////////////////////////////
__host__ __device__ void cumulantChimeraK17(Distribution27& distribution, real omega, real* forces);

}
}
#endif

#ifndef LBM_KERNEL_PARAMETER_H
#define LBM_KERNEL_PARAMETER_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/Core/DataTypes.h>


namespace vf::lbm
{

struct Distribution27
{
    real f[27];

    __host__ __device__ real getDensity_() const;
};


__host__ __device__ real abs_internal(real value);


struct KernelParameter
{
    Distribution27& distribution;
    real omega;
    real* forces;
};


}

#endif

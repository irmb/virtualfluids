#ifndef LBM_BGK_H
#define LBM_BGK_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include <basics/DataTypes.h>

#include "KernelParameter.h"

namespace vf
{
namespace lbm
{

__host__ __device__ void bgk(KernelParameter parameter);

}
}
#endif

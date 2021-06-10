#include "KernelParameter.h"


#include "MacroscopicQuantities.h"


namespace vf
{
namespace lbm
{



inline __host__ __device__ real Distribution27::getDensity_() const
{
    return getDensity(f);
}



__host__ __device__ real abs_internal(real value)
{
#ifdef __CUDA_ARCH__
    return ::abs(value);
#else
    return std::abs(value);
#endif
}


}
}
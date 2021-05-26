#ifndef LBM_CUMULANT_CHIMERA_K15_H
#define LBM_CUMULANT_CHIMERA_K15_H

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

//////////////////////////////////////////////////////////////////////////
//! Cumulant K17 Kernel is based on \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
//! and \ref
//! <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
//////////////////////////////////////////////////////////////////////////
__host__ __device__ void cumulantChimeraK15(Distribution27& distribution, real omega, real* forces);

}
}
#endif

#ifndef LBM_CUMULANT_CHIMERA_K17_H
#define LBM_CUMULANT_CHIMERA_K17_H

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#include "CumulantChimeraParameter.h"

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
__host__ __device__ void cumulantChimeraK17(CumulantChimeraParameter parameter);

}
}
#endif

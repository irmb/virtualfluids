#ifndef LBM_BACKWARDCHIMERA_H
#define LBM_BACKWARDCHIMERA_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#endif

#include "Core/DataTypes.h"

namespace VF
{
namespace LBM
{

////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardInverseChimeraWithK
//! Transformation from distributions to central moments according to Eq. (6)-(14) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a> Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
__host__ __device__ void forwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv,
                                                             real v2, real Kinverse, real K);
////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardInverseChimeraWithK
//! Transformation from central moments to distributions according to Eq. (57)-(65) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a> ] Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
__host__ __device__ void backwardInverseChimeraWithK(real &mfa, real &mfb, real &mfc, real vv,
                                                              real v2, real Kinverse, real K);
////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardChimera
//! Transformation from distributions to central moments according to Eq. (6)-(14) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a> for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations. Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
__host__ __device__ void forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2);
////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardChimera
//! Transformation from central moments to distributions according to Eq. (57)-(65) in
//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a> for \f$ K_{abc}=0 \f$. This is to avoid unnessary floating point operations. Modified for lower round-off errors.
////////////////////////////////////////////////////////////////////////////////
__host__ __device__ void backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2);

}
}

#endif

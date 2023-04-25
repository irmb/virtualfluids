//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file CalcMac27.cu
//! \ingroup GPU
//! \author Martin Schoenherr, Soeren Peters
//======================================================================================
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "lbm/MacroscopicQuantities.h"

#include "Kernel/Utilities/DistributionHelper.cuh"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMac27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    const unsigned int tx = threadIdx.x;    // Thread index = lokaler i index
    const unsigned int by = blockIdx.x;     // Block index x
    const unsigned int bz = blockIdx.y;     // Block index y
    const unsigned int x = tx + STARTOFFX;  // Globaler x-Index
    const unsigned int y = by + STARTOFFY;  // Globaler y-Index
    const unsigned int z = bz + STARTOFFZ;  // Globaler z-Index
 
    const unsigned nx = blockDim.x + 2 * STARTOFFX;
    const unsigned ny = gridDim.x + 2 * STARTOFFY;
 
    const unsigned int k = nx*(ny*z + y) + x; // Zugriff auf arrays im device
 
 
    if(k >= numberOfLBnodes)
        return;
 
    if(!isValidFluidNode(geoD[k]))
       return;
 
    rhoD[k] = c0o1;
    vxD[k]  = c0o1;
    vyD[k]  = c0o1;
    vzD[k]  = c0o1;
 
    DistributionWrapper distr_wrapper(distributions, numberOfLBnodes, isEvenTimestep, k, neighborX, neighborY, neighborZ);
    const auto& distribution = distr_wrapper.distribution;
 
    rhoD[k] = vf::lbm::getDensity(distribution.f);
    vxD[k] = vf::lbm::getIncompressibleVelocityX1(distribution.f);
    vyD[k] = vf::lbm::getIncompressibleVelocityX2(distribution.f);
    vzD[k] = vf::lbm::getIncompressibleVelocityX3(distribution.f);
}





////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMacSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();
   
    //////////////////////////////////////////////////////////////////////////
    if(nodeIndex<numberOfLBnodes)
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
       
        //////////////////////////////////////////////////////////////////////////
        //index
        unsigned int kzero= nodeIndex;
        unsigned int ke   = nodeIndex;
        unsigned int kw   = neighborX[nodeIndex];
        unsigned int kn   = nodeIndex;
        unsigned int ks   = neighborY[nodeIndex];
        unsigned int kt   = nodeIndex;
        unsigned int kb   = neighborZ[nodeIndex];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = nodeIndex;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = nodeIndex;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = nodeIndex;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
       
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            rhoD[nodeIndex] = 
                (dist.f[DIR_P00])[ke  ]+ (dist.f[DIR_M00])[kw  ]+ 
                (dist.f[DIR_0P0])[kn  ]+ (dist.f[DIR_0M0])[ks  ]+
                (dist.f[DIR_00P])[kt  ]+ (dist.f[DIR_00M])[kb  ]+
                (dist.f[DIR_PP0])[kne ]+ (dist.f[DIR_MM0])[ksw ]+
                (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_P0P])[kte ]+ (dist.f[DIR_M0M])[kbw ]+
                (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_0PP])[ktn ]+ (dist.f[DIR_0MM])[kbs ]+
                (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ]+
                (dist.f[DIR_000])[kzero]+ 
                (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]+ (dist.f[DIR_MMM])[kbsw]+ 
                (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw];
           
            vxD[nodeIndex] =
                (dist.f[DIR_P00])[ke  ]- (dist.f[DIR_M00])[kw  ]+ 
                (dist.f[DIR_PP0])[kne ]- (dist.f[DIR_MM0])[ksw ]+
                (dist.f[DIR_PM0])[kse ]- (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_P0P])[kte ]- (dist.f[DIR_M0M])[kbw ]+
                (dist.f[DIR_P0M])[kbe ]- (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_PPP])[ktne]- (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]- (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]+ 
                (dist.f[DIR_PMM])[kbse]- (dist.f[DIR_MPM])[kbnw];
           
            vyD[nodeIndex] =
                (dist.f[DIR_0P0])[kn  ]- (dist.f[DIR_0M0])[ks  ]+
                (dist.f[DIR_PP0])[kne ]- (dist.f[DIR_MM0])[ksw ]-
                (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_0PP])[ktn ]- (dist.f[DIR_0MM])[kbs ]+
                (dist.f[DIR_0PM])[kbn ]- (dist.f[DIR_0MP])[kts ]+
                (dist.f[DIR_PPP])[ktne]- (dist.f[DIR_MMP])[ktsw]- 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]- 
                (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw];
           
            vzD[nodeIndex] =
                (dist.f[DIR_00P])[kt  ]- (dist.f[DIR_00M])[kb  ]+
                (dist.f[DIR_P0P])[kte ]- (dist.f[DIR_M0M])[kbw ]-
                (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_0PP])[ktn ]- (dist.f[DIR_0MM])[kbs ]-
                (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ]+
                (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]- 
                (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]- 
                (dist.f[DIR_PMM])[kbse]- (dist.f[DIR_MPM])[kbnw];
           
            pressD[nodeIndex] =
                ((dist.f[DIR_P00])[ke  ]+ (dist.f[DIR_M00])[kw  ]+ 
                (dist.f[DIR_0P0])[kn  ]+ (dist.f[DIR_0M0])[ks  ]+
                (dist.f[DIR_00P])[kt  ]+ (dist.f[DIR_00M])[kb  ]+
                2.f*(
                (dist.f[DIR_PP0])[kne ]+ (dist.f[DIR_MM0])[ksw ]+
                (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_P0P])[kte ]+ (dist.f[DIR_M0M])[kbw ]+
                (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_0PP])[ktn ]+ (dist.f[DIR_0MM])[kbs ]+
                (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ])+
                3.f*(
                (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]+ (dist.f[DIR_MMM])[kbsw]+ 
                (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+c0o1*rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]; // times zero for incompressible case   
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
       }
    }
}
////////////////////////////////////////////////////////////////////////////////
































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMacCompSP27(
    real *vxD,
    real *vyD,
    real *vzD,
    real *rhoD,
    real *pressD,
    unsigned int *geoD,
    unsigned int *neighborX,
    unsigned int *neighborY,
    unsigned int *neighborZ,
    unsigned long long numberOfLBnodes,
    real *distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    if(nodeIndex >= numberOfLBnodes)
        return;

    pressD[nodeIndex] = c0o1;
    rhoD[nodeIndex]   = c0o1;
    vxD[nodeIndex]    = c0o1;
    vyD[nodeIndex]    = c0o1;
    vzD[nodeIndex]    = c0o1;

    if (!isValidFluidNode(geoD[nodeIndex]))
        return;

    DistributionWrapper distr_wrapper(distributions, numberOfLBnodes, isEvenTimestep, nodeIndex, neighborX, neighborY, neighborZ);
    const auto &distribution = distr_wrapper.distribution;

    rhoD[nodeIndex]   = vf::lbm::getDensity(distribution.f);
    vxD[nodeIndex]    = vf::lbm::getCompressibleVelocityX1(distribution.f, rhoD[nodeIndex]);
    vyD[nodeIndex]    = vf::lbm::getCompressibleVelocityX2(distribution.f, rhoD[nodeIndex]);
    vzD[nodeIndex]    = vf::lbm::getCompressibleVelocityX3(distribution.f, rhoD[nodeIndex]);
    pressD[nodeIndex] = vf::lbm::getPressure(distribution.f, rhoD[nodeIndex], vxD[nodeIndex], vyD[nodeIndex], vzD[nodeIndex]); 
}




































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMedSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    if( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
        
        //////////////////////////////////////////////////////////////////////////
        //index
        unsigned int kzero= nodeIndex;
        unsigned int ke   = nodeIndex;
        unsigned int kw   = neighborX[nodeIndex];
        unsigned int kn   = nodeIndex;
        unsigned int ks   = neighborY[nodeIndex];
        unsigned int kt   = nodeIndex;
        unsigned int kb   = neighborZ[nodeIndex];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = nodeIndex;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = nodeIndex;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = nodeIndex;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
        
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            rhoD[nodeIndex] =
                (dist.f[DIR_P00])[ke  ]+ (dist.f[DIR_M00])[kw  ]+ 
                (dist.f[DIR_0P0])[kn  ]+ (dist.f[DIR_0M0])[ks  ]+
                (dist.f[DIR_00P])[kt  ]+ (dist.f[DIR_00M])[kb  ]+
                (dist.f[DIR_PP0])[kne ]+ (dist.f[DIR_MM0])[ksw ]+
                (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_P0P])[kte ]+ (dist.f[DIR_M0M])[kbw ]+
                (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_0PP])[ktn ]+ (dist.f[DIR_0MM])[kbs ]+
                (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ]+
                (dist.f[DIR_000])[kzero]+ 
                (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]+ (dist.f[DIR_MMM])[kbsw]+ 
                (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw]+
                RHO;
            
            vxD[nodeIndex] =
                (dist.f[DIR_P00])[ke  ]- (dist.f[DIR_M00])[kw  ]+ 
                (dist.f[DIR_PP0])[kne ]- (dist.f[DIR_MM0])[ksw ]+
                (dist.f[DIR_PM0])[kse ]- (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_P0P])[kte ]- (dist.f[DIR_M0M])[kbw ]+
                (dist.f[DIR_P0M])[kbe ]- (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_PPP])[ktne]- (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]- (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]+ 
                (dist.f[DIR_PMM])[kbse]- (dist.f[DIR_MPM])[kbnw]+
                VX;
            
            vyD[nodeIndex] =
                (dist.f[DIR_0P0])[kn  ]- (dist.f[DIR_0M0])[ks  ]+
                (dist.f[DIR_PP0])[kne ]- (dist.f[DIR_MM0])[ksw ]-
                (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_0PP])[ktn ]- (dist.f[DIR_0MM])[kbs ]+
                (dist.f[DIR_0PM])[kbn ]- (dist.f[DIR_0MP])[kts ]+
                (dist.f[DIR_PPP])[ktne]- (dist.f[DIR_MMP])[ktsw]- 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]- 
                (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw]+
                VY;
            
            vzD[nodeIndex] =
                (dist.f[DIR_00P])[kt  ]- (dist.f[DIR_00M])[kb  ]+
                (dist.f[DIR_P0P])[kte ]- (dist.f[DIR_M0M])[kbw ]-
                (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_0PP])[ktn ]- (dist.f[DIR_0MM])[kbs ]-
                (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ]+
                (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]- 
                (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]- 
                (dist.f[DIR_PMM])[kbse]- (dist.f[DIR_MPM])[kbnw]+
                VZ;
            
            pressD[nodeIndex] =
                ((dist.f[DIR_P00])[ke  ]+ (dist.f[DIR_M00])[kw  ]+ 
                (dist.f[DIR_0P0])[kn  ]+ (dist.f[DIR_0M0])[ks  ]+
                (dist.f[DIR_00P])[kt  ]+ (dist.f[DIR_00M])[kb  ]+
                c2o1*(
                (dist.f[DIR_PP0])[kne ]+ (dist.f[DIR_MM0])[ksw ]+
                (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_P0P])[kte ]+ (dist.f[DIR_M0M])[kbw ]+
                (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_0PP])[ktn ]+ (dist.f[DIR_0MM])[kbs ]+
                (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ])+
                c3o1*(
                (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]+ (dist.f[DIR_MMM])[kbsw]+ 
                (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]+
                PRESS;    
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
        }
    }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMedCompSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    if( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
        
        //////////////////////////////////////////////////////////////////////////
        //index
        //unsigned int kzero= k;
        unsigned int ke   = nodeIndex;
        unsigned int kw   = neighborX[nodeIndex];
        unsigned int kn   = nodeIndex;
        unsigned int ks   = neighborY[nodeIndex];
        unsigned int kt   = nodeIndex;
        unsigned int kb   = neighborZ[nodeIndex];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = nodeIndex;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = nodeIndex;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = nodeIndex;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
        
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            real mfcbb = (dist.f[DIR_P00])[nodeIndex];//[ke   ];
            real mfabb = (dist.f[DIR_M00])[kw];//[kw   ];  
            real mfbcb = (dist.f[DIR_0P0])[nodeIndex];//[kn   ];
            real mfbab = (dist.f[DIR_0M0])[ks];//[ks   ];  
            real mfbbc = (dist.f[DIR_00P])[nodeIndex];//[kt   ];
            real mfbba = (dist.f[DIR_00M])[kb];//[kb   ];  
            real mfccb = (dist.f[DIR_PP0])[nodeIndex];//[kne  ];  
            real mfaab = (dist.f[DIR_MM0])[ksw];//[ksw  ];
            real mfcab = (dist.f[DIR_PM0])[ks];//[kse  ]; 
            real mfacb = (dist.f[DIR_MP0])[kw];//[knw  ]; 
            real mfcbc = (dist.f[DIR_P0P])[nodeIndex];//[kte  ];  
            real mfaba = (dist.f[DIR_M0M])[kbw];//[kbw  ];
            real mfcba = (dist.f[DIR_P0M])[kb];//[kbe  ]; 
            real mfabc = (dist.f[DIR_M0P])[kw];//[ktw  ]; 
            real mfbcc = (dist.f[DIR_0PP])[nodeIndex];//[ktn  ];  
            real mfbaa = (dist.f[DIR_0MM])[kbs];//[kbs  ];
            real mfbca = (dist.f[DIR_0PM])[kb];//[kbn  ]; 
            real mfbac = (dist.f[DIR_0MP])[ks];//[kts  ]; 
            real mfbbb = (dist.f[DIR_000])[nodeIndex];//[kzero];
            real mfccc = (dist.f[DIR_PPP])[nodeIndex];//[ktne ]; 
            real mfaac = (dist.f[DIR_MMP])[ksw];//[ktsw ]; 
            real mfcac = (dist.f[DIR_PMP])[ks];//[ktse ];
            real mfacc = (dist.f[DIR_MPP])[kw];//[ktnw ];
            real mfcca = (dist.f[DIR_PPM])[kb];//[kbne ];
            real mfaaa = (dist.f[DIR_MMM])[kbsw];//[kbsw ];
            real mfcaa = (dist.f[DIR_PMM])[kbs];//[kbse ]; 
            real mfaca = (dist.f[DIR_MPM])[kbw];//[kbnw ]; 
            ////////////////////////////////////////////////////////////////////////////////////
            real drho = 
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

            real rho = c1o1 + drho;

            rhoD[nodeIndex] = drho + RHO;

            vxD[nodeIndex] = 
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                (mfcbb - mfabb)) / rho) + VX;
            vyD[nodeIndex] = 
                (((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                (mfbcb - mfbab)) / rho) + VY;
            vzD[nodeIndex] = 
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                (mfbbc - mfbba)) / rho) + VZ;

            pressD[nodeIndex]  =
                ((dist.f[DIR_P00])[ke  ]+ (dist.f[DIR_M00])[kw  ]+ 
                (dist.f[DIR_0P0])[kn  ]+ (dist.f[DIR_0M0])[ks  ]+
                (dist.f[DIR_00P])[kt  ]+ (dist.f[DIR_00M])[kb  ]+
                c2o1*(
                (dist.f[DIR_PP0])[kne ]+ (dist.f[DIR_MM0])[ksw ]+
                (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                (dist.f[DIR_P0P])[kte ]+ (dist.f[DIR_M0M])[kbw ]+
                (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                (dist.f[DIR_0PP])[ktn ]+ (dist.f[DIR_0MM])[kbs ]+
                (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ])+
                c3o1*(
                (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                (dist.f[DIR_PPM])[kbne]+ (dist.f[DIR_MMM])[kbsw]+ 
                (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]+
                PRESS;    
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
        }
    }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMedCompAD27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    real* concD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    real* distributionsAD,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    if ( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist, distAD;
        getPointersToDistributions(dist,   distributions,   numberOfLBnodes, isEvenTimestep);
        getPointersToDistributions(distAD, distributionsAD, numberOfLBnodes, isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////
        //index
        //unsigned int kzero = k;
        unsigned int ke = nodeIndex;
        unsigned int kw = neighborX[nodeIndex];
        unsigned int kn = nodeIndex;
        unsigned int ks = neighborY[nodeIndex];
        unsigned int kt = nodeIndex;
        unsigned int kb = neighborZ[nodeIndex];
        unsigned int ksw = neighborY[kw];
        unsigned int kne = nodeIndex;
        unsigned int kse = ks;
        unsigned int knw = kw;
        unsigned int kbw = neighborZ[kw];
        unsigned int kte = nodeIndex;
        unsigned int kbe = kb;
        unsigned int ktw = kw;
        unsigned int kbs = neighborZ[ks];
        unsigned int ktn = nodeIndex;
        unsigned int kbn = kb;
        unsigned int kts = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        real CONC  = concD[nodeIndex];
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        concD[nodeIndex]  = c0o1;
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
        
        if (geoD[nodeIndex] == GEO_FLUID)
        {
            real mfcbb = (dist.f[DIR_P00])[nodeIndex];//[ke   ];
            real mfabb = (dist.f[DIR_M00])[kw];//[kw   ];  
            real mfbcb = (dist.f[DIR_0P0])[nodeIndex];//[kn   ];
            real mfbab = (dist.f[DIR_0M0])[ks];//[ks   ];  
            real mfbbc = (dist.f[DIR_00P])[nodeIndex];//[kt   ];
            real mfbba = (dist.f[DIR_00M])[kb];//[kb   ];  
            real mfccb = (dist.f[DIR_PP0])[nodeIndex];//[kne  ];  
            real mfaab = (dist.f[DIR_MM0])[ksw];//[ksw  ];
            real mfcab = (dist.f[DIR_PM0])[ks];//[kse  ]; 
            real mfacb = (dist.f[DIR_MP0])[kw];//[knw  ]; 
            real mfcbc = (dist.f[DIR_P0P])[nodeIndex];//[kte  ];  
            real mfaba = (dist.f[DIR_M0M])[kbw];//[kbw  ];
            real mfcba = (dist.f[DIR_P0M])[kb];//[kbe  ]; 
            real mfabc = (dist.f[DIR_M0P])[kw];//[ktw  ]; 
            real mfbcc = (dist.f[DIR_0PP])[nodeIndex];//[ktn  ];  
            real mfbaa = (dist.f[DIR_0MM])[kbs];//[kbs  ];
            real mfbca = (dist.f[DIR_0PM])[kb];//[kbn  ]; 
            real mfbac = (dist.f[DIR_0MP])[ks];//[kts  ]; 
            real mfbbb = (dist.f[DIR_000])[nodeIndex];//[kzero];
            real mfccc = (dist.f[DIR_PPP])[nodeIndex];//[ktne ]; 
            real mfaac = (dist.f[DIR_MMP])[ksw];//[ktsw ]; 
            real mfcac = (dist.f[DIR_PMP])[ks];//[ktse ];
            real mfacc = (dist.f[DIR_MPP])[kw];//[ktnw ];
            real mfcca = (dist.f[DIR_PPM])[kb];//[kbne ];
            real mfaaa = (dist.f[DIR_MMM])[kbsw];//[kbsw ];
            real mfcaa = (dist.f[DIR_PMM])[kbs];//[kbse ]; 
            real mfaca = (dist.f[DIR_MPM])[kbw];//[kbnw ]; 
            ////////////////////////////////////////////////////////////////////////////////////
            real drho =
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) + mfbbb;
            real rho = c1o1 + drho;
            ////////////////////////////////////////////////////////////////////////////////////
            
            rhoD[nodeIndex] = drho + RHO;
            
            vxD[nodeIndex] =
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                    (mfcbb - mfabb)) / rho) + VX;
            
            vyD[nodeIndex] =
                (((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                    (mfbcb - mfbab)) / rho) + VY;
            
            vzD[nodeIndex] =
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                    (mfbbc - mfbba)) / rho) + VZ;
            
            pressD[nodeIndex] = 
                ((dist.f[DIR_P00])[ke] + (dist.f[DIR_M00])[kw] +
                 (dist.f[DIR_0P0])[kn] + (dist.f[DIR_0M0])[ks] +
                 (dist.f[DIR_00P])[kt] + (dist.f[DIR_00M])[kb] +
                 c2o1*(
                 (dist.f[DIR_PP0])[kne] + (dist.f[DIR_MM0])[ksw] +
                 (dist.f[DIR_PM0])[kse] + (dist.f[DIR_MP0])[knw] +
                 (dist.f[DIR_P0P])[kte] + (dist.f[DIR_M0M])[kbw] +
                 (dist.f[DIR_P0M])[kbe] + (dist.f[DIR_M0P])[ktw] +
                 (dist.f[DIR_0PP])[ktn] + (dist.f[DIR_0MM])[kbs] +
                 (dist.f[DIR_0PM])[kbn] + (dist.f[DIR_0MP])[kts]) +
                 c3o1*(
                 (dist.f[DIR_PPP])[ktne] + (dist.f[DIR_MMP])[ktsw] +
                 (dist.f[DIR_PMP])[ktse] + (dist.f[DIR_MPP])[ktnw] +
                 (dist.f[DIR_PPM])[kbne] + (dist.f[DIR_MMM])[kbsw] +
                 (dist.f[DIR_PMM])[kbse] + (dist.f[DIR_MPM])[kbnw]) -
                 rhoD[nodeIndex] - (vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1 + rhoD[nodeIndex])) * c1o2 + rhoD[nodeIndex] +
                 PRESS;
                 //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
            //////////////////////////////////////////////////////////////////////////
            mfcbb = (distAD.f[DIR_P00])[nodeIndex   ];
            mfabb = (distAD.f[DIR_M00])[kw  ];
            mfbcb = (distAD.f[DIR_0P0])[nodeIndex   ];
            mfbab = (distAD.f[DIR_0M0])[ks  ];
            mfbbc = (distAD.f[DIR_00P])[nodeIndex   ];
            mfbba = (distAD.f[DIR_00M])[kb  ];
            mfccb = (distAD.f[DIR_PP0])[nodeIndex   ];
            mfaab = (distAD.f[DIR_MM0])[ksw ];
            mfcab = (distAD.f[DIR_PM0])[ks  ];
            mfacb = (distAD.f[DIR_MP0])[kw  ];
            mfcbc = (distAD.f[DIR_P0P])[nodeIndex   ];
            mfaba = (distAD.f[DIR_M0M])[kbw ];
            mfcba = (distAD.f[DIR_P0M])[kb  ];
            mfabc = (distAD.f[DIR_M0P])[kw  ];
            mfbcc = (distAD.f[DIR_0PP])[nodeIndex   ];
            mfbaa = (distAD.f[DIR_0MM])[kbs ];
            mfbca = (distAD.f[DIR_0PM])[kb  ];
            mfbac = (distAD.f[DIR_0MP])[ks  ];
            mfbbb = (distAD.f[DIR_000])[nodeIndex   ];
            mfccc = (distAD.f[DIR_PPP])[nodeIndex   ];
            mfaac = (distAD.f[DIR_MMP])[ksw ];
            mfcac = (distAD.f[DIR_PMP])[ks  ];
            mfacc = (distAD.f[DIR_MPP])[kw  ];
            mfcca = (distAD.f[DIR_PPM])[kb  ];
            mfaaa = (distAD.f[DIR_MMM])[kbsw];
            mfcaa = (distAD.f[DIR_PMM])[kbs ];
            mfaca = (distAD.f[DIR_MPM])[kbw ];
            //////////////////////////////////////////////////////////////////////////
            concD[nodeIndex] = 
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa)   + (mfaac + mfcca))) +
                 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba)   + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) +  mfbbb + CONC;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMacMedSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned int tdiff,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    if(nodeIndex<numberOfLBnodes)
    {
        //////////////////////////////////////////////////////////////////////////
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
       
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            rhoD[nodeIndex]    =   RHO   / tdiff;
            vxD[nodeIndex]     =   VX    / tdiff;
            vyD[nodeIndex]     =   VY    / tdiff;
            vzD[nodeIndex]     =   VZ    / tdiff;
            pressD[nodeIndex]  =   PRESS / tdiff;    
        }
    }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBResetMedianValuesSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    if ( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex] = c0o1;
        vxD[nodeIndex] = c0o1;
        vyD[nodeIndex] = c0o1;
        vzD[nodeIndex] = c0o1;
    }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBResetMedianValuesAD27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    real* concD,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    if (nodeIndex < numberOfLBnodes)
    {
        concD[nodeIndex]  = c0o1;
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
    }
}
////////////////////////////////////////////////////////////////////////////////





















////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMeasurePoints(
    real* vxMP,
    real* vyMP,
    real* vzMP,
    real* rhoMP,
    unsigned int* kMP,
    unsigned int numberOfPointskMP,
    unsigned int MPClockCycle,
    unsigned int t,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    if( nodeIndex < numberOfPointskMP )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////
        //index
        unsigned int kzero= kMP[nodeIndex];//k;
        unsigned int ke   = kzero;
        unsigned int kw   = neighborX[kzero];
        unsigned int kn   = kzero;
        unsigned int ks   = neighborY[kzero];
        unsigned int kt   = kzero;
        unsigned int kb   = neighborZ[kzero];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = kzero;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = kzero;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = kzero;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = kzero;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
	    unsigned int kMac = nodeIndex*MPClockCycle + t;
	    //////////////////////////////////////////////////////////////////////////
        
        if(geoD[kzero] == GEO_FLUID)
        {
            rhoMP[kMac]= (dist.f[DIR_P00])[ke  ]+ (dist.f[DIR_M00])[kw  ]+ 
                         (dist.f[DIR_0P0])[kn  ]+ (dist.f[DIR_0M0])[ks  ]+
                         (dist.f[DIR_00P])[kt  ]+ (dist.f[DIR_00M])[kb  ]+
                         (dist.f[DIR_PP0])[kne ]+ (dist.f[DIR_MM0])[ksw ]+
                         (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                         (dist.f[DIR_P0P])[kte ]+ (dist.f[DIR_M0M])[kbw ]+
                         (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                         (dist.f[DIR_0PP])[ktn ]+ (dist.f[DIR_0MM])[kbs ]+
                         (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ]+
                         (dist.f[DIR_000])[kzero]+ 
                         (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                         (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                         (dist.f[DIR_PPM])[kbne]+ (dist.f[DIR_MMM])[kbsw]+ 
                         (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw];
           
            vxMP[kMac] = (dist.f[DIR_P00])[ke  ]- (dist.f[DIR_M00])[kw  ]+ 
                         (dist.f[DIR_PP0])[kne ]- (dist.f[DIR_MM0])[ksw ]+
                         (dist.f[DIR_PM0])[kse ]- (dist.f[DIR_MP0])[knw ]+
                         (dist.f[DIR_P0P])[kte ]- (dist.f[DIR_M0M])[kbw ]+
                         (dist.f[DIR_P0M])[kbe ]- (dist.f[DIR_M0P])[ktw ]+
                         (dist.f[DIR_PPP])[ktne]- (dist.f[DIR_MMP])[ktsw]+ 
                         (dist.f[DIR_PMP])[ktse]- (dist.f[DIR_MPP])[ktnw]+ 
                         (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]+ 
                         (dist.f[DIR_PMM])[kbse]- (dist.f[DIR_MPM])[kbnw];
           
            vyMP[kMac] = (dist.f[DIR_0P0])[kn  ]- (dist.f[DIR_0M0])[ks  ]+
                         (dist.f[DIR_PP0])[kne ]- (dist.f[DIR_MM0])[ksw ]-
                         (dist.f[DIR_PM0])[kse ]+ (dist.f[DIR_MP0])[knw ]+
                         (dist.f[DIR_0PP])[ktn ]- (dist.f[DIR_0MM])[kbs ]+
                         (dist.f[DIR_0PM])[kbn ]- (dist.f[DIR_0MP])[kts ]+
                         (dist.f[DIR_PPP])[ktne]- (dist.f[DIR_MMP])[ktsw]- 
                         (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]+ 
                         (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]- 
                         (dist.f[DIR_PMM])[kbse]+ (dist.f[DIR_MPM])[kbnw];
           
            vzMP[kMac] = (dist.f[DIR_00P])[kt  ]- (dist.f[DIR_00M])[kb  ]+
                         (dist.f[DIR_P0P])[kte ]- (dist.f[DIR_M0M])[kbw ]-
                         (dist.f[DIR_P0M])[kbe ]+ (dist.f[DIR_M0P])[ktw ]+
                         (dist.f[DIR_0PP])[ktn ]- (dist.f[DIR_0MM])[kbs ]-
                         (dist.f[DIR_0PM])[kbn ]+ (dist.f[DIR_0MP])[kts ]+
                         (dist.f[DIR_PPP])[ktne]+ (dist.f[DIR_MMP])[ktsw]+ 
                         (dist.f[DIR_PMP])[ktse]+ (dist.f[DIR_MPP])[ktnw]- 
                         (dist.f[DIR_PPM])[kbne]- (dist.f[DIR_MMM])[kbsw]- 
                         (dist.f[DIR_PMM])[kbse]- (dist.f[DIR_MPM])[kbnw];
        }
    }
}
////////////////////////////////////////////////////////////////////////////////





































////////////////////////////////////////////////////////////////////////////////
__global__ void LBSetOutputWallVelocitySP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* vxWall,
    real* vyWall,
    real* vzWall,
    int numberOfWallNodes, 
    int* kWallNodes, 
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* DD,
    bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   //////////////////////////////////////////////////////////////////////////
   if(nodeIndex<numberOfWallNodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KWN  = kWallNodes[nodeIndex];
      //////////////////////////////////////////////////////////////////////////
      vxD[KWN] = 0.0;//vxWall[k];
      vyD[KWN] = 0.0;//vyWall[k];
      vzD[KWN] = 0.0;//vzWall[k];
   }
}






























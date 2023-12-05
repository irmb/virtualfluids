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
//! \author Martin Schoenherr, Soeren Peters
//======================================================================================
#include "MacroscopicQuantities.cuh"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "LBM/GPUHelperFunctions/KernelUtilities.h"
#include "LBM/LB.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

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
                (dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[d000])[kzero]+ 
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vxD[nodeIndex] =
                (dist.f[dP00])[ke  ]- (dist.f[dM00])[kw  ]+ 
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]- (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]- (dist.f[dM0P])[ktw ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]- (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
           
            vyD[nodeIndex] =
                (dist.f[d0P0])[kn  ]- (dist.f[d0M0])[ks  ]+
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]-
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]- (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]- 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vzD[nodeIndex] =
                (dist.f[d00P])[kt  ]- (dist.f[d00M])[kb  ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]-
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]-
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]- 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
           
            pressD[nodeIndex] =
                ((dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                2.f*(
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ])+
                3.f*(
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+c0o1*rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]; // times zero for incompressible case   
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
       }
    }
}

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

    Distributions27 dist;
    vf::gpu::getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
    vf::gpu::ListIndices listIndices(nodeIndex, neighborX, neighborY, neighborZ);

    real distribution[27];
    vf::gpu::getPreCollisionDistribution(distribution, dist, listIndices);

    rhoD[nodeIndex] = vf::lbm::getDensity(distribution);
    vxD[nodeIndex] = vf::lbm::getCompressibleVelocityX1(distribution, rhoD[nodeIndex]);
    vyD[nodeIndex] = vf::lbm::getCompressibleVelocityX2(distribution, rhoD[nodeIndex]);
    vzD[nodeIndex] = vf::lbm::getCompressibleVelocityX3(distribution, rhoD[nodeIndex]);
    pressD[nodeIndex] = vf::lbm::getPressure(distribution, rhoD[nodeIndex], vxD[nodeIndex], vyD[nodeIndex], vzD[nodeIndex]);
}

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
                (dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[d000])[kzero]+ 
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw]+
                RHO;
            
            vxD[nodeIndex] =
                (dist.f[dP00])[ke  ]- (dist.f[dM00])[kw  ]+ 
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]- (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]- (dist.f[dM0P])[ktw ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]- (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw]+
                VX;
            
            vyD[nodeIndex] =
                (dist.f[d0P0])[kn  ]- (dist.f[d0M0])[ks  ]+
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]-
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]- (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]- 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw]+
                VY;
            
            vzD[nodeIndex] =
                (dist.f[d00P])[kt  ]- (dist.f[d00M])[kb  ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]-
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]-
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]- 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw]+
                VZ;
            
            pressD[nodeIndex] =
                ((dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                c2o1*(
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ])+
                c3o1*(
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]+
                PRESS;    
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
        }
    }
}

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
            real mfcbb = (dist.f[dP00])[nodeIndex];//[ke   ];
            real mfabb = (dist.f[dM00])[kw];//[kw   ];  
            real mfbcb = (dist.f[d0P0])[nodeIndex];//[kn   ];
            real mfbab = (dist.f[d0M0])[ks];//[ks   ];  
            real mfbbc = (dist.f[d00P])[nodeIndex];//[kt   ];
            real mfbba = (dist.f[d00M])[kb];//[kb   ];  
            real mfccb = (dist.f[dPP0])[nodeIndex];//[kne  ];  
            real mfaab = (dist.f[dMM0])[ksw];//[ksw  ];
            real mfcab = (dist.f[dPM0])[ks];//[kse  ]; 
            real mfacb = (dist.f[dMP0])[kw];//[knw  ]; 
            real mfcbc = (dist.f[dP0P])[nodeIndex];//[kte  ];  
            real mfaba = (dist.f[dM0M])[kbw];//[kbw  ];
            real mfcba = (dist.f[dP0M])[kb];//[kbe  ]; 
            real mfabc = (dist.f[dM0P])[kw];//[ktw  ]; 
            real mfbcc = (dist.f[d0PP])[nodeIndex];//[ktn  ];  
            real mfbaa = (dist.f[d0MM])[kbs];//[kbs  ];
            real mfbca = (dist.f[d0PM])[kb];//[kbn  ]; 
            real mfbac = (dist.f[d0MP])[ks];//[kts  ]; 
            real mfbbb = (dist.f[d000])[nodeIndex];//[kzero];
            real mfccc = (dist.f[dPPP])[nodeIndex];//[ktne ]; 
            real mfaac = (dist.f[dMMP])[ksw];//[ktsw ]; 
            real mfcac = (dist.f[dPMP])[ks];//[ktse ];
            real mfacc = (dist.f[dMPP])[kw];//[ktnw ];
            real mfcca = (dist.f[dPPM])[kb];//[kbne ];
            real mfaaa = (dist.f[dMMM])[kbsw];//[kbsw ];
            real mfcaa = (dist.f[dPMM])[kbs];//[kbse ]; 
            real mfaca = (dist.f[dMPM])[kbw];//[kbnw ]; 
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
                ((dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                c2o1*(
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ])+
                c3o1*(
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]+
                PRESS;    
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
        }
    }
}

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
            real mfcbb = (dist.f[dP00])[nodeIndex];//[ke   ];
            real mfabb = (dist.f[dM00])[kw];//[kw   ];  
            real mfbcb = (dist.f[d0P0])[nodeIndex];//[kn   ];
            real mfbab = (dist.f[d0M0])[ks];//[ks   ];  
            real mfbbc = (dist.f[d00P])[nodeIndex];//[kt   ];
            real mfbba = (dist.f[d00M])[kb];//[kb   ];  
            real mfccb = (dist.f[dPP0])[nodeIndex];//[kne  ];  
            real mfaab = (dist.f[dMM0])[ksw];//[ksw  ];
            real mfcab = (dist.f[dPM0])[ks];//[kse  ]; 
            real mfacb = (dist.f[dMP0])[kw];//[knw  ]; 
            real mfcbc = (dist.f[dP0P])[nodeIndex];//[kte  ];  
            real mfaba = (dist.f[dM0M])[kbw];//[kbw  ];
            real mfcba = (dist.f[dP0M])[kb];//[kbe  ]; 
            real mfabc = (dist.f[dM0P])[kw];//[ktw  ]; 
            real mfbcc = (dist.f[d0PP])[nodeIndex];//[ktn  ];  
            real mfbaa = (dist.f[d0MM])[kbs];//[kbs  ];
            real mfbca = (dist.f[d0PM])[kb];//[kbn  ]; 
            real mfbac = (dist.f[d0MP])[ks];//[kts  ]; 
            real mfbbb = (dist.f[d000])[nodeIndex];//[kzero];
            real mfccc = (dist.f[dPPP])[nodeIndex];//[ktne ]; 
            real mfaac = (dist.f[dMMP])[ksw];//[ktsw ]; 
            real mfcac = (dist.f[dPMP])[ks];//[ktse ];
            real mfacc = (dist.f[dMPP])[kw];//[ktnw ];
            real mfcca = (dist.f[dPPM])[kb];//[kbne ];
            real mfaaa = (dist.f[dMMM])[kbsw];//[kbsw ];
            real mfcaa = (dist.f[dPMM])[kbs];//[kbse ]; 
            real mfaca = (dist.f[dMPM])[kbw];//[kbnw ]; 
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
                ((dist.f[dP00])[ke] + (dist.f[dM00])[kw] +
                 (dist.f[d0P0])[kn] + (dist.f[d0M0])[ks] +
                 (dist.f[d00P])[kt] + (dist.f[d00M])[kb] +
                 c2o1*(
                 (dist.f[dPP0])[kne] + (dist.f[dMM0])[ksw] +
                 (dist.f[dPM0])[kse] + (dist.f[dMP0])[knw] +
                 (dist.f[dP0P])[kte] + (dist.f[dM0M])[kbw] +
                 (dist.f[dP0M])[kbe] + (dist.f[dM0P])[ktw] +
                 (dist.f[d0PP])[ktn] + (dist.f[d0MM])[kbs] +
                 (dist.f[d0PM])[kbn] + (dist.f[d0MP])[kts]) +
                 c3o1*(
                 (dist.f[dPPP])[ktne] + (dist.f[dMMP])[ktsw] +
                 (dist.f[dPMP])[ktse] + (dist.f[dMPP])[ktnw] +
                 (dist.f[dPPM])[kbne] + (dist.f[dMMM])[kbsw] +
                 (dist.f[dPMM])[kbse] + (dist.f[dMPM])[kbnw]) -
                 rhoD[nodeIndex] - (vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1 + rhoD[nodeIndex])) * c1o2 + rhoD[nodeIndex] +
                 PRESS;
                 //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
            //////////////////////////////////////////////////////////////////////////
            mfcbb = (distAD.f[dP00])[nodeIndex   ];
            mfabb = (distAD.f[dM00])[kw  ];
            mfbcb = (distAD.f[d0P0])[nodeIndex   ];
            mfbab = (distAD.f[d0M0])[ks  ];
            mfbbc = (distAD.f[d00P])[nodeIndex   ];
            mfbba = (distAD.f[d00M])[kb  ];
            mfccb = (distAD.f[dPP0])[nodeIndex   ];
            mfaab = (distAD.f[dMM0])[ksw ];
            mfcab = (distAD.f[dPM0])[ks  ];
            mfacb = (distAD.f[dMP0])[kw  ];
            mfcbc = (distAD.f[dP0P])[nodeIndex   ];
            mfaba = (distAD.f[dM0M])[kbw ];
            mfcba = (distAD.f[dP0M])[kb  ];
            mfabc = (distAD.f[dM0P])[kw  ];
            mfbcc = (distAD.f[d0PP])[nodeIndex   ];
            mfbaa = (distAD.f[d0MM])[kbs ];
            mfbca = (distAD.f[d0PM])[kb  ];
            mfbac = (distAD.f[d0MP])[ks  ];
            mfbbb = (distAD.f[d000])[nodeIndex   ];
            mfccc = (distAD.f[dPPP])[nodeIndex   ];
            mfaac = (distAD.f[dMMP])[ksw ];
            mfcac = (distAD.f[dPMP])[ks  ];
            mfacc = (distAD.f[dMPP])[kw  ];
            mfcca = (distAD.f[dPPM])[kb  ];
            mfaaa = (distAD.f[dMMM])[kbsw];
            mfcaa = (distAD.f[dPMM])[kbs ];
            mfaca = (distAD.f[dMPM])[kbw ];
            //////////////////////////////////////////////////////////////////////////
            concD[nodeIndex] = 
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa)   + (mfaac + mfcca))) +
                 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba)   + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) +  mfbbb + CONC;
        }
    }
}

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

__global__ void LBResetMeanValuesSP27(
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

__global__ void LBResetMeanValuesAD27(
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
            rhoMP[kMac]= (dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                         (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                         (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                         (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                         (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                         (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                         (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                         (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                         (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                         (dist.f[d000])[kzero]+ 
                         (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                         (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                         (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                         (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vxMP[kMac] = (dist.f[dP00])[ke  ]- (dist.f[dM00])[kw  ]+ 
                         (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]+
                         (dist.f[dPM0])[kse ]- (dist.f[dMP0])[knw ]+
                         (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]+
                         (dist.f[dP0M])[kbe ]- (dist.f[dM0P])[ktw ]+
                         (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]+ 
                         (dist.f[dPMP])[ktse]- (dist.f[dMPP])[ktnw]+ 
                         (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]+ 
                         (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
           
            vyMP[kMac] = (dist.f[d0P0])[kn  ]- (dist.f[d0M0])[ks  ]+
                         (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]-
                         (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                         (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]+
                         (dist.f[d0PM])[kbn ]- (dist.f[d0MP])[kts ]+
                         (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]- 
                         (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                         (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                         (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vzMP[kMac] = (dist.f[d00P])[kt  ]- (dist.f[d00M])[kb  ]+
                         (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]-
                         (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                         (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]-
                         (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                         (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                         (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]- 
                         (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                         (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
        }
    }
}

__global__ void LBSetOutputWallVelocitySP27(real* vxD, real* vyD, real* vzD, real* vxWall, real* vyWall, real* vzWall,
                                            int numberOfWallNodes, int* kWallNodes, real* rhoD, real* pressD,
                                            unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY,
                                            unsigned int* neighborZ, unsigned long long numberOfLBnodes, real* DD,
                                            bool isEvenTimestep)
{
    const unsigned nodeIndex = getNodeIndex();

    if (nodeIndex < numberOfWallNodes) {
        unsigned int KWN = kWallNodes[nodeIndex];

        vxD[KWN] = 0.0; // vxWall[k];
        vyD[KWN] = 0.0; // vyWall[k];
        vzD[KWN] = 0.0; // vzWall[k];
    }
}

void CalcMacSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                 unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                 unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMacSP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                               numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalcMacSP27 execution failed");
}

void CalcMacCompSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                     unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                     unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMacCompSP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                                   numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalcMacCompSP27 execution failed");
}

void CalcMedSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                 unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                 unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMedSP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                               numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalcMedSP27 execution failed");
}

void CalcMedCompSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                     unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                     unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMedCompSP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                                   numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalcMedCompSP27 execution failed");
}

void CalcMedCompAD27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, real* concD, unsigned int* geoD,
                     unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                     unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD, real* DD_AD,
                     bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMedCompAD27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, concD, geoD, neighborX, neighborY, neighborZ,
                                                   numberOfLBnodes, DD, DD_AD, isEvenTimestep);
    getLastCudaError("LBCalcMedCompAD27 execution failed");
}

void CalcMacMedSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                    unsigned int* neighborY, unsigned int* neighborZ, unsigned int tdiff, unsigned long long numberOfLBnodes,
                    unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcMacMedSP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ, tdiff,
                                                  numberOfLBnodes, isEvenTimestep);
    getLastCudaError("LBCalcMacMedSP27 execution failed");
}

void ResetMeanValuesSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned long long numberOfLBnodes,
                         unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBResetMeanValuesSP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, numberOfLBnodes, isEvenTimestep);
    getLastCudaError("LBResetMeanValuesSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ResetMeanValuesAD27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, real* concD,
                         unsigned long long numberOfLBnodes, unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBResetMeanValuesAD27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, concD, numberOfLBnodes, isEvenTimestep);
    getLastCudaError("LBResetMeanValuesAD27 execution failed");
}

void LBCalcMeasurePoints27(real* vxMP, real* vyMP, real* vzMP, real* rhoMP, unsigned int* kMP,
                           unsigned int numberOfPointskMP, unsigned int MPClockCycle, unsigned int t, unsigned int* geoD,
                           unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                           unsigned long long numberOfLBnodes, real* DD, unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskMP);

    LBCalcMeasurePoints<<<grid.grid, grid.threads>>>(vxMP, vyMP, vzMP, rhoMP, kMP, numberOfPointskMP, MPClockCycle, t, geoD,
                                                     neighborX, neighborY, neighborZ, numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalcMeasurePoints execution failed");
}

void SetOutputWallVelocitySP27(unsigned int numberOfThreads, real* vxD, real* vyD, real* vzD, real* vxWall, real* vyWall,
                               real* vzWall, int numberOfWallNodes, int* kWallNodes, real* rhoD, real* pressD,
                               unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                               unsigned long long numberOfLBnodes, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfWallNodes);

    LBSetOutputWallVelocitySP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, vxWall, vyWall, vzWall, numberOfWallNodes,
                                                             kWallNodes, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                                             numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBSetOutputWallVelocitySP27 execution failed");
}

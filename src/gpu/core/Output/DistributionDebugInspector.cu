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
//! \file DistributionDebugInspector.cu
//! \ingroup Output
//! \author Henrik Asmuth, Henry Korb
//======================================================================================
#include "DistributionDebugInspector.h"

#include "Parameter/Parameter.h"
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "Utilities/KernelUtilities.h"

#include <cuda_helper/CudaGrid.h>
#include <cuda.h>

#include <iostream>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void printFs(
    real* distributions,
    bool isEvenTimestep,
    unsigned long long numberOfFluidNodes,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    uint* typeOfGridNode,
    real* coordX,
    real* coordY,
    real* coordZ,
    real minX,
    real maxX,
    real minY,
    real maxY,
    real minZ,
    real maxZ)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned k_000 = getNodeIndex();

    if (k_000 >= numberOfFluidNodes || typeOfGridNode[k_000]!=GEO_FLUID ) 
        return;

    real coordNodeX = coordX[k_000];
    real coordNodeY = coordY[k_000];
    real coordNodeZ = coordZ[k_000];

    if( coordNodeX>=minX && coordNodeX<=maxX &&
        coordNodeY>=minY && coordNodeY<=maxY &&
        coordNodeZ>=minZ && coordNodeZ<=maxZ    )
        {
            Distributions27 dist;
            getPointersToDistributions(dist, distributions, numberOfFluidNodes, isEvenTimestep);
            ////////////////////////////////////////////////////////////////////////////////
            //! - Set neighbor indices (necessary for indirect addressing)
            uint k_M00 = neighborX[k_000];
            uint k_0M0 = neighborY[k_000];
            uint k_00M = neighborZ[k_000];
            uint k_MM0 = neighborY[k_M00];
            uint k_M0M = neighborZ[k_M00];
            uint k_0MM = neighborZ[k_0M0];
            uint k_MMM = neighborZ[k_MM0];
            ////////////////////////////////////////////////////////////////////////////////////
            //! - Set local distributions
            //!
            real f_000 = (dist.f[d000])[k_000];
            real f_P00 = (dist.f[dP00])[k_000];
            real f_M00 = (dist.f[dM00])[k_M00];
            real f_0P0 = (dist.f[d0P0])[k_000];
            real f_0M0 = (dist.f[d0M0])[k_0M0];
            real f_00P = (dist.f[d00P])[k_000];
            real f_00M = (dist.f[d00M])[k_00M];
            real f_PP0 = (dist.f[dPP0])[k_000];
            real f_MM0 = (dist.f[dMM0])[k_MM0];
            real f_PM0 = (dist.f[dPM0])[k_0M0];
            real f_MP0 = (dist.f[dMP0])[k_M00];
            real f_P0P = (dist.f[dP0P])[k_000];
            real f_M0M = (dist.f[dM0M])[k_M0M];
            real f_P0M = (dist.f[dP0M])[k_00M];
            real f_M0P = (dist.f[dM0P])[k_M00];
            real f_0PP = (dist.f[d0PP])[k_000];
            real f_0MM = (dist.f[d0MM])[k_0MM];
            real f_0PM = (dist.f[d0PM])[k_00M];
            real f_0MP = (dist.f[d0MP])[k_0M0];
            real f_PPP = (dist.f[dPPP])[k_000];
            real f_MPP = (dist.f[dMPP])[k_M00];
            real f_PMP = (dist.f[dPMP])[k_0M0];
            real f_MMP = (dist.f[dMMP])[k_MM0];
            real f_PPM = (dist.f[dPPM])[k_00M];
            real f_MPM = (dist.f[dMPM])[k_M0M];
            real f_PMM = (dist.f[dPMM])[k_0MM];
            real f_MMM = (dist.f[dMMM])[k_MMM];

            real drho = ((((f_PPP + f_MMM) + (f_MPM + f_PMP)) + ((f_MPP + f_PMM) + (f_MMP + f_PPM))) +
                        (((f_0MP + f_0PM) + (f_0MM + f_0PP)) + ((f_M0P + f_P0M) + (f_M0M + f_P0P)) +
                        ((f_MP0 + f_PM0) + (f_MM0 + f_PP0))) +
                        ((f_M00 + f_P00) + (f_0M0 + f_0P0) + (f_00M + f_00P))) +
                            f_000;

            real oneOverRho = c1o1 / (c1o1 + drho);

            real vvx = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_PMM - f_MPP) + (f_PPM - f_MMP))) +
                        (((f_P0M - f_M0P) + (f_P0P - f_M0M)) + ((f_PM0 - f_MP0) + (f_PP0 - f_MM0))) + (f_P00 - f_M00)) *
                    oneOverRho;
            real vvy = ((((f_PPP - f_MMM) + (f_MPM - f_PMP)) + ((f_MPP - f_PMM) + (f_PPM - f_MMP))) +
                        (((f_0PM - f_0MP) + (f_0PP - f_0MM)) + ((f_MP0 - f_PM0) + (f_PP0 - f_MM0))) + (f_0P0 - f_0M0)) *
                    oneOverRho;
            real vvz = ((((f_PPP - f_MMM) + (f_PMP - f_MPM)) + ((f_MPP - f_PMM) + (f_MMP - f_PPM))) +
                        (((f_0MP - f_0PM) + (f_0PP - f_0MM)) + ((f_M0P - f_P0M) + (f_P0P - f_M0M))) + (f_00P - f_00M)) *
                    oneOverRho;

            printf("Node %u \t (%f\t%f\t%f)\n rho: %f\t velo: %f\t %f \t %f\n\n" , k_000, coordNodeX, coordNodeY, coordNodeZ, drho, vvx, vvy, vvz);
            printf("Node %u \t (%f\t%f\t%f)\n f_M00\t%f\t f_000\t%f\t f_P00\t%f\n f_MP0\t%f\t f_0P0\t%f\t f_PP0\t%f\n f_MM0\t%f\t f_0M0\t%f\t f_PM0\t%f\n f_M0P\t%f\t f_00P\t%f\t f_P0P\t%f\n f_M0M\t%f\t f_00M\t%f\t f_P0M\t%f\n f_MPP\t%f\t f_0PP\t%f\t f_PPP\t%f\n f_MPM\t%f\t f_0PM\t%f\t f_PPM\t%f\n f_MMP\t%f\t f_0MP\t%f\t f_PMP\t%f\n f_MMM\t%f\t f_0MM\t%f\t f_PMM\t%f\n\n\n" , k_000, coordNodeX, coordNodeY, coordNodeZ, f_M00, f_000, f_P00,f_MP0, f_0P0, f_PP0, f_MM0, f_0M0, f_PM0, f_M0P, f_00P, f_P0P, f_M0M, f_00M, f_P0M, f_MPP, f_0PP, f_PPP, f_MPM, f_0PM, f_PPM, f_MMP, f_0MP, f_PMP, f_MMM, f_0MM, f_PMM);

        }

}




void DistributionDebugInspector::inspect(std::shared_ptr<Parameter> para, uint level, uint t)
{
    if(this->inspectionLevel!=level)
        return;

    std::cout << tag << ": distributions on level " << level << " at t " << t <<  std::endl;

    vf::cuda::CudaGrid cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    printFs <<< cudaGrid.grid, cudaGrid.threads >>>(    para->getParD(level)->distributions.f[0],
                                                        para->getParD(level)->isEvenTimestep,
                                                        para->getParD(level)->numberOfNodes,
                                                        para->getParD(level)->neighborX,
                                                        para->getParD(level)->neighborY,
                                                        para->getParD(level)->neighborZ,
                                                        para->getParD(level)->typeOfGridNode,
                                                        para->getParD(level)->coordinateX,
                                                        para->getParD(level)->coordinateY,
                                                        para->getParD(level)->coordinateZ,
                                                        minX,
                                                        maxX,
                                                        minY,
                                                        maxY,
                                                        minZ,
                                                        maxZ);

}
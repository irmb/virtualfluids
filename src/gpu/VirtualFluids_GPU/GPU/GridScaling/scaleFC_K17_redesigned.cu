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
//! \file scaleFC_K17_redesigned.cu
//! \ingroup GPU/GridScaling
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================

#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"
#include "Kernel/Utilities/DistributionHelper.cuh"
#include "Kernel/ChimeraTransformation.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

__device__ __inline__ void calculateMomentsOnSourceNodes(
    Distributions27& distFine,
    real& omegaF,
    unsigned int& k_000,
    unsigned int& k_M00,
    unsigned int& k_0M0,
    unsigned int& k_00M,
    unsigned int& k_MM0,
    unsigned int& k_M0M,
    unsigned int& k_0MM,
    unsigned int& k_MMM,
    real& drho,
    real& velocityX,
    real& velocityY,
    real& velocityZ,
    real& kxyFromfcNEQ,
    real& kyzFromfcNEQ,
    real& kxzFromfcNEQ,
    real& kxxMyyFromfcNEQ,
    real& kxxMzzFromfcNEQ
    ){
        real f_000 = (distFine.f[DIR_000])[k_000]; 
        real f_P00 = (distFine.f[DIR_P00])[k_000];
        real f_M00 = (distFine.f[DIR_M00])[k_M00];
        real f_0P0 = (distFine.f[DIR_0P0])[k_000];
        real f_0M0 = (distFine.f[DIR_0M0])[k_0M0];
        real f_00P = (distFine.f[DIR_00P])[k_000];
        real f_00M = (distFine.f[DIR_00M])[k_00M];
        real f_PP0 = (distFine.f[DIR_PP0])[k_000];
        real f_MM0 = (distFine.f[DIR_MM0])[k_MM0];
        real f_PM0 = (distFine.f[DIR_PM0])[k_0M0];
        real f_MP0 = (distFine.f[DIR_MP0])[k_M00];
        real f_P0P = (distFine.f[DIR_P0P])[k_000];
        real f_M0M = (distFine.f[DIR_M0M])[k_M0M];
        real f_P0M = (distFine.f[DIR_P0M])[k_00M];
        real f_M0P = (distFine.f[DIR_M0P])[k_M00];
        real f_0PP = (distFine.f[DIR_0PP])[k_000];
        real f_0MM = (distFine.f[DIR_0MM])[k_0MM];
        real f_0PM = (distFine.f[DIR_0PM])[k_00M];
        real f_0MP = (distFine.f[DIR_0MP])[k_0M0];
        real f_PPP = (distFine.f[DIR_PPP])[k_000];
        real f_MPP = (distFine.f[DIR_MPP])[k_M00];
        real f_PMP = (distFine.f[DIR_PMP])[k_0M0];
        real f_MMP = (distFine.f[DIR_MMP])[k_MM0];
        real f_PPM = (distFine.f[DIR_PPM])[k_00M];
        real f_MPM = (distFine.f[DIR_MPM])[k_M0M];
        real f_PMM = (distFine.f[DIR_PMM])[k_0MM];
        real f_MMM = (distFine.f[DIR_MMM])[k_MMM];

        drho = f_P00 + f_M00 + f_0P0 + f_0M0 + f_00P + f_00M + f_PP0 + f_MM0 + f_PM0 + f_MP0 + f_P0P + f_M0M + f_P0M + f_M0P + f_0PP +
                   f_0MM + f_0PM + f_0MP + f_000 + f_PPP + f_MMP + f_PMP + f_MPP + f_PPM + f_MMM + f_PMM + f_MPM;
        velocityX = (((f_PPP - f_MMM) + (f_PMP - f_MPM) + (f_PPM - f_MMP) + (f_PMM - f_MPP)) +
                  (((f_PP0 - f_MM0) + (f_P0P - f_M0M)) + ((f_PM0 - f_MP0) + (f_P0M - f_M0P))) + (f_P00 - f_M00)) /
                  (c1o1 + drho);
        velocityY = (((f_PPP - f_MMM) + (f_MPP - f_PMM) + (f_PPM - f_MMP) + (f_MPM - f_PMP)) +
                  (((f_PP0 - f_MM0) + (f_0PP - f_0MM)) + ((f_0PM - f_0MP) + (f_MP0 - f_PM0))) + (f_0P0 - f_0M0)) /
                  (c1o1 + drho);
        velocityZ = (((f_PPP - f_MMM) + (f_MPP - f_PMM) + (f_PMP - f_MPM) + (f_MMP - f_PPM)) +
                  (((f_P0P - f_M0M) + (f_0PP - f_0MM)) + ((f_M0P - f_P0M) + (f_0MP - f_0PM))) + (f_00P - f_00M)) /
                  (c1o1 + drho);

        kxyFromfcNEQ =
            -c3o1 * omegaF *
            ((f_MM0 + f_MMM + f_MMP - f_MP0 - f_MPM - f_MPP - f_PM0 - f_PMM - f_PMP + f_PP0 + f_PPM + f_PPP) /
                 (c1o1 + drho) -
             ((velocityX * velocityY)));
        kyzFromfcNEQ =
            -c3o1 * omegaF *
            ((f_0MM + f_PMM + f_MMM - f_0MP - f_PMP - f_MMP - f_0PM - f_PPM - f_MPM + f_0PP + f_PPP + f_MPP) /
                 (c1o1 + drho) -
             ((velocityY * velocityZ)));
        kxzFromfcNEQ =
            -c3o1 * omegaF *
            ((f_M0M + f_MMM + f_MPM - f_M0P - f_MMP - f_MPP - f_P0M - f_PMM - f_PPM + f_P0P + f_PMP + f_PPP) /
                 (c1o1 + drho) -
             ((velocityX * velocityZ)));
        kxxMyyFromfcNEQ =
            -c3o2 * omegaF *
            ((f_M0M + f_M00 + f_M0P - f_0MM - f_0M0 - f_0MP - f_0PM - f_0P0 - f_0PP + f_P0M + f_P00 + f_P0P) / (c1o1 + drho) -
             ((velocityX * velocityX - velocityY * velocityY)));
        kxxMzzFromfcNEQ =
            -c3o2 * omegaF *
            ((f_MM0 + f_M00 + f_MP0 - f_0MM - f_0MP - f_00M - f_00P - f_0PM - f_0PP + f_PM0 + f_P00 + f_PP0) / (c1o1 + drho) -
             ((velocityX * velocityX - velocityZ * velocityZ)));
}

__global__ void scaleFC_K17_redesigned(
    real *distributionsCoarse,
    real *distributionsFine,
    unsigned int *neighborXcoarse,
    unsigned int *neighborYcoarse,
    unsigned int *neighborZcoarse,
    unsigned int *neighborXfine,
    unsigned int *neighborYfine,
    unsigned int *neighborZfine,
    unsigned int numberOfLBnodesCoarse,
    unsigned int numberOfLBnodesFine,
    bool isEvenTimestep,
    unsigned int *indicesCoarse000,
    unsigned int *indicesFineMMM,
    unsigned int numberOfInterfaceNodes,
    real omegaCoarse,
    real omegaFine,
    OffFC offsetFC)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get the thread index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned k_000 = vf::gpu::getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
    //! timestep is based on the esoteric twist algorithm \ref <a
    //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    Distributions27 distFine   = vf::gpu::getDistributionReferences27(distributionsFine,   numberOfLBnodesFine,   true);
    Distributions27 distCoarse = vf::gpu::getDistributionReferences27(distributionsCoarse, numberOfLBnodesCoarse, isEvenTimestep);

    // real *feC, *fwC, *fnC, *fsC, *ftC, *fbC, *fneC, *fswC, *fseC, *fnwC, *fteC, *fbwC, *fbeC, *ftwC, *ftnC, *fbsC,
    //     *fbnC, *ftsC, *fzeroC, *ftneC, *ftswC, *ftseC, *ftnwC, *fbneC, *fbswC, *fbseC, *fbnwC;

    // if (isEvenTimestep == true) {
    //     feC    = &distributionsCoarse[DIR_P00 * numberOfLBnodesCoarse];
    //     fwC    = &distributionsCoarse[DIR_M00 * numberOfLBnodesCoarse];
    //     fnC    = &distributionsCoarse[DIR_0P0 * numberOfLBnodesCoarse];
    //     fsC    = &distributionsCoarse[DIR_0M0 * numberOfLBnodesCoarse];
    //     ftC    = &distributionsCoarse[DIR_00P * numberOfLBnodesCoarse];
    //     fbC    = &distributionsCoarse[DIR_00M * numberOfLBnodesCoarse];
    //     fneC   = &distributionsCoarse[DIR_PP0 * numberOfLBnodesCoarse];
    //     fswC   = &distributionsCoarse[DIR_MM0 * numberOfLBnodesCoarse];
    //     fseC   = &distributionsCoarse[DIR_PM0 * numberOfLBnodesCoarse];
    //     fnwC   = &distributionsCoarse[DIR_MP0 * numberOfLBnodesCoarse];
    //     fteC   = &distributionsCoarse[DIR_P0P * numberOfLBnodesCoarse];
    //     fbwC   = &distributionsCoarse[DIR_M0M * numberOfLBnodesCoarse];
    //     fbeC   = &distributionsCoarse[DIR_P0M * numberOfLBnodesCoarse];
    //     ftwC   = &distributionsCoarse[DIR_M0P * numberOfLBnodesCoarse];
    //     ftnC   = &distributionsCoarse[DIR_0PP * numberOfLBnodesCoarse];
    //     fbsC   = &distributionsCoarse[DIR_0MM * numberOfLBnodesCoarse];
    //     fbnC   = &distributionsCoarse[DIR_0PM * numberOfLBnodesCoarse];
    //     ftsC   = &distributionsCoarse[DIR_0MP * numberOfLBnodesCoarse];
    //     fzeroC = &distributionsCoarse[DIR_000 * numberOfLBnodesCoarse];
    //     ftneC  = &distributionsCoarse[DIR_PPP * numberOfLBnodesCoarse];
    //     ftswC  = &distributionsCoarse[DIR_MMP * numberOfLBnodesCoarse];
    //     ftseC  = &distributionsCoarse[DIR_PMP * numberOfLBnodesCoarse];
    //     ftnwC  = &distributionsCoarse[DIR_MPP * numberOfLBnodesCoarse];
    //     fbneC  = &distributionsCoarse[DIR_PPM * numberOfLBnodesCoarse];
    //     fbswC  = &distributionsCoarse[DIR_MMM * numberOfLBnodesCoarse];
    //     fbseC  = &distributionsCoarse[DIR_PMM * numberOfLBnodesCoarse];
    //     fbnwC  = &distributionsCoarse[DIR_MPM * numberOfLBnodesCoarse];
    // } else {
    //     fwC    = &distributionsCoarse[DIR_P00 * numberOfLBnodesCoarse];
    //     feC    = &distributionsCoarse[DIR_M00 * numberOfLBnodesCoarse];
    //     fsC    = &distributionsCoarse[DIR_0P0 * numberOfLBnodesCoarse];
    //     fnC    = &distributionsCoarse[DIR_0M0 * numberOfLBnodesCoarse];
    //     fbC    = &distributionsCoarse[DIR_00P * numberOfLBnodesCoarse];
    //     ftC    = &distributionsCoarse[DIR_00M * numberOfLBnodesCoarse];
    //     fswC   = &distributionsCoarse[DIR_PP0 * numberOfLBnodesCoarse];
    //     fneC   = &distributionsCoarse[DIR_MM0 * numberOfLBnodesCoarse];
    //     fnwC   = &distributionsCoarse[DIR_PM0 * numberOfLBnodesCoarse];
    //     fseC   = &distributionsCoarse[DIR_MP0 * numberOfLBnodesCoarse];
    //     fbwC   = &distributionsCoarse[DIR_P0P * numberOfLBnodesCoarse];
    //     fteC   = &distributionsCoarse[DIR_M0M * numberOfLBnodesCoarse];
    //     ftwC   = &distributionsCoarse[DIR_P0M * numberOfLBnodesCoarse];
    //     fbeC   = &distributionsCoarse[DIR_M0P * numberOfLBnodesCoarse];
    //     fbsC   = &distributionsCoarse[DIR_0PP * numberOfLBnodesCoarse];
    //     ftnC   = &distributionsCoarse[DIR_0MM * numberOfLBnodesCoarse];
    //     ftsC   = &distributionsCoarse[DIR_0PM * numberOfLBnodesCoarse];
    //     fbnC   = &distributionsCoarse[DIR_0MP * numberOfLBnodesCoarse];
    //     fzeroC = &distributionsCoarse[DIR_000 * numberOfLBnodesCoarse];
    //     fbswC  = &distributionsCoarse[DIR_PPP * numberOfLBnodesCoarse];
    //     fbneC  = &distributionsCoarse[DIR_MMP * numberOfLBnodesCoarse];
    //     fbnwC  = &distributionsCoarse[DIR_PMP * numberOfLBnodesCoarse];
    //     fbseC  = &distributionsCoarse[DIR_MPP * numberOfLBnodesCoarse];
    //     ftswC  = &distributionsCoarse[DIR_PPM * numberOfLBnodesCoarse];
    //     ftneC  = &distributionsCoarse[DIR_MMM * numberOfLBnodesCoarse];
    //     ftnwC  = &distributionsCoarse[DIR_PMM * numberOfLBnodesCoarse];
    //     ftseC  = &distributionsCoarse[DIR_MPM * numberOfLBnodesCoarse];
    // }

    ////////////////////////////////////////////////////////////////////////////////
    real eps_new = c2o1;
    real omegaF  = omegaFine;
    real omegaC  = omegaCoarse;

    real xoff, yoff, zoff;
    real xoff_sq, yoff_sq, zoff_sq;

    real press;
    real drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP;
    real drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP;
    real drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP;
    real drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP;
    real drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM;
    real drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM;
    real drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM;
    real drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM;

    real f_P00, f_M00, f_0P0, f_0M0, f_00P, f_00M, f_PP0, f_MM0, f_PM0, f_MP0, f_P0P, f_M0M, f_P0M, f_M0P, f_0PP, f_0MM, f_0PM, f_0MP, f_000,
        f_PPP, f_MMP, f_PMP, f_MPP, f_PPM, f_MMM, f_PMM, f_MPM;

    real kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP;
    real kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP;
    real kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP;
    real kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP;
    real kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM;
    real kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM;
    real kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM;
    real kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM;

    real a_000, a_100, a_010, a_001, a_200, a_020, a_002, a_110, a_101, a_011;
    real b_000, b_100, b_010, b_001, b_200, b_020, b_002, b_110, b_101, b_011;
    real c_000, c_100, c_010, c_001, c_200, c_020, c_002, c_110, c_101, c_011;
    real d_000, d_100, d_010, d_001, d_110, d_101, d_011;

    if (k_000 < numberOfInterfaceNodes) {
        //////////////////////////////////////////////////////////////////////////
        xoff    = offsetFC.xOffFC[k_000];
        yoff    = offsetFC.yOffFC[k_000];
        zoff    = offsetFC.zOffFC[k_000];
        xoff_sq = xoff * xoff;
        yoff_sq = yoff * yoff;
        zoff_sq = zoff * zoff;

        //////////////////////////////////////////////////////////////////////////
        // source node BSW = MMM
        //////////////////////////////////////////////////////////////////////////
        // index of the base node
        unsigned int k_base_000 = indicesFineMMM[k_000];
        unsigned int k_base_M00 = neighborXfine [k_base_000];
        unsigned int k_base_0M0 = neighborYfine [k_base_000];
        unsigned int k_base_00M = neighborZfine [k_base_000];
        unsigned int k_base_MM0 = neighborYfine [k_base_M00];
        unsigned int k_base_M0M = neighborZfine [k_base_M00];
        unsigned int k_base_0MM = neighborZfine [k_base_0M0];
        unsigned int k_base_MMM = neighborZfine [k_base_MM0];
        //////////////////////////////////////////////////////////////////////////
        // index
        unsigned int k_000 = k_base_000;
        unsigned int k_M00 = k_base_M00;
        unsigned int k_0M0 = k_base_0M0;
        unsigned int k_00M = k_base_00M;
        unsigned int k_MM0 = k_base_MM0;
        unsigned int k_M0M = k_base_M0M;
        unsigned int k_0MM = k_base_0MM;
        unsigned int k_MMM = k_base_MMM;

        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM,
            kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM);

        //////////////////////////////////////////////////////////////////////////
        // source node TSW = MMP
        //////////////////////////////////////////////////////////////////////////
        // index
        k_000 = k_00M;
        k_M00 = k_M0M;
        k_0M0 = k_0MM;
        k_00M = neighborZfine[k_00M];
        k_MM0 = k_MMM;
        k_M0M = neighborZfine[k_M0M];
        k_0MM = neighborZfine[k_0MM];
        k_MMM = neighborZfine[k_MMM];

        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP,
            kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP);

        //////////////////////////////////////////////////////////////////////////
        // source node TSE = PMP
        //////////////////////////////////////////////////////////////////////////
        // index
        k_000 = k_M00;
        k_M00 = neighborXfine[k_M00];
        k_0M0 = k_MM0;
        k_00M = k_M0M;
        k_MM0 = neighborXfine[k_MM0];
        k_M0M = neighborXfine[k_M0M];
        k_0MM = k_MMM;
        k_MMM = neighborXfine[k_MMM];

        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP,
            kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP);

        //////////////////////////////////////////////////////////////////////////
        // source node BSE = PMM 
        //////////////////////////////////////////////////////////////////////////
        // index
        k_00M = k_000;
        k_M0M = k_M00;
        k_0MM = k_0M0;
        k_MMM = k_MM0;
        k_000 = k_base_M00;
        k_M00 = neighborXfine[k_base_M00];
        k_0M0 = k_base_MM0;
        k_MM0 = neighborXfine[k_base_MM0];

        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM,
            kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM);

        //////////////////////////////////////////////////////////////////////////
        // source node BNW = MPM
        //////////////////////////////////////////////////////////////////////////
        // index of the base node
        k_base_000 = k_base_0M0;
        k_base_M00 = k_base_MM0;
        k_base_0M0 = neighborYfine[k_base_0M0];
        k_base_00M = k_base_0MM;
        k_base_MM0 = neighborYfine[k_base_MM0];
        k_base_M0M = k_base_MMM;
        k_base_0MM = neighborYfine[k_base_0MM];
        k_base_MMM = neighborYfine[k_base_MMM];
        //////////////////////////////////////////////////////////////////////////
        // index
        k_000 = k_base_000;
        k_M00 = k_base_M00;
        k_0M0 = k_base_0M0;
        k_00M = k_base_00M;
        k_MM0 = k_base_MM0;
        k_M0M = k_base_M0M;
        k_0MM = k_base_0MM;
        k_MMM = k_base_MMM;

        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM,
            kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM);

        //////////////////////////////////////////////////////////////////////////
        // source node TNW = MPP
        //////////////////////////////////////////////////////////////////////////
        // index
        k_000 = k_00M;
        k_M00 = k_M0M;
        k_0M0 = k_0MM;
        k_00M = neighborZfine[k_00M];
        k_MM0 = k_MMM;
        k_M0M = neighborZfine[k_M0M];
        k_0MM = neighborZfine[k_0MM];
        k_MMM = neighborZfine[k_MMM];
        
        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP,
            kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP);

        //////////////////////////////////////////////////////////////////////////
        // source node TNE = PPP
        //////////////////////////////////////////////////////////////////////////
        // index
        k_000 = k_M00;
        k_M00 = neighborXfine[k_M00];
        k_0M0 = k_MM0;
        k_00M = k_M0M;
        k_MM0 = neighborXfine[k_MM0];
        k_M0M = neighborXfine[k_M0M];
        k_0MM = k_MMM;
        k_MMM = neighborXfine[k_MMM];

        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP,
            kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP);

        //////////////////////////////////////////////////////////////////////////
        // source node BNE = PPM
        //////////////////////////////////////////////////////////////////////////
        // index
        k_00M = k_000;
        k_M0M = k_M00;
        k_0MM = k_0M0;
        k_MMM = k_MM0;
        k_000 = k_base_M00;
        k_M00 = neighborXfine[k_base_M00];
        k_0M0 = k_base_MM0;
        k_MM0 = neighborXfine[k_base_MM0];
        
        calculateMomentsOnSourceNodes( distFine, omegaF,
            k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM,
            kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM);

        //////////////////////////////////////////////////////////////////////////
        // 3. scaling coefficients
        //////////////////////////////////////////////////////////////////////////
        a_000 = (-kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_PPP + kxxMyyFromfcNEQ_MPM + kxxMyyFromfcNEQ_MPP -
                kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_PMP + kxxMyyFromfcNEQ_MMM + kxxMyyFromfcNEQ_MMP -
                kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_PPP + kxxMzzFromfcNEQ_MPM + kxxMzzFromfcNEQ_MPP -
                kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_PMP + kxxMzzFromfcNEQ_MMM + kxxMzzFromfcNEQ_MMP -
                c2o1 * kxyFromfcNEQ_PPM - c2o1 * kxyFromfcNEQ_PPP - c2o1 * kxyFromfcNEQ_MPM - c2o1 * kxyFromfcNEQ_MPP +
                c2o1 * kxyFromfcNEQ_PMM + c2o1 * kxyFromfcNEQ_PMP + c2o1 * kxyFromfcNEQ_MMM + c2o1 * kxyFromfcNEQ_MMP +
                c2o1 * kxzFromfcNEQ_PPM - c2o1 * kxzFromfcNEQ_PPP + c2o1 * kxzFromfcNEQ_MPM - c2o1 * kxzFromfcNEQ_MPP +
                c2o1 * kxzFromfcNEQ_PMM - c2o1 * kxzFromfcNEQ_PMP + c2o1 * kxzFromfcNEQ_MMM - c2o1 * kxzFromfcNEQ_MMP +
                c8o1 * vx1_PPM + c8o1 * vx1_PPP + c8o1 * vx1_MPM + c8o1 * vx1_MPP + c8o1 * vx1_PMM + c8o1 * vx1_PMP +
                c8o1 * vx1_MMM + c8o1 * vx1_MMP + c2o1 * vx2_PPM + c2o1 * vx2_PPP - c2o1 * vx2_MPM - c2o1 * vx2_MPP -
                c2o1 * vx2_PMM - c2o1 * vx2_PMP + c2o1 * vx2_MMM + c2o1 * vx2_MMP - c2o1 * vx3_PPM + c2o1 * vx3_PPP +
                c2o1 * vx3_MPM - c2o1 * vx3_MPP - c2o1 * vx3_PMM + c2o1 * vx3_PMP + c2o1 * vx3_MMM - c2o1 * vx3_MMP) /
                c64o1;
        b_000 = (c2o1 * kxxMyyFromfcNEQ_PPM + c2o1 * kxxMyyFromfcNEQ_PPP + c2o1 * kxxMyyFromfcNEQ_MPM +
                c2o1 * kxxMyyFromfcNEQ_MPP - c2o1 * kxxMyyFromfcNEQ_PMM - c2o1 * kxxMyyFromfcNEQ_PMP -
                c2o1 * kxxMyyFromfcNEQ_MMM - c2o1 * kxxMyyFromfcNEQ_MMP - kxxMzzFromfcNEQ_PPM - kxxMzzFromfcNEQ_PPP -
                kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_MPP + kxxMzzFromfcNEQ_PMM + kxxMzzFromfcNEQ_PMP +
                kxxMzzFromfcNEQ_MMM + kxxMzzFromfcNEQ_MMP - c2o1 * kxyFromfcNEQ_PPM - c2o1 * kxyFromfcNEQ_PPP +
                c2o1 * kxyFromfcNEQ_MPM + c2o1 * kxyFromfcNEQ_MPP - c2o1 * kxyFromfcNEQ_PMM - c2o1 * kxyFromfcNEQ_PMP +
                c2o1 * kxyFromfcNEQ_MMM + c2o1 * kxyFromfcNEQ_MMP + c2o1 * kyzFromfcNEQ_PPM - c2o1 * kyzFromfcNEQ_PPP +
                c2o1 * kyzFromfcNEQ_MPM - c2o1 * kyzFromfcNEQ_MPP + c2o1 * kyzFromfcNEQ_PMM - c2o1 * kyzFromfcNEQ_PMP +
                c2o1 * kyzFromfcNEQ_MMM - c2o1 * kyzFromfcNEQ_MMP + c2o1 * vx1_PPM + c2o1 * vx1_PPP - c2o1 * vx1_MPM -
                c2o1 * vx1_MPP - c2o1 * vx1_PMM - c2o1 * vx1_PMP + c2o1 * vx1_MMM + c2o1 * vx1_MMP + c8o1 * vx2_PPM +
                c8o1 * vx2_PPP + c8o1 * vx2_MPM + c8o1 * vx2_MPP + c8o1 * vx2_PMM + c8o1 * vx2_PMP + c8o1 * vx2_MMM +
                c8o1 * vx2_MMP - c2o1 * vx3_PPM + c2o1 * vx3_PPP - c2o1 * vx3_MPM + c2o1 * vx3_MPP + c2o1 * vx3_PMM -
                c2o1 * vx3_PMP + c2o1 * vx3_MMM - c2o1 * vx3_MMP) /
                c64o1;
        c_000 = (kxxMyyFromfcNEQ_PPM - kxxMyyFromfcNEQ_PPP + kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_MPP +
                kxxMyyFromfcNEQ_PMM - kxxMyyFromfcNEQ_PMP + kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_MMP -
                c2o1 * kxxMzzFromfcNEQ_PPM + c2o1 * kxxMzzFromfcNEQ_PPP - c2o1 * kxxMzzFromfcNEQ_MPM +
                c2o1 * kxxMzzFromfcNEQ_MPP - c2o1 * kxxMzzFromfcNEQ_PMM + c2o1 * kxxMzzFromfcNEQ_PMP -
                c2o1 * kxxMzzFromfcNEQ_MMM + c2o1 * kxxMzzFromfcNEQ_MMP - c2o1 * kxzFromfcNEQ_PPM -
                c2o1 * kxzFromfcNEQ_PPP + c2o1 * kxzFromfcNEQ_MPM + c2o1 * kxzFromfcNEQ_MPP - c2o1 * kxzFromfcNEQ_PMM -
                c2o1 * kxzFromfcNEQ_PMP + c2o1 * kxzFromfcNEQ_MMM + c2o1 * kxzFromfcNEQ_MMP - c2o1 * kyzFromfcNEQ_PPM -
                c2o1 * kyzFromfcNEQ_PPP - c2o1 * kyzFromfcNEQ_MPM - c2o1 * kyzFromfcNEQ_MPP + c2o1 * kyzFromfcNEQ_PMM +
                c2o1 * kyzFromfcNEQ_PMP + c2o1 * kyzFromfcNEQ_MMM + c2o1 * kyzFromfcNEQ_MMP - c2o1 * vx1_PPM +
                c2o1 * vx1_PPP + c2o1 * vx1_MPM - c2o1 * vx1_MPP - c2o1 * vx1_PMM + c2o1 * vx1_PMP + c2o1 * vx1_MMM -
                c2o1 * vx1_MMP - c2o1 * vx2_PPM + c2o1 * vx2_PPP - c2o1 * vx2_MPM + c2o1 * vx2_MPP + c2o1 * vx2_PMM -
                c2o1 * vx2_PMP + c2o1 * vx2_MMM - c2o1 * vx2_MMP + c8o1 * vx3_PPM + c8o1 * vx3_PPP + c8o1 * vx3_MPM +
                c8o1 * vx3_MPP + c8o1 * vx3_PMM + c8o1 * vx3_PMP + c8o1 * vx3_MMM + c8o1 * vx3_MMP) /
                c64o1;
        a_100  = (vx1_PPM + vx1_PPP - vx1_MPM - vx1_MPP + vx1_PMM + vx1_PMP - vx1_MMM - vx1_MMP) / c4o1;
        b_100  = (vx2_PPM + vx2_PPP - vx2_MPM - vx2_MPP + vx2_PMM + vx2_PMP - vx2_MMM - vx2_MMP) / c4o1;
        c_100  = (vx3_PPM + vx3_PPP - vx3_MPM - vx3_MPP + vx3_PMM + vx3_PMP - vx3_MMM - vx3_MMP) / c4o1;
        a_200 = (kxxMyyFromfcNEQ_PPM + kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MPM - kxxMyyFromfcNEQ_MPP +
                kxxMyyFromfcNEQ_PMM + kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MMM - kxxMyyFromfcNEQ_MMP +
                kxxMzzFromfcNEQ_PPM + kxxMzzFromfcNEQ_PPP - kxxMzzFromfcNEQ_MPM - kxxMzzFromfcNEQ_MPP +
                kxxMzzFromfcNEQ_PMM + kxxMzzFromfcNEQ_PMP - kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_MMP + c2o1 * vx2_PPM +
                c2o1 * vx2_PPP - c2o1 * vx2_MPM - c2o1 * vx2_MPP - c2o1 * vx2_PMM - c2o1 * vx2_PMP + c2o1 * vx2_MMM +
                c2o1 * vx2_MMP - c2o1 * vx3_PPM + c2o1 * vx3_PPP + c2o1 * vx3_MPM - c2o1 * vx3_MPP - c2o1 * vx3_PMM +
                c2o1 * vx3_PMP + c2o1 * vx3_MMM - c2o1 * vx3_MMP) /
                c16o1;
        b_200 = (kxyFromfcNEQ_PPM + kxyFromfcNEQ_PPP - kxyFromfcNEQ_MPM - kxyFromfcNEQ_MPP + kxyFromfcNEQ_PMM +
                kxyFromfcNEQ_PMP - kxyFromfcNEQ_MMM - kxyFromfcNEQ_MMP - c2o1 * vx1_PPM - c2o1 * vx1_PPP +
                c2o1 * vx1_MPM + c2o1 * vx1_MPP + c2o1 * vx1_PMM + c2o1 * vx1_PMP - c2o1 * vx1_MMM - c2o1 * vx1_MMP) /
                c8o1;
        c_200 = (kxzFromfcNEQ_PPM + kxzFromfcNEQ_PPP - kxzFromfcNEQ_MPM - kxzFromfcNEQ_MPP + kxzFromfcNEQ_PMM +
                kxzFromfcNEQ_PMP - kxzFromfcNEQ_MMM - kxzFromfcNEQ_MMP + c2o1 * vx1_PPM - c2o1 * vx1_PPP -
                c2o1 * vx1_MPM + c2o1 * vx1_MPP + c2o1 * vx1_PMM - c2o1 * vx1_PMP - c2o1 * vx1_MMM + c2o1 * vx1_MMP) /
                c8o1;
        a_010  = (vx1_PPM + vx1_PPP + vx1_MPM + vx1_MPP - vx1_PMM - vx1_PMP - vx1_MMM - vx1_MMP) / c4o1;
        b_010  = (vx2_PPM + vx2_PPP + vx2_MPM + vx2_MPP - vx2_PMM - vx2_PMP - vx2_MMM - vx2_MMP) / c4o1;
        c_010  = (vx3_PPM + vx3_PPP + vx3_MPM + vx3_MPP - vx3_PMM - vx3_PMP - vx3_MMM - vx3_MMP) / c4o1;
        a_020 = (kxyFromfcNEQ_PPM + kxyFromfcNEQ_PPP + kxyFromfcNEQ_MPM + kxyFromfcNEQ_MPP - kxyFromfcNEQ_PMM -
                kxyFromfcNEQ_PMP - kxyFromfcNEQ_MMM - kxyFromfcNEQ_MMP - c2o1 * vx2_PPM - c2o1 * vx2_PPP +
                c2o1 * vx2_MPM + c2o1 * vx2_MPP + c2o1 * vx2_PMM + c2o1 * vx2_PMP - c2o1 * vx2_MMM - c2o1 * vx2_MMP) /
                c8o1;
        b_020 = (-c2o1 * kxxMyyFromfcNEQ_PPM - c2o1 * kxxMyyFromfcNEQ_PPP - c2o1 * kxxMyyFromfcNEQ_MPM -
                c2o1 * kxxMyyFromfcNEQ_MPP + c2o1 * kxxMyyFromfcNEQ_PMM + c2o1 * kxxMyyFromfcNEQ_PMP +
                c2o1 * kxxMyyFromfcNEQ_MMM + c2o1 * kxxMyyFromfcNEQ_MMP + kxxMzzFromfcNEQ_PPM + kxxMzzFromfcNEQ_PPP +
                kxxMzzFromfcNEQ_MPM + kxxMzzFromfcNEQ_MPP - kxxMzzFromfcNEQ_PMM - kxxMzzFromfcNEQ_PMP -
                kxxMzzFromfcNEQ_MMM - kxxMzzFromfcNEQ_MMP + c2o1 * vx1_PPM + c2o1 * vx1_PPP - c2o1 * vx1_MPM -
                c2o1 * vx1_MPP - c2o1 * vx1_PMM - c2o1 * vx1_PMP + c2o1 * vx1_MMM + c2o1 * vx1_MMP - c2o1 * vx3_PPM +
                c2o1 * vx3_PPP - c2o1 * vx3_MPM + c2o1 * vx3_MPP + c2o1 * vx3_PMM - c2o1 * vx3_PMP + c2o1 * vx3_MMM -
                c2o1 * vx3_MMP) /
                c16o1;
        c_020 = (kyzFromfcNEQ_PPM + kyzFromfcNEQ_PPP + kyzFromfcNEQ_MPM + kyzFromfcNEQ_MPP - kyzFromfcNEQ_PMM -
                kyzFromfcNEQ_PMP - kyzFromfcNEQ_MMM - kyzFromfcNEQ_MMP + c2o1 * vx2_PPM - c2o1 * vx2_PPP +
                c2o1 * vx2_MPM - c2o1 * vx2_MPP - c2o1 * vx2_PMM + c2o1 * vx2_PMP - c2o1 * vx2_MMM + c2o1 * vx2_MMP) /
                c8o1;
        a_001  = (-vx1_PPM + vx1_PPP - vx1_MPM + vx1_MPP - vx1_PMM + vx1_PMP - vx1_MMM + vx1_MMP) / c4o1;
        b_001  = (-vx2_PPM + vx2_PPP - vx2_MPM + vx2_MPP - vx2_PMM + vx2_PMP - vx2_MMM + vx2_MMP) / c4o1;
        c_001  = (-vx3_PPM + vx3_PPP - vx3_MPM + vx3_MPP - vx3_PMM + vx3_PMP - vx3_MMM + vx3_MMP) / c4o1;
        a_002 = (-kxzFromfcNEQ_PPM + kxzFromfcNEQ_PPP - kxzFromfcNEQ_MPM + kxzFromfcNEQ_MPP - kxzFromfcNEQ_PMM +
                kxzFromfcNEQ_PMP - kxzFromfcNEQ_MMM + kxzFromfcNEQ_MMP + c2o1 * vx3_PPM - c2o1 * vx3_PPP -
                c2o1 * vx3_MPM + c2o1 * vx3_MPP + c2o1 * vx3_PMM - c2o1 * vx3_PMP - c2o1 * vx3_MMM + c2o1 * vx3_MMP) /
                c8o1;
        b_002 = (-kyzFromfcNEQ_PPM + kyzFromfcNEQ_PPP - kyzFromfcNEQ_MPM + kyzFromfcNEQ_MPP - kyzFromfcNEQ_PMM +
                kyzFromfcNEQ_PMP - kyzFromfcNEQ_MMM + kyzFromfcNEQ_MMP + c2o1 * vx3_PPM - c2o1 * vx3_PPP +
                c2o1 * vx3_MPM - c2o1 * vx3_MPP - c2o1 * vx3_PMM + c2o1 * vx3_PMP - c2o1 * vx3_MMM + c2o1 * vx3_MMP) /
                c8o1;
        c_002 = (-kxxMyyFromfcNEQ_PPM + kxxMyyFromfcNEQ_PPP - kxxMyyFromfcNEQ_MPM + kxxMyyFromfcNEQ_MPP -
                kxxMyyFromfcNEQ_PMM + kxxMyyFromfcNEQ_PMP - kxxMyyFromfcNEQ_MMM + kxxMyyFromfcNEQ_MMP +
                c2o1 * kxxMzzFromfcNEQ_PPM - c2o1 * kxxMzzFromfcNEQ_PPP + c2o1 * kxxMzzFromfcNEQ_MPM -
                c2o1 * kxxMzzFromfcNEQ_MPP + c2o1 * kxxMzzFromfcNEQ_PMM - c2o1 * kxxMzzFromfcNEQ_PMP +
                c2o1 * kxxMzzFromfcNEQ_MMM - c2o1 * kxxMzzFromfcNEQ_MMP - c2o1 * vx1_PPM + c2o1 * vx1_PPP +
                c2o1 * vx1_MPM - c2o1 * vx1_MPP - c2o1 * vx1_PMM + c2o1 * vx1_PMP + c2o1 * vx1_MMM - c2o1 * vx1_MMP -
                c2o1 * vx2_PPM + c2o1 * vx2_PPP - c2o1 * vx2_MPM + c2o1 * vx2_MPP + c2o1 * vx2_PMM - c2o1 * vx2_PMP +
                c2o1 * vx2_MMM - c2o1 * vx2_MMP) /
                c16o1;
        a_110 = (vx1_PPM + vx1_PPP - vx1_MPM - vx1_MPP - vx1_PMM - vx1_PMP + vx1_MMM + vx1_MMP) / c2o1;
        b_110 = (vx2_PPM + vx2_PPP - vx2_MPM - vx2_MPP - vx2_PMM - vx2_PMP + vx2_MMM + vx2_MMP) / c2o1;
        c_110 = (vx3_PPM + vx3_PPP - vx3_MPM - vx3_MPP - vx3_PMM - vx3_PMP + vx3_MMM + vx3_MMP) / c2o1;
        a_101 = (-vx1_PPM + vx1_PPP + vx1_MPM - vx1_MPP - vx1_PMM + vx1_PMP + vx1_MMM - vx1_MMP) / c2o1;
        b_101 = (-vx2_PPM + vx2_PPP + vx2_MPM - vx2_MPP - vx2_PMM + vx2_PMP + vx2_MMM - vx2_MMP) / c2o1;
        c_101 = (-vx3_PPM + vx3_PPP + vx3_MPM - vx3_MPP - vx3_PMM + vx3_PMP + vx3_MMM - vx3_MMP) / c2o1;
        a_011 = (-vx1_PPM + vx1_PPP - vx1_MPM + vx1_MPP + vx1_PMM - vx1_PMP + vx1_MMM - vx1_MMP) / c2o1;
        b_011 = (-vx2_PPM + vx2_PPP - vx2_MPM + vx2_MPP + vx2_PMM - vx2_PMP + vx2_MMM - vx2_MMP) / c2o1;
        c_011 = (-vx3_PPM + vx3_PPP - vx3_MPM + vx3_MPP + vx3_PMM - vx3_PMP + vx3_MMM - vx3_MMP) / c2o1;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        real kxyAverage    = c0o1;
        real kyzAverage    = c0o1;
        real kxzAverage    = c0o1;
        real kxxMyyAverage = c0o1;
        real kxxMzzAverage = c0o1;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // drho
        real LaplaceRho = ((xoff != c0o1) || (yoff != c0o1) || (zoff != c0o1))
                          ? c0o1
                          : -c3o1 * (a_100 * a_100 + b_010 * b_010 + c_001 * c_001) - c6o1 * (b_100 * a_010 + c_100 * a_001 + c_010 * b_001);
        d_000 = ( drho_PPM + drho_PPP + drho_MPM + drho_MPP + drho_PMM + drho_PMP + drho_MMM + drho_MMP - c2o1 * LaplaceRho) * c1o8;
        d_100 = ( drho_PPM + drho_PPP - drho_MPM - drho_MPP + drho_PMM + drho_PMP - drho_MMM - drho_MMP) * c1o4;
        d_010 = ( drho_PPM + drho_PPP + drho_MPM + drho_MPP - drho_PMM - drho_PMP - drho_MMM - drho_MMP) * c1o4;
        d_001 = (-drho_PPM + drho_PPP - drho_MPM + drho_MPP - drho_PMM + drho_PMP - drho_MMM + drho_MMP) * c1o4;
        d_110 = ( drho_PPM + drho_PPP - drho_MPM - drho_MPP - drho_PMM - drho_PMP + drho_MMM + drho_MMP) * c1o2;
        d_101 = (-drho_PPM + drho_PPP + drho_MPM - drho_MPP - drho_PMM + drho_PMP + drho_MMM - drho_MMP) * c1o2;
        d_011 = (-drho_PPM + drho_PPP - drho_MPM + drho_MPP + drho_PMM - drho_PMP + drho_MMM - drho_MMP) * c1o2;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // x------x
        // |      |
        // |   ---+--->X
        // |      |  \
        // x------x   \
        //          offset-vector
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        a_000 = a_000 + xoff * a_100 + yoff * a_010 + zoff * a_001 + xoff_sq * a_200 + yoff_sq * a_020 + zoff_sq * a_002 +
                xoff * yoff * a_110 + xoff * zoff * a_101 + yoff * zoff * a_011;
        a_100 = a_100 + c2o1 * xoff * a_200 + yoff * a_110 + zoff * a_101;
        a_010 = a_010 + c2o1 * yoff * a_020 + xoff * a_110 + zoff * a_011;
        a_001 = a_001 + c2o1 * zoff * a_002 + xoff * a_101 + yoff * a_011;
        b_000 = b_000 + xoff * b_100 + yoff * b_010 + zoff * b_001 + xoff_sq * b_200 + yoff_sq * b_020 + zoff_sq * b_002 +
                xoff * yoff * b_110 + xoff * zoff * b_101 + yoff * zoff * b_011;
        b_100 = b_100 + c2o1 * xoff * b_200 + yoff * b_110 + zoff * b_101;
        b_010 = b_010 + c2o1 * yoff * b_020 + xoff * b_110 + zoff * b_011;
        b_001 = b_001 + c2o1 * zoff * b_002 + xoff * b_101 + yoff * b_011;
        c_000 = c_000 + xoff * c_100 + yoff * c_010 + zoff * c_001 + xoff_sq * c_200 + yoff_sq * c_020 + zoff_sq * c_002 +
                xoff * yoff * c_110 + xoff * zoff * c_101 + yoff * zoff * c_011;
        c_100 = c_100 + c2o1 * xoff * c_200 + yoff * c_110 + zoff * c_101;
        c_010 = c_010 + c2o1 * yoff * c_020 + xoff * c_110 + zoff * c_011;
        c_001 = c_001 + c2o1 * zoff * c_002 + xoff * c_101 + yoff * c_011;
        d_000 = d_000 + xoff * d_100 + yoff * d_010 + zoff * d_001 + xoff * yoff * d_110 + xoff * zoff * d_101 + yoff * zoff * d_011;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        real m_111 = c0o1;
        real m_211 = c0o1;
        real m_011 = c0o1;
        real m_121 = c0o1;
        real m_101 = c0o1;
        real m_112 = c0o1;
        real m_110 = c0o1;
        real m_221 = c0o1;
        real m_001 = c0o1;
        real m_201 = c0o1;
        real m_021 = c0o1;
        real m_212 = c0o1;
        real m_010 = c0o1;
        real m_210 = c0o1;
        real m_012 = c0o1;
        real m_122 = c0o1;
        real m_100 = c0o1;
        real m_120 = c0o1;
        real m_102 = c0o1;
        real m_222 = c0o1;
        real m_022 = c0o1;
        real m_202 = c0o1;
        real m_002 = c0o1;
        real m_220 = c0o1;
        real m_020 = c0o1;
        real m_200 = c0o1;
        real m_000 = c0o1;

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Define aliases to use the same variable for the distributions (f's):
        //!
        real& f_000 = m_111;
        real& f_P00 = m_211;
        real& f_M00 = m_011;
        real& f_0P0 = m_121;
        real& f_0M0 = m_101;
        real& f_00P = m_112;
        real& f_00M = m_110;
        real& f_PP0 = m_221;
        real& f_MM0 = m_001;
        real& f_PM0 = m_201;
        real& f_MP0 = m_021;
        real& f_P0P = m_212;
        real& f_M0M = m_010;
        real& f_P0M = m_210;
        real& f_M0P = m_012;
        real& f_0PP = m_122;
        real& f_0MM = m_100;
        real& f_0PM = m_120;
        real& f_0MP = m_102;
        real& f_PPP = m_222;
        real& f_MPP = m_022;
        real& f_PMP = m_202;
        real& f_MMP = m_002;
        real& f_PPM = m_220;
        real& f_MPM = m_020;
        real& f_PMM = m_200;
        real& f_MMM = m_000;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        real m0, m1, m2, vvx, vvy, vvz, vx2, vy2, vz2, oMdrho;
        real mxxPyyPzz, mxxMyy, mxxMzz, mxxyPyzz, mxxyMyzz, mxxzPyyz, mxxzMyyz, mxyyPxzz, mxyyMxzz;
        real NeqOn = c1o1; // zero;//one;   //.... one = on ..... zero = off
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // Position Coarse 0., 0., 0.
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // x = 0.;
        // y = 0.;
        // z = 0.;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        press = d_000;
        vvx   = a_000;
        vvy   = b_000;
        vvz   = c_000;

        m_000 = press; // if drho is interpolated directly

        vx2    = vvx * vvx;
        vy2    = vvy * vvy;
        vz2    = vvz * vvz;
        oMdrho = c1o1;

        // linear combinations for second order moments
        mxxPyyPzz = m_000;

        mxxMyy = -c2o3 * ((a_100 - b_010) + kxxMyyAverage) * eps_new / omegaC * (c1o1 + press);
        mxxMzz = -c2o3 * ((a_100 - c_001) + kxxMzzAverage) * eps_new / omegaC * (c1o1 + press);

        m_011 = -c1o3 * ((b_001 + c_010) + kyzAverage) * eps_new / omegaC * (c1o1 + press);
        m_101 = -c1o3 * ((a_001 + c_100) + kxzAverage) * eps_new / omegaC * (c1o1 + press);
        m_110 = -c1o3 * ((a_010 + b_100) + kxyAverage) * eps_new / omegaC * (c1o1 + press);

        // linear combinations back
        m_200 = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz) * NeqOn;
        m_020 = c1o3 * (-c2o1 * mxxMyy + mxxMzz + mxxPyyPzz) * NeqOn;
        m_002 = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * NeqOn;

        // linear combinations for third order moments
        m_111 = c0o1;

        mxxyPyzz = c0o1;
        mxxyMyzz = c0o1;
        mxxzPyyz = c0o1;
        mxxzMyyz = c0o1;
        mxyyPxzz = c0o1;
        mxyyMxzz = c0o1;

        // linear combinations back
        m_210 = (mxxyMyzz + mxxyPyzz) * c1o2;
        m_012 = (-mxxyMyzz + mxxyPyzz) * c1o2;
        m_201 = (mxxzMyyz + mxxzPyyz) * c1o2;
        m_021 = (-mxxzMyyz + mxxzPyyz) * c1o2;
        m_120 = (mxyyMxzz + mxyyPxzz) * c1o2;
        m_102 = (-mxyyMxzz + mxyyPxzz) * c1o2;

        // fourth order moments
        m_022 = m_000 * c1o9;
        m_202 = m_022;
        m_220 = m_022;

        // fifth order moments

        // sixth order moments
        m_222 = m_000 * c1o27;

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
        //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
        //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
        //! ]</b></a>
        //!
        ////////////////////////////////////////////////////////////////////////////////////
        // X - Dir
        backwardInverseChimeraWithK(m_000, m_100, m_200, vvx, vx2, c1o1, c1o1);
        backwardChimera(            m_010, m_110, m_210, vvx, vx2);
        backwardInverseChimeraWithK(m_020, m_120, m_220, vvx, vx2, c3o1, c1o3);
        backwardChimera(            m_001, m_101, m_201, vvx, vx2);
        backwardChimera(            m_011, m_111, m_211, vvx, vx2);
        backwardChimera(            m_021, m_121, m_221, vvx, vx2);
        backwardInverseChimeraWithK(m_002, m_102, m_202, vvx, vx2, c3o1, c1o3);
        backwardChimera(            m_012, m_112, m_212, vvx, vx2);
        backwardInverseChimeraWithK(m_022, m_122, m_222, vvx, vx2, c9o1, c1o9);

        ////////////////////////////////////////////////////////////////////////////////////
        // Y - Dir
        backwardInverseChimeraWithK(m_000, m_010, m_020, vvy, vy2, c6o1, c1o6);
        backwardChimera(            m_001, m_011, m_021, vvy, vy2);
        backwardInverseChimeraWithK(m_002, m_012, m_022, vvy, vy2, c18o1, c1o18);
        backwardInverseChimeraWithK(m_100, m_110, m_120, vvy, vy2, c3o2, c2o3);
        backwardChimera(            m_101, m_111, m_121, vvy, vy2);
        backwardInverseChimeraWithK(m_102, m_112, m_122, vvy, vy2, c9o2, c2o9);
        backwardInverseChimeraWithK(m_200, m_210, m_220, vvy, vy2, c6o1, c1o6);
        backwardChimera(            m_201, m_211, m_221, vvy, vy2);
        backwardInverseChimeraWithK(m_202, m_212, m_222, vvy, vy2, c18o1, c1o18);

        ////////////////////////////////////////////////////////////////////////////////////
        // Z - Dir
        backwardInverseChimeraWithK(m_000, m_001, m_002, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(m_010, m_011, m_012, vvz, vz2, c9o1,  c1o9);
        backwardInverseChimeraWithK(m_020, m_021, m_022, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(m_100, m_101, m_102, vvz, vz2, c9o1,  c1o9);
        backwardInverseChimeraWithK(m_110, m_111, m_112, vvz, vz2, c9o4,  c4o9);
        backwardInverseChimeraWithK(m_120, m_121, m_122, vvz, vz2, c9o1,  c1o9);
        backwardInverseChimeraWithK(m_200, m_201, m_202, vvz, vz2, c36o1, c1o36);
        backwardInverseChimeraWithK(m_210, m_211, m_212, vvz, vz2, c9o1,  c1o9);
        backwardInverseChimeraWithK(m_220, m_221, m_222, vvz, vz2, c36o1, c1o36);


        ////////////////////////////////////////////////////////////////////////////////////
        // index 0
        k_000 = indicesCoarse000[k_000];
        k_M00 = neighborXcoarse[k_000];
        k_0M0 = neighborYcoarse[k_000];
        k_00M = neighborZcoarse[k_000];
        k_MM0 = neighborYcoarse[k_M00];
        k_M0M = neighborZcoarse[k_M00];
        k_0MM = neighborZcoarse[k_0M0];
        k_MMM = neighborZcoarse[k_MM0];
        ////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Write distributions: style of reading and writing the distributions from/to
        //! stored arrays dependent on timestep is based on the esoteric twist algorithm
        //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        (distCoarse.f[DIR_000])[k_000] = f_000;
        (distCoarse.f[DIR_P00])[k_000] = f_P00;
        (distCoarse.f[DIR_M00])[k_M00] = f_M00;
        (distCoarse.f[DIR_0P0])[k_000] = f_0P0;
        (distCoarse.f[DIR_0M0])[k_0M0] = f_0M0;
        (distCoarse.f[DIR_00P])[k_000] = f_00P;
        (distCoarse.f[DIR_00M])[k_00M] = f_00M;
        (distCoarse.f[DIR_PP0])[k_000] = f_PP0;
        (distCoarse.f[DIR_MM0])[k_MM0] = f_MM0;
        (distCoarse.f[DIR_PM0])[k_0M0] = f_PM0;
        (distCoarse.f[DIR_MP0])[k_M00] = f_MP0;
        (distCoarse.f[DIR_P0P])[k_000] = f_P0P;
        (distCoarse.f[DIR_M0M])[k_M0M] = f_M0M;
        (distCoarse.f[DIR_P0M])[k_00M] = f_P0M;
        (distCoarse.f[DIR_M0P])[k_M00] = f_M0P;
        (distCoarse.f[DIR_0PP])[k_000] = f_0PP;
        (distCoarse.f[DIR_0MM])[k_0MM] = f_0MM;
        (distCoarse.f[DIR_0PM])[k_00M] = f_0PM;
        (distCoarse.f[DIR_0MP])[k_0M0] = f_0MP;
        (distCoarse.f[DIR_PPP])[k_000] = f_PPP;
        (distCoarse.f[DIR_MPP])[k_M00] = f_MPP;
        (distCoarse.f[DIR_PMP])[k_0M0] = f_PMP;
        (distCoarse.f[DIR_MMP])[k_MM0] = f_MMP;
        (distCoarse.f[DIR_PPM])[k_00M] = f_PPM;
        (distCoarse.f[DIR_MPM])[k_M0M] = f_MPM;
        (distCoarse.f[DIR_PMM])[k_0MM] = f_PMM;
        (distCoarse.f[DIR_MMM])[k_MMM] = f_MMM;

        ////////////////////////////////////////////////////////////////////////////////////
    }
}

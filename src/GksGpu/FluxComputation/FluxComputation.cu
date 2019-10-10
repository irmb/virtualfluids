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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file FluxComputation.cu
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#include "FluxComputation.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"

#include "CellProperties/CellProperties.cuh"

#include "FluxComputation/Moments.cuh"
#include "FluxComputation/Reconstruction.cuh"
#include "FluxComputation/Transformation.cuh"
#include "FluxComputation/ExpansionCoefficients.cuh"
#include "FluxComputation/AssembleFlux.cuh"
#include "FluxComputation/ApplyFlux.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//! \brief This is a CUDA Kernel that computes the face index and calls \ref fluxFunction for this index
__global__                 void fluxKernel  ( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint numberOfEntities );

//! \brief This function performs the flux computation
//!
//! The \ref reconstructFiniteDifferences function computes the flow state 
//! on the cell face as well as the gradient of the flow state.
//! The flow state is computed by directly averaging the conserved variables
//! of the positive and negative neighboring cells and afterwards transformed
//! to primitive variables.
//! The normal derivatives are also clearly defined by positive and negative
//! cells and are computed by central difference.
//! For the tangential derivative the orientation (in the sense of normal to 
//! which coordinate direction) of the face is utilized to choose the correct 
//! data access pattern.
//! While the positive and negative cell indices are stored per face, the cells
//! for the tangential derivatives are obtained pointer chasing, i.e. utilizing
//! the cell to cell connectivity stored per cell. The gradients of the 
//! conserved variables are divided by the density for implementation reasons.
//! 
//! In order to unify the flux computation for all three possible directions,
//! the conserved variables in both gradients and flow state are rotated into
//! a local frame of reference with the function \ref transformGlobalToLocal.
//! 
//! Subsequently the expansion coefficients according to Appendix C in 
//! <a href="https://doi.org/10.1142/9324"><b>[ Kun Xu, (2015), DOI: 10.1142/9324 ]</b></a>
//! are computed by \ref computeExpansionCoefficients.
//! 
//! The next major stage in the flux computation is related to the time 
//! derivative, which is evaluated from the spatial derivatives. First a set
//! of several moments of the equilibrium distribution is evaluated explicitly.
//! The formulas for the moments can be found for example in 
//! <a href="https://doi.org/10.1142/9324"><b>[ Kun Xu (2015), DOI: 10.1142/9324 ]</b></a>. 
//! The moments are subsequently weighted by the spatial expansion coefficients
//! and assembled to form the time derivative according to Eq. (17) in 
//! <a href="https://doi.org/10.1016/j.ijthermalsci.2018.10.004"><b>[ Lenz et al. (2019), DOI: 10.1016/j.ijthermalsci.2018.10.004 ]</b></a>
//! in function \ref computeTimeDerivative. The temporal expansion coefficients 
//! are then computed with the same function used for the spatial expansion 
//! coefficients.
//! 
//! With all expansion coefficients in place, they are used to assemble the 
//! flux according to Eq. (19) in 
//! <a href="https://doi.org/10.1016/j.ijthermalsci.2018.10.004"><b>[ Lenz et al. (2019), DOI: 10.1016/j.ijthermalsci.2018.10.004 ]</b></a>
//! in the function \ref assembleFlux. This function includes the Prandtl 
//! number fix from Eq. (27) in 
//! <a href="https://doi.org/10.1006/jcph.2001.6790"><b>[ Kun Xu (2001), DOI: 10.1006/jcph.2001.6790 ]</b></a>.
//! 
//! At this stage the flux exists as conserved variables in the local frame of 
//! reference. Transformation to the global frame follows in the function 
//! \ref transformLocalToGlobal.
//! 
//! Concluding the flux computation the fluxes are applied to the positive and
//! negative cells \ref DataBase::dataUpdate variables. In this process several 
//! \ref CellProperties, which are stored as bitmaps, are evaluated, such as 
//! impenetrable walls. Finally, the functions \ref applyFluxToNegCell and 
//! \ref applyFluxToPosCell add the fluxes to \ref DataBase::dataUpdate of 
//! negative and positive cells respectively. As there is the possibility of race
//! conditions in the case, where multiple faces are processed simultaneously.
//! In order to prevent this, CUDA atomics are used to protect the write operations.
//! 
//! The fluxes are multiplied by the face area \f$\Delta x^2\f$. Hence, they have to 
//! interpreted as absolute amounts of conserved quantities, in opposition to the
//! storage of cell average conserved quantities.
__host__ __device__ inline void fluxFunction( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, uint level, bool evaluateCommFaces )
{
    CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfInnerFaces, 64, CudaUtility::computeStream);

    runKernel(fluxKernel,
                fluxFunction,
                dataBase->getDeviceType(), grid,
                dataBase->toStruct(),
                parameters,
                'x',
                dataBase->perLevelCount[level].startOfFacesX);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void fluxKernel(DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    fluxFunction( dataBase, parameters, direction, startIndex, index );
}

__host__ __device__ inline void fluxFunction(DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint index)
{
    uint faceIndex = startIndex + index;

    real K = parameters.K;

    direction = dataBase.faceOrientation[ faceIndex ];

    parameters.D1 = parameters.D;
    parameters.D2 = parameters.D;

    //////////////////////////////////////////////////////////////////////////

    PrimitiveVariables facePrim;

    //////////////////////////////////////////////////////////////////////////

    real ax[LENGTH_CELL_DATA];
    real ay[LENGTH_CELL_DATA];
    real az[LENGTH_CELL_DATA];
    real at[LENGTH_CELL_DATA];

#pragma unroll
    for( uint i = 0; i < LENGTH_CELL_DATA; i++ )
    { 
        ax[i] = c0o1; 
        ay[i] = c0o1; 
        az[i] = c0o1; 
        at[i] = c0o1;
    }

    //////////////////////////////////////////////////////////////////////////

    {
        ConservedVariables gradN, gradT1, gradT2;

        reconstructFiniteDifferences(faceIndex,
                                     dataBase,
                                     parameters,
                                     direction,
                                     gradN,
                                     gradT1,
                                     gradT2,
                                     facePrim,
                                     K);

        transformGlobalToLocal( gradN , direction );
        transformGlobalToLocal( gradT1, direction );
        transformGlobalToLocal( gradT2, direction );

        transformGlobalToLocal( facePrim, direction );

        computeExpansionCoefficients(facePrim, gradN , K, ax);
        computeExpansionCoefficients(facePrim, gradT1, K, ay);
        computeExpansionCoefficients(facePrim, gradT2, K, az);

        //////////////////////////////////////////////////////////////////////////
    }

    {
        ConservedVariables flux;

        {
            real momentU [ NUMBER_OF_MOMENTS ]; 
            real momentV [ NUMBER_OF_MOMENTS ]; 
            real momentW [ NUMBER_OF_MOMENTS ]; 
            real momentXi[ NUMBER_OF_MOMENTS ];

            computeMoments( facePrim, K, momentU, momentV, momentW, momentXi );

            Vec3 force = parameters.force;

            transformGlobalToLocal(force, direction);

            {
                ConservedVariables timeGrad;
                computeTimeDerivative( facePrim, 
                                       momentU, 
                                       momentV, 
                                       momentW, 
                                       momentXi, 
                                       ax, ay, az,
                                       force,
                                       timeGrad );

                computeExpansionCoefficients( facePrim, timeGrad, K, at );
            }
            {
                real timeCoefficients[4];
                computeTimeCoefficients( facePrim, parameters, timeCoefficients );

                real heatFlux;
                assembleFlux( facePrim, 
                              momentU, momentV, momentW, momentXi,
                              ax, ay, az, at, 
                              timeCoefficients, 
                              parameters,
                              force,
                              flux,
                              heatFlux );

                transformLocalToGlobal( flux, direction );
            }
        }

        //////////////////////////////////////////////////////////////////////////

        {
            uint negCellIdx = dataBase.faceToCell[ NEG_CELL(faceIndex, dataBase.numberOfFaces) ];
            uint posCellIdx = dataBase.faceToCell[ POS_CELL(faceIndex, dataBase.numberOfFaces) ];

        #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
            atomicAdd( &( dataBase.diffusivity[ negCellIdx ] ), parameters.D * parameters.dx * parameters.dx * parameters.dt );
            atomicAdd( &( dataBase.diffusivity[ posCellIdx ] ), parameters.D * parameters.dx * parameters.dx * parameters.dt );
        #endif

            CellProperties negCellProperties = dataBase.cellProperties[ negCellIdx ];
            CellProperties posCellProperties = dataBase.cellProperties[ posCellIdx ];

            //if( isCellProperties( negCellProperties, CELL_PROPERTIES_IS_FLUX_BC ) || 
            //    isCellProperties( posCellProperties, CELL_PROPERTIES_IS_FLUX_BC ) )
            //    return;

            if( isCellProperties( negCellProperties, CELL_PROPERTIES_WALL ) || 
                isCellProperties( posCellProperties, CELL_PROPERTIES_WALL ) )
            {
                flux.rho    = c0o1;
            #ifdef USE_PASSIVE_SCALAR
                flux.rhoS_1 = c0o1;
                flux.rhoS_2 = c0o1;
            #endif //USE_PASSIVE_SCALAR
            }

            if( isCellProperties( negCellProperties, CELL_PROPERTIES_IS_INSULATED ) || 
                isCellProperties( posCellProperties, CELL_PROPERTIES_IS_INSULATED ) )
            {
                flux.rhoE   = c0o1;
            }

            uint negCellParentIdx = dataBase.parentCell[ negCellIdx ];
            uint posCellParentIdx = dataBase.parentCell[ posCellIdx ];

            //if( !( negCellParentIdx != INVALID_INDEX ) != !( posCellParentIdx != INVALID_INDEX ) ) // XOR
            if( ( negCellParentIdx == INVALID_INDEX ) != ( posCellParentIdx == INVALID_INDEX ) ) // XOR
            {
                if( !isCellProperties( negCellProperties, CELL_PROPERTIES_GHOST ) && 
                    !isCellProperties( posCellProperties, CELL_PROPERTIES_GHOST ) )
                {
                    if (negCellParentIdx != INVALID_INDEX)
                    {
                        applyFluxToNegCell(dataBase, negCellParentIdx, flux, direction, parameters);
                    }

                    if (posCellParentIdx != INVALID_INDEX)
                    {
                        applyFluxToPosCell(dataBase, posCellParentIdx, flux, direction, parameters);
                    }
                }
            }

            applyFluxToNegCell(dataBase, negCellIdx, flux, direction, parameters);
            applyFluxToPosCell(dataBase, posCellIdx, flux, direction, parameters);
        }
    }
}

#include "CellUpdate.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>

#include "PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/ThermalDependencies.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu {

inline __host__ __device__ real getTurbulentViscosityDeardorff(const DataBaseStruct& dataBase, const Parameters& parameters, const uint cellIndex, const ConservedVariables& cons )
{
    // See FDS 6 Technical Reference Guide, Section 4.2.3

    PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

    ConservedVariables neighborCons;
    PrimitiveVariables neighborPrim;

    real kSGS = c0o1;

    {
        real uHead = c1o2 * prim.U;

        {
            uint neighborCellIndex = dataBase.cellToCell[CELL_TO_CELL(cellIndex, 0, dataBase.numberOfCells)];
            readCellData(cellIndex, dataBase, neighborCons);
            neighborPrim = toPrimitiveVariables(neighborCons, parameters.K);

            uHead += c1o4 * neighborPrim.U;
        }
        {
            uint neighborCellIndex = dataBase.cellToCell[CELL_TO_CELL(cellIndex, 1, dataBase.numberOfCells)];
            readCellData(cellIndex, dataBase, neighborCons);
            neighborPrim = toPrimitiveVariables(neighborCons, parameters.K);

            uHead += c1o4 * neighborPrim.U;
        }

        kSGS += c1o2 * ( prim.U - uHead ) * ( prim.U - uHead );
    }

    {
        real vHead = c1o2 * prim.V;

        {
            uint neighborCellIndex = dataBase.cellToCell[CELL_TO_CELL(cellIndex, 2, dataBase.numberOfCells)];
            readCellData(cellIndex, dataBase, neighborCons);
            neighborPrim = toPrimitiveVariables(neighborCons, parameters.K);

            vHead += c1o4 * neighborPrim.V;
        }
        {
            uint neighborCellIndex = dataBase.cellToCell[CELL_TO_CELL(cellIndex, 3, dataBase.numberOfCells)];
            readCellData(cellIndex, dataBase, neighborCons);
            neighborPrim = toPrimitiveVariables(neighborCons, parameters.K);

            vHead += c1o4 * neighborPrim.V;
        }

        kSGS += c1o2 * ( prim.V - vHead ) * ( prim.V - vHead );
    }

    {
        real wHead = c1o2 * prim.W;

        {
            uint neighborCellIndex = dataBase.cellToCell[CELL_TO_CELL(cellIndex, 4, dataBase.numberOfCells)];
            readCellData(cellIndex, dataBase, neighborCons);
            neighborPrim = toPrimitiveVariables(neighborCons, parameters.K);

            wHead += c1o4 * neighborPrim.W;
        }
        {
            uint neighborCellIndex = dataBase.cellToCell[CELL_TO_CELL(cellIndex, 5, dataBase.numberOfCells)];
            readCellData(cellIndex, dataBase, neighborCons);
            neighborPrim = toPrimitiveVariables(neighborCons, parameters.K);

            wHead += c1o4 * neighborPrim.W;
        }

        kSGS += c1o2 * ( prim.W - wHead ) * ( prim.W - wHead );
    }

    //real turbulentViscosity = prim.rho * parameters.dx * c1o10 * sqrt(kSGS) / 0.3;

    dataBase.diffusivity[cellIndex] = (realAccumulator) kSGS;

    //printf("%f", kSGS);

    return kSGS;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void chemicalReaction(DataBaseStruct dataBase, Parameters parameters, uint cellIndex, ConservedVariables& cons)
{
    // see FDS 5 Technical reference guide, section 6.1.4 for combustion model
#ifdef USE_PASSIVE_SCALAR
    if (parameters.enableReaction)
    {
        CellProperties cellProperties = dataBase.cellProperties[ cellIndex ];

        if( isCellProperties( cellProperties, CELL_PROPERTIES_FINE_GHOST ) ) return;

        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

        //////////////////////////////////////////////////////////////////////////

        //real diffusivity = getTurbulentViscosityDeardorff(dataBase, parameters, cellIndex, cons);
        //real diffusivity = dataBase.diffusivity[ cellIndex ];
        real diffusivity = dataBase.diffusivity[ cellIndex ] / ( c6o1 * parameters.dx * parameters.dx * parameters.dt );
        dataBase.diffusivity[ cellIndex ] = c0o1;

        //////////////////////////////////////////////////////////////////////////

        real mixingTimeScale = real(0.1) * parameters.dx * parameters.dx / diffusivity;

        //real kSGS = getTurbulentViscosityDeardorff(dataBase, parameters, cellIndex, cons);

        //real mixingTimeScale_d = parameters.dx * parameters.dx / parameters.D;

        //real mixingTimeScale_u = real(0.4) * parameters.dx / sqrt( c2o3 * kSGS );

        //real mixingTimeScale_g = sqrt( c2o1 * parameters.dx / fabs( parameters.force.z ) );

        //real mixingTimeScale = fminf( mixingTimeScale_d, mixingTimeScale_u );
        //mixingTimeScale      = fminf( mixingTimeScale_g, mixingTimeScale   );

        //////////////////////////////////////////////////////////////////////////

        real Y_F = prim.S_1;
        real Y_P = prim.S_2;

        real Y_A = c1o1 - Y_F - Y_P;

        ///////////////////////////////////////////////////////////////////////////////

        real Y_O2 = rX * ( M_O2 / M_A ) * Y_A;

        ///////////////////////////////////////////////////////////////////////////////

        real s = M_F / ( c2o1 * M_O2 );

        real heatReleaseRate = cons.rho * fminf(Y_F, s * Y_O2) / mixingTimeScale * ( parameters.heatOfReaction / M_F );

        //////////////////////////////////////////////////////////////////////////

        if( heatReleaseRate < c0o1 )
            heatReleaseRate = c0o1;

        //////////////////////////////////////////////////////////////////////////

        if( parameters.useHeatReleaseRateLimiter )
        if( heatReleaseRate > parameters.heatReleaseRateLimiter )
            heatReleaseRate = parameters.heatReleaseRateLimiter;

        //////////////////////////////////////////////////////////////////////////

        real drhoY_F = heatReleaseRate * parameters.dt / ( parameters.heatOfReaction / M_F );

        //real r = c1o1 + ( c1o2 / rX ) * ( M_A / M_F );
        real r = c1o1 + ( c2o1 / rX ) * ( M_A / M_F );

        cons.rhoS_1 -=     drhoY_F;
        cons.rhoS_2 += r * drhoY_F;
        cons.rhoE   += heatReleaseRate * parameters.dt;

        //////////////////////////////////////////////////////////////////////////
    }

#endif // USE_PASSIVE_SCALAR
}

} // namespace GksGpu
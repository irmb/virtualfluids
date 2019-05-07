#include "CellUpdate.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/ThermalDependencies.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void chemicalReaction(DataBaseStruct dataBase, Parameters parameters, uint cellIndex, ConservedVariables& cons)
{

#ifdef USE_PASSIVE_SCALAR
    if (true)
    {
        CellProperties cellProperties = dataBase.cellProperties[ cellIndex ];

        if( !isCellProperties( cellProperties, CELL_PROPERTIES_FINE_GHOST ) )
        //if( !isCellProperties( cellProperties, CELL_PROPERTIES_GHOST ) )
        {
            PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

            //////////////////////////////////////////////////////////////////////////

            real Y_F = prim.S_1;
            real Y_P = prim.S_2;

            real Y_A = one - Y_F - Y_P;

            real M = one / ( Y_A / M_A
                           + Y_F / M_F
                           + Y_P / M_P );

            real X_A = Y_A * M / M_A;
            real X_F = Y_F * M / M_F;
            real X_P = Y_P * M / M_P;

            ///////////////////////////////////////////////////////////////////////////////

            real X_O2 = real(0.21) * X_A;

            ///////////////////////////////////////////////////////////////////////////////

            {
                //const real heatOfReaction = real(802310.0); // kJ / kmol
                const real heatOfReaction = real(8000.0); // kJ / kmol
                //const real heatOfReaction = two * real(1000.0); // kJ / kmol

                //////////////////////////////////////////////////////////////////////////

                PrimitiveVariables limitPrim = prim;

                real r = 1.001;

                limitPrim.lambda /= r;

                ConservedVariables limitCons = toConservedVariables(limitPrim, parameters.K);

                real maxHeatRelease = limitCons.rhoE - cons.rhoE;

                real dX_F_max = maxHeatRelease * M / cons.rho / heatOfReaction;

                //////////////////////////////////////////////////////////////////////////

                real dX_F = fminf(X_F, c1o2 * X_O2);

                //////////////////////////////////////////////////////////////////////////

                dX_F = fminf(dX_F_max, dX_F);

                if( dX_F < zero ) dX_F = zero;

                //////////////////////////////////////////////////////////////////////////

                real dn_F = cons.rho * dX_F / M;

                real releasedHeat = dn_F * heatOfReaction;

                ///////////////////////////////////////////////////////////////////////////////

                //real X_F_new = X_F - dX_F;
                //real X_P_new = X_P + dX_F;

                real X_A_new = X_A - two * dX_F / real(0.21);
                real X_F_new = X_F - dX_F;

                real X_P_new = one - X_A_new - X_F_new;

                real Z1 = X_F_new * M_F / M;
                real Z2 = X_P_new * M_P / M;

                //////////////////////////////////////////////////////////////////////////

                //if( Z1 < zero ) { Z2 -= Z1; Z1 = zero; }
                //if( Z2 < zero ) { Z1 -= Z2; Z2 = zero; }

                //if( Z1 > one  ) { Z2 += Z1 - one; Z1 = one; }
                //if( Z2 > one  ) { Z1 += Z2 - one; Z2 = one; }

                if( Z1 < zero ) Z1 = zero;
                if( Z2 < zero ) Z2 = zero;

                if( Z1 > one  ) Z1 = one;
                if( Z2 > one  ) Z2 = one;

                if( Z1 + Z2 > one )
                {
                    real faktor = (Z1 + Z2);

                    Z1 /= faktor;
                    Z2 /= faktor;
                }

                ///////////////////////////////////////////////////////////////////////////////

                ConservedVariables testCons = cons;

                testCons.rhoE += releasedHeat;

                PrimitiveVariables testPrim = toPrimitiveVariables(testCons, parameters.K);

                //////////////////////////////////////////////////////////////////////////

                //if( getT( testPrim ) < 20 )
                {
                    cons.rhoS_1 = Z1 * cons.rho;
                    cons.rhoS_2 = Z2 * cons.rho;
                    cons.rhoE += releasedHeat;
                }
            }
        }

        if( cons.rhoS_1 < zero ) cons.rhoS_1 = zero;
        if( cons.rhoS_2 < zero ) cons.rhoS_2 = zero;

        if( cons.rhoS_1 > cons.rho  ) cons.rhoS_1 = cons.rho;
        if( cons.rhoS_2 > cons.rho  ) cons.rhoS_2 = cons.rho;

        if( cons.rhoS_1 + cons.rhoS_2 > cons.rho )
        {
            real faktor = (cons.rhoS_1 + cons.rhoS_2) / cons.rho;

            cons.rhoS_1 /= faktor;
            cons.rhoS_2 /= faktor;
        }
    }

#endif // USE_PASSIVE_SCALAR
}

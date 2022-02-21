#include "PlanarAverageProbe.h"

#include <cuda/CudaGrid.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"

void PlanarAverageProbe::findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                            std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                            std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                            int level)
{
    real dx = abs(para->getParH(level)->coordX_SP[1]-para->getParH(level)->coordX_SP[para->getParH(level)->neighborX_SP[1]]);
    
    real *pointCoordsInplane1_par, *pointCoordsInplane2_par, *pointCoordsNormal_par;
    std::vector<real> *pointCoordsInplane1, *pointCoordsInplane2, *pointCoordsNormal;
    
    if(this->planeNormal == 'x'){  
                                    pointCoordsNormal       = &pointCoordsX_level; 
                                    pointCoordsInplane1     = &pointCoordsY_level; 
                                    pointCoordsInplane2     = &pointCoordsZ_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordX_SP; 
                                    pointCoordsInplane1_par = para->getParH(level)->coordY_SP; 
                                    pointCoordsInplane2_par = para->getParH(level)->coordZ_SP;
                                }
    if(this->planeNormal == 'y'){  
                                    pointCoordsNormal       = &pointCoordsY_level; 
                                    pointCoordsInplane1     = &pointCoordsX_level; 
                                    pointCoordsInplane2     = &pointCoordsZ_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordY_SP; 
                                    pointCoordsInplane1_par = para->getParH(level)->coordX_SP; 
                                    pointCoordsInplane2_par = para->getParH(level)->coordZ_SP;
                                }
    if(this->planeNormal == 'z'){  
                                    pointCoordsNormal       = &pointCoordsZ_level; 
                                    pointCoordsInplane1     = &pointCoordsX_level; 
                                    pointCoordsInplane2     = &pointCoordsY_level;
                                    pointCoordsNormal_par   = para->getParH(level)->coordZ_SP; 
                                    pointCoordsInplane1_par = para->getParH(level)->coordX_SP; 
                                    pointCoordsInplane2_par = para->getParH(level)->coordY_SP;
                                }

    // Find all points along the normal direction and add to normalCoords
    for(uint j=1; j<para->getParH(level)->size_Mat_SP; j++ )
    {
        if(para->getParH(level)->geoSP[j] == GEO_FLUID)
        {   
            if( std::find(pointCoordsNormal->begin(), pointCoordsNormal->end(), pointCoordsNormal_par[j]) == pointCoordsNormal->end())  
            {
                pointCoordsNormal->push_back( pointCoordsNormal_par[j] );
                pointCoordsInplane1->push_back(999999.);
                pointCoordsInplane2->push_back(999999.);
            }
        }
    }
    std::sort(pointCoordsNormal->begin(), pointCoordsNormal->end());
    
    // Find all pointCoords in the first plane 
    for(uint j=1; j<para->getParH(level)->size_Mat_SP; j++ )
    {
        if( para->getParH(level)->geoSP[j] == GEO_FLUID && pointCoordsNormal_par[j] == pointCoordsNormal->at(0)) 
        {
            // pointCoordsNormal->push_back( pointCoordsNormal_par[j] ); //not needed in current state, might become relevant for two-point correlations
            // pointCoordsInplane1->push_back( pointCoordsInplane1_par[j] );
            // pointCoordsInplane2->push_back( pointCoordsInplane2_par[j] );
            probeIndices_level.push_back(j);
        }
    }
}

void PlanarAverageProbe::calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, int level)
{


}

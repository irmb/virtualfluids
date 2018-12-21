#ifndef  TurbulenceAnalyzerStruct_H
#define  TurbulenceAnalyzerStruct_H

#include "VirtualFluidsDefinitions.h"
#include "Core/DataTypes.h"

struct TurbulenceAnalyzerStruct
{
    uint counter;

    real* U ;
    real* V ;
    real* W ;
    
    real* UU;
    real* VV;
    real* WW;
    
    real* UV;
    real* UW;
    real* VW;
    
    real* T ;
    real* p ;
};

#endif

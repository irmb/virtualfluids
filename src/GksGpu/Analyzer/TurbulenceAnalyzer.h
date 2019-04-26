#ifndef  TurbulenceAnalyzer_H
#define  TurbulenceAnalyzer_H

#include <vector>
#include <string>
#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "FlowStateData/FlowStateData.cuh"

struct DataBase;
struct Parameters;

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
    real* TT;

    real* p ;
};

class VF_PUBLIC TurbulenceAnalyzer
{
private:

    SPtr<DataBase> dataBase;

    uint analyzeStartIter;

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
    real* TT;
    real* p ;

public:

    uint counter;

    std::vector<real> h_U ;
    std::vector<real> h_V ;
    std::vector<real> h_W ;

    std::vector<real> h_UU;
    std::vector<real> h_VV;
    std::vector<real> h_WW;

    std::vector<real> h_UV;
    std::vector<real> h_UW;
    std::vector<real> h_VW;

    std::vector<real> h_T ;
    std::vector<real> h_TT;
    std::vector<real> h_p ;

public:

    ~TurbulenceAnalyzer();

    TurbulenceAnalyzer( SPtr<DataBase> dataBase, uint analyzeStartIter );

    bool run( uint iter, Parameters parameters );

    void writeToFile( std::string filename );

    TurbulenceAnalyzerStruct toStruct();

    void download();
};

#endif

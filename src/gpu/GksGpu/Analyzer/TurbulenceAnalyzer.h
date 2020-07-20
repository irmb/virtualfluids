#ifndef  TurbulenceAnalyzer_H
#define  TurbulenceAnalyzer_H

#include <vector>
#include <string>
#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "FlowStateData/FlowStateData.cuh"

namespace GksGpu {

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

class VIRTUALFLUIDS_GPU_EXPORT TurbulenceAnalyzer
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

    bool collect_U ;
    bool collect_V ;
    bool collect_W ;
    
    bool collect_UU;
    bool collect_VV;
    bool collect_WW;
    
    bool collect_UV;
    bool collect_UW;
    bool collect_VW;
    
    bool collect_T ;
    bool collect_TT;
    bool collect_p ;

public:

    ~TurbulenceAnalyzer();

    TurbulenceAnalyzer( SPtr<DataBase> dataBase, uint analyzeStartIter );

    void free();

    void allocate();

    bool run( uint iter, Parameters parameters );

    void writeRestartFile( std::string filename );

    void readRestartFile( std::string filename );

    TurbulenceAnalyzerStruct toStruct();

    void download(bool normalize = true);

    void upload();
};

} // namespace GksGpu

#endif

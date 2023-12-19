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
//! \file PrecursorWriter.h
//! \author Henry Korb, Henrik Asmuth
//! \date 05/12/2022
//! \brief Probe writing planes of data to be used as inflow data in successor simulation using PrecursorBC
//!
//! The probe writes out yz-planes at a specific x position ( \param xPos ) of either velocity or distributions 
//! that can be read by PrecursorBC as inflow data.
//=======================================================================================


#ifndef PRECURSORPROBE_H_
#define PRECURSORPROBE_H_

#include <future>
#include <string>
#include <vector>

#include <basics/PointerDefinitions.h>

#include <logger/Logger.h>

#include "Calculation/Calculation.h"
#include "PreCollisionInteractor.h"
#include "WbWriterVtkXmlImageBinary.h"

class Parameter;
class CudaMemoryManager;
class GridProvider;

enum class OutputVariable {
   //! - Velocities
    Velocities,
    //! - Distributions
    Distributions
};

static constexpr uint PrecP00 = 0;
static constexpr uint PrecPP0 = 1;
static constexpr uint PrecPM0 = 2;
static constexpr uint PrecP0P = 3;
static constexpr uint PrecP0M = 4;
static constexpr uint PrecPPP = 5;
static constexpr uint PrecPMP = 6;
static constexpr uint PrecPPM = 7;
static constexpr uint PrecPMM = 8;

struct PrecursorStruct
{
    uint numberOfPointsInBC, numberOfPointsInData, numberOfTimestepsPerFile, numberOfFilesWritten, numberOfTimestepsBuffered;
    uint *indicesH, *indicesD;
    real *dataH, *dataD;
    real *bufferH, *bufferD;
    uint numberOfQuantities;
    UbTupleInt4 extent;
    UbTupleFloat2 origin;
    UbTupleFloat3 spacing;
    int* indicesOnPlane;
    cudaStream_t stream;
};

class PrecursorWriter : public PreCollisionInteractor
{
public:
    PrecursorWriter(
        const std::string _fileName,
        const std::string _outputPath,
        real _xPos,
        real _yMin, real _yMax,
        real _zMin, real _zMax,
        uint _tStartOut,
        uint _tSave,
        OutputVariable _outputVariable,
        uint _maxTimestepsPerFile=uint(1e4)
    ): 
    fileName(_fileName), 
    outputPath(_outputPath), 
    xPos(_xPos),
    yMin(_yMin),
    yMax(_yMax),
    zMin(_zMin),
    zMax(_zMax),
    tStartOut(_tStartOut), 
    tSave(_tSave),
    outputVariable(_outputVariable),
    maxtimestepsPerFile(_maxTimestepsPerFile),
    PreCollisionInteractor()
    {
        nodedatanames = determineNodeDataNames();
        writeFuture = std::async([](){});
    };

    ~PrecursorWriter();

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override;

    OutputVariable getOutputVariable(){ return this->outputVariable; }

    SPtr<PrecursorStruct> getPrecursorStruct(int level){return precursorStructs[level];}
    static std::string makeFileName(std::string fileName, int level, int id, uint part);

    void setWritePrecision(uint _writePrecision){ this->writePrecision=_writePrecision;}
    
private:
    WbWriterVtkXmlImageBinary* getWriter(){ return WbWriterVtkXmlImageBinary::getInstance(); };
    void write(int level, uint numberOfTimestepsBuffered);

    std::vector<std::string> determineNodeDataNames()
    {
        switch (outputVariable)
        {
        case OutputVariable::Velocities:
            return {"vx", "vy", "vz"};
            break;       
        case OutputVariable::Distributions:
            return {"fP00", "fPP0", "fPM0", "fP0P", "fP0M", "fPPP", "fPMP", "fPPM", "fPMM"};
            break;
        
        default:
            throw std::runtime_error("Invalid OutputVariable for PrecursorWriter");
            break;
        }
    }

private:
    std::vector<SPtr<PrecursorStruct>> precursorStructs;
    std::string fileName, outputPath;
    std::vector<std::string> nodedatanames;
    std::vector<std::string> celldatanames;
    uint tStartOut, tSave, maxtimestepsPerFile;
    real xPos, yMin, yMax, zMin, zMax;
    OutputVariable outputVariable;
    std::future<void> writeFuture;
    uint writePrecision = 8;
};

#endif //PRECURSORPROBE_H_

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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_PreCollisionInteractor PreCollisionInteractor
//! \ingroup gpu_core core
//! \{
//! \author Henry Korb, Henrik Asmuth
//! \date 05/12/2022
//=======================================================================================


#ifndef PRECURSORPROBE_H_
#define PRECURSORPROBE_H_

#include <future>
#include <string>
#include <vector>

#include <basics/PointerDefinitions.h>
#include <basics/DataTypes.h>

#include <logger/Logger.h>

#include "Calculation/Calculation.h"
#include "PreCollisionInteractor.h"
#include "WbWriterVtkXmlImageBinary.h"

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

//! \brief Probe writing planes of data to be used as inflow data in successor simulation using PrecursorBC
//!
//! The probe writes out yz-planes at a specific x position ( \param xPos ) of either velocity or distributions 
//! that can be read by PrecursorBC as inflow data.
//!
class PrecursorWriter : public PreCollisionInteractor
{
public:
    PrecursorWriter(
        const std::string fileName,
        const std::string outputPath,
        real xPos,
        real yMin, real yMax,
        real zMin, real zMax,
        uint tStartOut,
        uint tSave,
        OutputVariable outputVariable,
        uint maxTimestepsPerFile=uint(1e4)
    ): 
    fileName(fileName), 
    outputPath(outputPath), 
    xPos(xPos),
    yMin(yMin),
    yMax(yMax),
    zMin(zMin),
    zMax(zMax),
    tStartOut(tStartOut), 
    tSave(tSave),
    outputVariable(outputVariable),
    maxtimestepsPerFile(maxTimestepsPerFile)
    {
        nodedatanames = determineNodeDataNames();
        writeFuture = std::async([](){});
    };

    ~PrecursorWriter();

    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override;

    OutputVariable getOutputVariable(){ return this->outputVariable; }

    SPtr<PrecursorStruct> getPrecursorStruct(int level){return precursorStructs[level];}
    static std::string makeFileName(std::string fileName, int level, int id, uint numberOfFilesWritten);

    void setWritePrecision(uint writePrecision){ this->writePrecision=writePrecision;}
    
private:
    void init() override;
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

//! \}

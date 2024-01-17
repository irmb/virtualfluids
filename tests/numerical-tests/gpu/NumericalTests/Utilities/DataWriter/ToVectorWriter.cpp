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
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "ToVectorWriter.h"

#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/Output/FileWriter.h"

#include "Utilities/Structs/VectorWriterInformationStruct.h"

#include "Cuda/CudaMemoryManager.h"

ToVectorWriter::ToVectorWriter(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength)
{
    this->startTimeVectorWriter = vectorWriterInfo->startTimeVectorWriter;
    this->writeVTKFiles = vectorWriterInfo->writeVTKFiles;
    this->startTimeVTKWriter = vectorWriterInfo->startTimeVTKDataWriter;

    this->vtkFileWriter = std::shared_ptr<FileWriter> (new FileWriter());

    this->timeStepLength = timeStepLength;
}

void ToVectorWriter::writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager)
{
    if (startTimeVectorWriter == 0) {
        for (int level = para->getCoarse(); level <= para->getFine(); level++)
            cudaManager->cudaCopyPrint(level);
        writeTimestep(para, 0);
    }
        
    if (writeVTKFiles && startTimeVTKWriter == 0)
        vtkFileWriter->writeTimestep(para, 0);
}

void ToVectorWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int t)
{
    if (startTimeVectorWriter <= t)
    {
        for (int level = para->getCoarse(); level <= para->getFine(); level++)                
            writeTimestep(para, t, level);
    }
    if (writeVTKFiles && startTimeVTKWriter < t)
        vtkFileWriter->writeTimestep(para, t);
}

ToVectorWriter::ToVectorWriter()
{
}

//! \}

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
//! \addtogroup NumericalTestPostProcessing
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "MathematicaFile.h"

#include "Utilities/MathematicaFunction/MathematicaFunktion.h"

#include <ctime>
#include <filesystem>
#include <iomanip>
#include <sstream>

std::shared_ptr<MathematicaFile> MathematicaFile::getNewInstance(std::string filePath)
{
    return std::shared_ptr<MathematicaFile>(new MathematicaFile(filePath));
}

void MathematicaFile::addMathematicaFunction(std::shared_ptr<MathematicaFunction> aMathFunc)
{
    if (!fileFinished)
        myFile << aMathFunc->getFunction() << std::endl;
}

void MathematicaFile::finishFile()
{
    if (!fileFinished) {
        fileFinished = true;
        myFile.close();

        std::filesystem::rename(filePathTxtFile, filePathMathematicaFile);
        std::system(filePathMathematicaFile.c_str());
    }
}

MathematicaFile::MathematicaFile(std::string filePath)
{
    fileFinished = false;

    std::ostringstream basicFilePath, textFile, mathematicaFile;
    basicFilePath << filePath << "/NumericalTestPostProcessing_" << calcDateAndTime();
    textFile << basicFilePath.str() << ".txt";
    mathematicaFile << basicFilePath.str() << ".nb";
    filePathTxtFile = textFile.str();
    filePathMathematicaFile = mathematicaFile.str();

    myFile.open(filePathTxtFile);
}

std::string MathematicaFile::calcDateAndTime()
{
    std::ostringstream oss;
    time_t now = time(NULL);
    struct tm nowLocal = *localtime(&now);
    oss << std::setfill('0') << nowLocal.tm_year + 1900 << std::setw(2) << nowLocal.tm_mon + 1 << std::setw(2) << nowLocal.tm_mday << "_" << std::setw(2) << nowLocal.tm_hour << std::setw(2) << nowLocal.tm_min << std::setw(2) << nowLocal.tm_sec;
    return oss.str();
}

MathematicaFile::MathematicaFile()
{
}

//! \}

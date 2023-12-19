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
#include "ConfigFileReaderNT.h"

#include <basics/config/ConfigurationFile.h>
#include "StringUtilities/StringUtil.h"

#include <memory>
#include <fstream>
#include <string>

#define VAL(str) #str
#define TOSTRING(str) VAL(str)

using ConfigFilePtr = std::shared_ptr<vf::basics::ConfigurationFile>;
using ConfigDataPtr = std::shared_ptr<ConfigDataStruct>;


std::ifstream openConfigFile(const std::string aFilePath)
{
    std::ifstream stream;
    stream.open(aFilePath.c_str(), std::ios::in);
    if (stream.fail())
        throw "can not open config file!\n";

    return stream;
}

bool checkConfigFile(ConfigFilePtr input)
{
    std::vector<double> u0TGVux               = StringUtil::toDoubleVector(input->getValue<std::string>("ux_TGV_Ux"));
    std::vector<double> amplitudeTGVux        = StringUtil::toDoubleVector(input->getValue<std::string>("Amplitude_TGV_Ux"));
    std::vector<int> basisTimeStepLengthTGVux = StringUtil::toIntVector(input->getValue<std::string>("BasisTimeStepLength_TGV_Ux"));

    std::vector<double> v0TGVuz               = StringUtil::toDoubleVector(input->getValue<std::string>("uz_TGV_Uz"));
    std::vector<double> amplitudeTGVuz        = StringUtil::toDoubleVector(input->getValue<std::string>("Amplitude_TGV_Uz"));
    std::vector<int> basisTimeStepLengthTGVuz = StringUtil::toIntVector(input->getValue<std::string>("BasisTimeStepLength_TGV_Uz"));

    std::vector<double> v0SW               = StringUtil::toDoubleVector(input->getValue<std::string>("v0_SW"));
    std::vector<double> u0SW               = StringUtil::toDoubleVector(input->getValue<std::string>("u0_SW"));
    std::vector<int> basisTimeStepLengthSW = StringUtil::toIntVector(input->getValue<std::string>("BasisTimeStepLength_SW"));

    if (u0TGVux.size() != amplitudeTGVux.size() || u0TGVux.size() != basisTimeStepLengthTGVux.size()) {
        std::cout << "Length u0_TGV_U0 is unequal to Lenght Amplitude_TGV_U0 or BasisTimeStepLength_TGV_U0!"
                  << std::endl
                  << std::flush;
        return false;
    } else if (v0TGVuz.size() != amplitudeTGVuz.size() || v0TGVuz.size() != basisTimeStepLengthTGVuz.size()) {
        std::cout << "Length v0_TGV_V0 is unequal to Lenght Amplitude_TGV_V0 or BasisTimeStepLength_TGV_V0!"
                  << std::endl
                  << std::flush;
        return false;
    } else if (u0SW.size() != v0SW.size() || u0SW.size() != basisTimeStepLengthSW.size()) {
        std::cout << "Length u0_SW is unequal to Lenght v0_SW!" << std::endl << std::flush;
        return false;
    } else {
        return true;
    }
}

std::shared_ptr<BasicSimulationParameterStruct>
makeBasicSimulationParameter(ConfigFilePtr input)
{
    std::shared_ptr<BasicSimulationParameterStruct> basicSimPara =
        std::shared_ptr<BasicSimulationParameterStruct>(new BasicSimulationParameterStruct);

    basicSimPara->numberOfTimeSteps = StringUtil::toInt(input->getValue<std::string>("NumberOfTimeSteps"));
    basicSimPara->devices           = StringUtil::toUintVector(input->getValue<std::string>("Devices"));
    return basicSimPara;
}

std::vector<std::shared_ptr<TaylorGreenVortexUxParameterStruct>>
makeTaylorGreenVortexUxParameter(const std::string pathNumericalTests, 
                                 ConfigFilePtr input,
                                 std::shared_ptr<BasicSimulationParameterStruct> basicSimParameter)
{
    std::vector<int> basisTimeStepLength = StringUtil::toIntVector(input->getValue<std::string>("BasisTimeStepLength_TGV_Ux"));
    std::vector<double> amplitude        = StringUtil::toDoubleVector(input->getValue<std::string>("Amplitude_TGV_Ux"));
    std::vector<double> u0               = StringUtil::toDoubleVector(input->getValue<std::string>("ux_TGV_Ux"));
    int l0                               = StringUtil::toInt(input->getValue<std::string>("l0_TGV_Ux"));
    basicSimParameter->l0                = l0;

    std::vector<std::shared_ptr<TaylorGreenVortexUxParameterStruct>> parameter;
    for (int i = 0; i < u0.size(); i++) {
        std::shared_ptr<TaylorGreenVortexUxParameterStruct> aParameter =
            std::shared_ptr<TaylorGreenVortexUxParameterStruct>(new TaylorGreenVortexUxParameterStruct);
        aParameter->basicSimulationParameter = basicSimParameter;

        aParameter->ux                  = u0.at(i);
        aParameter->amplitude           = amplitude.at(i);
        aParameter->basicTimeStepLength = basisTimeStepLength.at(i);
        aParameter->l0                  = l0;
        aParameter->rho0                = StringUtil::toDouble(input->getValue<std::string>("Rho0"));
        aParameter->vtkFilePath         = pathNumericalTests + input->getValue<std::string>("FolderForVTKFileWriting");
        aParameter->dataToCalcTests     = StringUtil::toStringVector(input->getValue<std::string>("DataToCalcTests_TGV_Ux"));
        parameter.push_back(aParameter);
    }
    return parameter;
}

std::vector<std::shared_ptr<TaylorGreenVortexUzParameterStruct>>
makeTaylorGreenVortexUzParameter(const std::string pathNumericalTests,
                                 ConfigFilePtr input,
                                 std::shared_ptr<BasicSimulationParameterStruct> basicSimParameter)
{
    std::vector<int> basisTimeStepLength = StringUtil::toIntVector(input->getValue<std::string>("BasisTimeStepLength_TGV_Uz"));
    std::vector<double> amplitude        = StringUtil::toDoubleVector(input->getValue<std::string>("Amplitude_TGV_Uz"));
    std::vector<double> uz               = StringUtil::toDoubleVector(input->getValue<std::string>("uz_TGV_Uz"));
    int l0                               = StringUtil::toInt(input->getValue<std::string>("l0_TGV_Uz"));
    basicSimParameter->l0                = l0;

    std::vector<std::shared_ptr<TaylorGreenVortexUzParameterStruct>> parameter;
    for (int i = 0; i < uz.size(); i++) {
        std::shared_ptr<TaylorGreenVortexUzParameterStruct> aParameter =
            std::shared_ptr<TaylorGreenVortexUzParameterStruct>(new TaylorGreenVortexUzParameterStruct);
        aParameter->basicSimulationParameter = basicSimParameter;
        aParameter->uz                       = uz.at(i);
        aParameter->amplitude                = amplitude.at(i);
        aParameter->basicTimeStepLength      = basisTimeStepLength.at(i);
        aParameter->l0                       = l0;
        aParameter->rho0                     = StringUtil::toDouble(input->getValue<std::string>("Rho0"));
        aParameter->vtkFilePath              = pathNumericalTests + input->getValue<std::string>("FolderForVTKFileWriting");
        aParameter->dataToCalcTests          = StringUtil::toStringVector(input->getValue<std::string>("DataToCalcTests_TGV_Uz"));
        parameter.push_back(aParameter);
    }
    return parameter;
}
std::vector<std::shared_ptr<ShearWaveParameterStruct>>
makeShearWaveParameter(const std::string pathNumericalTests,
                       ConfigFilePtr input,
                       std::shared_ptr<BasicSimulationParameterStruct> basicSimParameter)
{
    std::vector<int> basisTimeStepLength = StringUtil::toIntVector(input->getValue<std::string>("BasisTimeStepLength_SW"));
    std::vector<double> uz               = StringUtil::toDoubleVector(input->getValue<std::string>("v0_SW"));
    std::vector<double> ux               = StringUtil::toDoubleVector(input->getValue<std::string>("u0_SW"));
    int l0                               = StringUtil::toInt(input->getValue<std::string>("l0_SW"));
    basicSimParameter->l0                = l0;

    std::vector<std::shared_ptr<ShearWaveParameterStruct>> parameter;
    for (int i = 0; i < uz.size(); i++) {
        std::shared_ptr<ShearWaveParameterStruct> aParameter =
            std::shared_ptr<ShearWaveParameterStruct>(new ShearWaveParameterStruct);
        aParameter->basicSimulationParameter = basicSimParameter;
        aParameter->uz                       = uz.at(i);
        aParameter->ux                       = ux.at(i);
        aParameter->basicTimeStepLength      = basisTimeStepLength.at(i);
        aParameter->l0                       = l0;
        aParameter->rho0                     = StringUtil::toDouble(input->getValue<std::string>("Rho0"));
        aParameter->vtkFilePath              = pathNumericalTests + input->getValue<std::string>("FolderForVTKFileWriting");
        aParameter->dataToCalcTests          = StringUtil::toStringVector(input->getValue<std::string>("DataToCalcTests_SW"));
        parameter.push_back(aParameter);
    }
    return parameter;
}

std::shared_ptr<NyTestParameterStruct> makeNyTestParameter(ConfigFilePtr input)
{
    std::shared_ptr<BasicTestParameterStruct> basicTestParameter =
        std::shared_ptr<BasicTestParameterStruct>(new BasicTestParameterStruct);
    basicTestParameter->runTest              = StringUtil::toBool(input->getValue<std::string>("NyTest"));
    basicTestParameter->ySliceForCalculation = StringUtil::toInt(input->getValue<std::string>("ySliceForCalculation"));

    std::shared_ptr<NyTestParameterStruct> testParameter =
        std::shared_ptr<NyTestParameterStruct>(new NyTestParameterStruct);
    testParameter->basicTestParameter       = basicTestParameter;
    testParameter->endTimeStepCalculation   = StringUtil::toInt(input->getValue<std::string>("EndTimeStepCalculation_Ny"));
    testParameter->minOrderOfAccuracy       = StringUtil::toDouble(input->getValue<std::string>("MinOrderOfAccuracy_Ny"));
    testParameter->startTimeStepCalculation = StringUtil::toInt(input->getValue<std::string>("StartTimeStepCalculation_Ny"));

    return testParameter;
}

std::shared_ptr<PhiTestParameterStruct> makePhiTestParameter(ConfigFilePtr input)
{
    std::shared_ptr<BasicTestParameterStruct> basicTestParameter =
        std::shared_ptr<BasicTestParameterStruct>(new BasicTestParameterStruct);
    basicTestParameter->runTest              = StringUtil::toBool(input->getValue<std::string>("PhiTest"));
    basicTestParameter->ySliceForCalculation = StringUtil::toInt(input->getValue<std::string>("ySliceForCalculation"));

    std::shared_ptr<PhiTestParameterStruct> testParameter =
        std::shared_ptr<PhiTestParameterStruct>(new PhiTestParameterStruct);
    testParameter->basicTestParameter       = basicTestParameter;
    testParameter->endTimeStepCalculation   = StringUtil::toInt(input->getValue<std::string>("EndTimeStepCalculation_Phi"));
    testParameter->minOrderOfAccuracy       = StringUtil::toDouble(input->getValue<std::string>("MinOrderOfAccuracy_Phi"));
    testParameter->startTimeStepCalculation = StringUtil::toInt(input->getValue<std::string>("StartTimeStepCalculation_Phi"));

    return testParameter;
}

std::shared_ptr<L2NormTestParameterStruct>
makeL2NormTestParameter(ConfigFilePtr input)
{
    std::shared_ptr<BasicTestParameterStruct> basicTestParameter =
        std::shared_ptr<BasicTestParameterStruct>(new BasicTestParameterStruct);
    basicTestParameter->runTest              = StringUtil::toBool(input->getValue<std::string>("L2NormTest"));
    basicTestParameter->ySliceForCalculation = StringUtil::toInt(input->getValue<std::string>("ySliceForCalculation"));

    std::shared_ptr<L2NormTestParameterStruct> testParameter =
        std::shared_ptr<L2NormTestParameterStruct>(new L2NormTestParameterStruct);
    testParameter->basicTestParameter = basicTestParameter;
    testParameter->basicTimeStep      = StringUtil::toInt(input->getValue<std::string>("BasicTimeStep_L2"));
    testParameter->divergentTimeStep  = StringUtil::toInt(input->getValue<std::string>("DivergentTimeStep_L2"));
    testParameter->normalizeData      = StringUtil::toStringVector(input->getValue<std::string>("NormalizeData_L2Norm"));
    testParameter->maxDiff            = StringUtil::toDoubleVector(input->getValue<std::string>("MaxL2NormDiff"));

    return testParameter;
}

std::vector<std::string> readKernelList(ConfigFilePtr input)
{
    if (StringUtil::toBool(input->getValue<std::string>("L2NormBetweenKernelsTest"))) {
        std::vector<std::string> kernelList = StringUtil::toStringVector(input->getValue<std::string>("KernelsToTest"));
        std::string beginKernel             = input->getValue<std::string>("BasicKernel_L2NormBetweenKernels");
        bool basicKernelInKernelList        = false;
        for (int i = 0; i < kernelList.size(); i++) {
            if (kernelList.at(i) == beginKernel)
                basicKernelInKernelList = true;
        }
        if (!basicKernelInKernelList)
            kernelList.push_back(beginKernel);

        std::vector<std::string> kernelNames = kernelList;

        while (kernelNames.at(0) != beginKernel) {
            kernelNames.push_back(kernelNames.at(0));
            std::vector<std::string>::iterator it = kernelNames.begin();
            kernelNames.erase(it);
        }
        std::vector<std::string> kernels;
        for (int i = 0; i < kernelNames.size(); i++)
            kernels.push_back(kernelNames.at(i));
        return kernels;
    } else {
        std::vector<std::string> kernelList = StringUtil::toStringVector(input->getValue<std::string>("KernelsToTest"));
        std::vector<std::string> kernels;
        for (int i = 0; i < kernelList.size(); i++)
            kernels.push_back(kernelList.at(i));

        return kernels;
    }
}

std::shared_ptr<L2NormTestBetweenKernelsParameterStruct>
makeL2NormTestBetweenKernelsParameter(ConfigFilePtr input)
{
    std::shared_ptr<BasicTestParameterStruct> basicTestParameter =
        std::shared_ptr<BasicTestParameterStruct>(new BasicTestParameterStruct);
    basicTestParameter->runTest              = StringUtil::toBool(input->getValue<std::string>("L2NormBetweenKernelsTest"));
    basicTestParameter->ySliceForCalculation = StringUtil::toInt(input->getValue<std::string>("ySliceForCalculation"));

    std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testParameter =
        std::shared_ptr<L2NormTestBetweenKernelsParameterStruct>(new L2NormTestBetweenKernelsParameterStruct);
    testParameter->basicTestParameter = basicTestParameter;
    testParameter->basicKernel        = input->getValue<std::string>("BasicKernel_L2NormBetweenKernels");
    testParameter->kernelsToTest      = readKernelList(input);
    testParameter->timeSteps          = StringUtil::toIntVector(input->getValue<std::string>("Timesteps_L2NormBetweenKernels"));
    testParameter->normalizeData      = StringUtil::toStringVector(input->getValue<std::string>("NormalizeData_L2Norm"));

    bool correct = false;
    for (int i = 0; i < testParameter->normalizeData.size(); i++)
        if (testParameter->normalizeData.at(i) == "Amplitude" || testParameter->normalizeData.at(i) == "BasicData")
            correct = true;

    if (!correct) {
        std::cout << "invalid input in ConfigFile." << std::endl
                  << "possible data for NormalizeWith Parameter in L2-Norm Test Between Kernels Parameter:" << std::endl
                  << "Amplitude, BasicData" << std::endl
                  << std::endl;
        exit(1);
    }

    return testParameter;
}

std::vector<std::shared_ptr<GridInformationStruct>>
makeGridInformation(const std::string pathNumericalTests, ConfigFilePtr input, std::string simName)
{
    int number = 32;
    std::vector<std::string> valueNames;
    std::vector<std::string> gridPaths;
    for (int i = 1; i <= 5; i++) {
        std::string aValueName = simName;
        aValueName += std::to_string(number);
        valueNames.push_back(aValueName);
        std::string aGridpath = "GridPath";
        aGridpath += std::to_string(number);
        gridPaths.push_back(aGridpath);
        number *= 2;
    }

    std::vector<double> lx;
    std::vector<double> lz;
    std::vector<std::string> gridPath;

    double nextNumber = 32.0;

    for (int i = 0; i < valueNames.size(); i++) {
        if (StringUtil::toBool(input->getValue<std::string>(valueNames.at(i)))) {
            lx.push_back(nextNumber);
            lz.push_back(nextNumber * 3.0 / 2.0);
            gridPath.push_back(pathNumericalTests + input->getValue<std::string>(gridPaths.at(i)));
            nextNumber *= 2;
        }
    }

    std::vector<std::shared_ptr<GridInformationStruct>> gridInformation;
    for (int i = 0; i < lx.size(); i++) {
        std::shared_ptr<GridInformationStruct> aGridInformation =
            std::shared_ptr<GridInformationStruct>(new GridInformationStruct);
        aGridInformation->numberOfGridLevels = StringUtil::toInt(input->getValue<std::string>("NumberOfGridLevels"));
        aGridInformation->maxLevel           = aGridInformation->numberOfGridLevels - 1;
        aGridInformation->gridPath           = gridPath.at(i);
        aGridInformation->lx                 = lx.at(i);
        aGridInformation->lz                 = lz.at(i);
        gridInformation.push_back(aGridInformation);
    }
    return gridInformation;
}

unsigned int calcStartStepForToVectorWriter(ConfigFilePtr input)
{
    std::vector<unsigned int> startStepsTests;
    startStepsTests.push_back(StringUtil::toInt(input->getValue<std::string>("BasicTimeStep_L2")));
    startStepsTests.push_back(StringUtil::toInt(input->getValue<std::string>("StartTimeStepCalculation_Ny")));
    startStepsTests.push_back(StringUtil::toInt(input->getValue<std::string>("StartTimeStepCalculation_Phi")));
    std::sort(startStepsTests.begin(), startStepsTests.end());

    return startStepsTests.at(0);
}

std::shared_ptr<VectorWriterInformationStruct>
makeVectorWriterInformationStruct(ConfigFilePtr input)
{
    std::shared_ptr<VectorWriterInformationStruct> vectorWriter =
        std::shared_ptr<VectorWriterInformationStruct>(new VectorWriterInformationStruct);
    vectorWriter->startTimeVectorWriter  = calcStartStepForToVectorWriter(input);
    vectorWriter->startTimeVTKDataWriter = StringUtil::toInt(input->getValue<std::string>("StartStepFileWriter"));
    vectorWriter->writeVTKFiles          = StringUtil::toBool(input->getValue<std::string>("WriteVTKFiles"));

    return vectorWriter;
}

std::shared_ptr<LogFileParameterStruct> makeLogFilePara(ConfigFilePtr input)
{
    std::shared_ptr<LogFileParameterStruct> logFilePara =
        std::shared_ptr<LogFileParameterStruct>(new LogFileParameterStruct);
    logFilePara->devices              = StringUtil::toIntVector(input->getValue<std::string>("Devices"));
    logFilePara->numberOfTimeSteps    = StringUtil::toInt(input->getValue<std::string>("NumberOfTimeSteps"));
    logFilePara->writeAnalyticalToVTK = StringUtil::toBool(input->getValue<std::string>("WriteAnalyResultsToVTK"));

    return logFilePara;
}

int calcNumberOfSimulationGroup(ConfigFilePtr input, std::string simName)
{
    int counter = 0;
    int number  = 32;
    std::vector<std::string> valueNames;
    for (int i = 1; i <= 5; i++) {
        std::string aValueName = simName;
        aValueName += std::to_string(number);
        valueNames.push_back(aValueName);
        number *= 2;
    }
    for (int i = 0; i < valueNames.size(); i++) {
        if (StringUtil::toBool(input->getValue<std::string>(valueNames.at(i))))
            counter++;
    }
    return counter;
}

int calcNumberOfSimulations(ConfigFilePtr input, ConfigDataPtr configData)
{
    int counter = 0;

    int tgvCounterU0 = calcNumberOfSimulationGroup(input, "TaylorGreenVortexUx");
    tgvCounterU0 *= int(StringUtil::toDoubleVector(input->getValue<std::string>("ux_TGV_Ux")).size());
    counter += tgvCounterU0;

    int tgvCounterV0 = calcNumberOfSimulationGroup(input, "TaylorGreenVortexUz");
    ;
    tgvCounterV0 *= int(StringUtil::toDoubleVector(input->getValue<std::string>("uz_TGV_Uz")).size());
    counter += tgvCounterV0;

    int swCounter = calcNumberOfSimulationGroup(input, "ShearWave");
    ;
    swCounter *= int(StringUtil::toDoubleVector(input->getValue<std::string>("u0_SW")).size());
    counter += swCounter;

    counter *= int(StringUtil::toDoubleVector(input->getValue<std::string>("Viscosity")).size());
    counter *= int(configData->kernelsToTest.size());

    return counter;
}

ConfigDataPtr vf::gpu::tests::readConfigFile(const std::string aFilePath, const std::string &pathNumericalTests)
{
    auto configData = std::make_shared<ConfigDataStruct>();
    auto input      = std::make_shared<vf::basics::ConfigurationFile>();
    input->load(aFilePath);

    if (!checkConfigFile(input))
        exit(1);

    std::cout << pathNumericalTests << "\n";

    configData->viscosity            = StringUtil::toDoubleVector(input->getValue<std::string>("Viscosity"));
    configData->kernelsToTest        = readKernelList(input);
    configData->writeAnalyticalToVTK = StringUtil::toBool(input->getValue<std::string>("WriteAnalyResultsToVTK"));
    configData->ySliceForCalculation = StringUtil::toInt(input->getValue<std::string>("ySliceForCalculation"));

    configData->logFilePath         = pathNumericalTests + input->getValue<std::string>("FolderLogFile");
    configData->numberOfSimulations = calcNumberOfSimulations(input, configData);

    auto basicSimPara = makeBasicSimulationParameter(input);

    configData->taylorGreenVortexUxParameter       = makeTaylorGreenVortexUxParameter(pathNumericalTests, input, basicSimPara);
    configData->taylorGreenVortexUxGridInformation = makeGridInformation(pathNumericalTests, input, "TaylorGreenVortexUx");

    configData->taylorGreenVortexUzParameter       = makeTaylorGreenVortexUzParameter(pathNumericalTests, input, basicSimPara);
    configData->taylorGreenVortexUzGridInformation = makeGridInformation(pathNumericalTests, input, "TaylorGreenVortexUz");

    configData->shearWaveParameter       = makeShearWaveParameter(pathNumericalTests, input, basicSimPara);
    configData->shearWaveGridInformation = makeGridInformation(pathNumericalTests, input, "ShearWave");

    configData->phiTestParameter                  = makePhiTestParameter(input);
    configData->nyTestParameter                   = makeNyTestParameter(input);
    configData->l2NormTestParameter               = makeL2NormTestParameter(input);
    configData->l2NormTestBetweenKernelsParameter = makeL2NormTestBetweenKernelsParameter(input);

    configData->vectorWriterInfo = makeVectorWriterInformationStruct(input);

    configData->logFilePara = makeLogFilePara(input);

    return configData;
}
//! \}

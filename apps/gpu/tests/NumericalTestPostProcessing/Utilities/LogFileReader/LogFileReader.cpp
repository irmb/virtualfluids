#include "LogFileReader.h"

#include "Simulation/ShearWave/LogFileData/ShearWaveLogFileDataImp.h"
#include "Simulation/TaylorGreenVortexUx/LogFileData/TaylorGreenVortexUxLogFileDataImp.h"
#include "Simulation/TaylorGreenVortexUz/LogFileData/TaylorGreenVortexUzLogFileDataImp.h"

#include "Tests/L2Norm/LogFileData/L2NormLogFileDataImp.h"
#include "Tests/L2NormBetweenKernels/LogFileData/L2NormBetweenKernelsLogFileDataImp.h"
#include "Tests/NyTest/LogFileData/NyLogFileDataImp.h"
#include "Tests/PhiTest/LogFileData/PhiLogFileDataImp.h"

#include "Utilities/LogFileData/LogFileDataImp.h"

#include <basics/config/ConfigurationFile.h>
#include "StringUtilities/StringUtil.h"

#include "Utilities/AlmostEquals.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string.h>

std::shared_ptr<LogFileReader> LogFileReader::getInstance()
{
    static std::shared_ptr<LogFileReader> uniqueInstance;
    if (!uniqueInstance)
        uniqueInstance = std::shared_ptr<LogFileReader>(new LogFileReader());
    return uniqueInstance;
}

std::shared_ptr<LogFileData> LogFileReader::readLogFileToLogFileData(std::string filePath)
{
    std::shared_ptr<LogFileDataImp> logFileData = LogFileDataImp::getNewInstance();

    auto input = std::make_shared<vf::basics::ConfigurationFile>();
    input->load(filePath);

    logFileData->setFilePath(filePath);
    logFileData->setDate(input->getValue<std::string>("Date"));
    logFileData->setTime(input->getValue<std::string>("Time"));
    logFileData->setGpuDevices(StringUtil::toStringVector(input->getValue<std::string>("GPU_Devices")));

    logFileData->setKernel(input->getValue<std::string>("Kernel"));
    logFileData->setNumberOfTimeSteps(StringUtil::toInt(input->getValue<std::string>("NumberOfTimeSteps")));
    logFileData->setViscosity(StringUtil::toDouble(input->getValue<std::string>("Viscosity")));
    logFileData->setBasisTimeStepLength(StringUtil::toInt(input->getValue<std::string>("BasisTimeStepLength")));

    logFileData->setSimName(input->getValue<std::string>("SimulationName"));




    std::ostringstream simSigniture;
    if (logFileData->getSimName() == "ShearWave") {
        std::vector<double> shearWaveLx = StringUtil::toDoubleVector(input->getValue<std::string>("Lx"));
        logFileData->setBasicGridLengths(shearWaveLx);
        std::vector<int> shearWaveL0;
        std::vector<double> shearWaveUx;
        std::vector<double> shearWaveUz;
        for (int i = 0; i < shearWaveLx.size(); i++) {
            std::ostringstream l0, ux, uz;
            l0 << "l0_" << shearWaveLx.at(i);
            ux << "ux_" << shearWaveLx.at(i);
            uz << "uz_" << shearWaveLx.at(i);
            shearWaveL0.push_back(StringUtil::toInt(input->getValue<std::string>(l0.str())));
            shearWaveUx.push_back(StringUtil::toDouble(input->getValue<std::string>(ux.str())));
            shearWaveUz.push_back(StringUtil::toDouble(input->getValue<std::string>(uz.str())));
        }
        std::shared_ptr<ShearWaveLogFileDataImp> swLogFileData = ShearWaveLogFileDataImp::getNewInstance();
        swLogFileData->setL0(shearWaveL0);
        swLogFileData->setUx(shearWaveUx);
        swLogFileData->setUz(shearWaveUz);
        logFileData->setShearWaveLogFileData(swLogFileData);
        simSigniture << logFileData->getKernel() << "ShearWaveViscosity" << logFileData->getViscosity() << "ux" << shearWaveUx.at(0) << "uz" << shearWaveUz.at(0);
        logFileData->setBasicSimulation(ShearWave);
    }
    if (logFileData->getSimName() == "TaylorGreenVortexUx") {
        std::vector<double> tgvUxLx = StringUtil::toDoubleVector(input->getValue<std::string>("Lx"));
        logFileData->setBasicGridLengths(tgvUxLx);
        std::vector<int> tgvUxL0;
        std::vector<double> tgvUxUx;
        std::vector<double> tgvUxAmp;
        for (int i = 0; i < tgvUxLx.size(); i++) {
            std::ostringstream l0, ux, amplitude;
            l0 << "l0_" << tgvUxLx.at(i);
            ux << "ux_" << tgvUxLx.at(i);
            amplitude << "Amplitude_" << tgvUxLx.at(i);
            tgvUxL0.push_back(StringUtil::toInt(input->getValue<std::string>(l0.str())));
            tgvUxUx.push_back(StringUtil::toDouble(input->getValue<std::string>(ux.str())));
            tgvUxAmp.push_back(StringUtil::toDouble(input->getValue<std::string>(amplitude.str())));
        }
        std::shared_ptr<TaylorGreenVortexUxLogFileDataImp> tgvUxLogFileData = TaylorGreenVortexUxLogFileDataImp::getNewInstance();
        tgvUxLogFileData->setL0(tgvUxL0);
        tgvUxLogFileData->setUx(tgvUxUx);
        tgvUxLogFileData->setAmplitude(tgvUxAmp);
        logFileData->setTaylorGreenVortexUxLogFileData(tgvUxLogFileData);
        simSigniture << logFileData->getKernel() << "TaylorGreenVortexUxViscosity" << logFileData->getViscosity() << "Ux" << tgvUxUx.at(0) << "Amp" << tgvUxAmp.at(0);
        logFileData->setBasicSimulation(TaylorGreenVortexUx);
    }
    if (logFileData->getSimName() == "TaylorGreenVortexUz") {
        std::vector<double> tgvUzLz = StringUtil::toDoubleVector(input->getValue<std::string>("Lx"));
        logFileData->setBasicGridLengths(tgvUzLz);
        std::vector<int> tgvUzL0;
        std::vector<double> tgvUzUz;
        std::vector<double> tgvUzAmp;
        for (int i = 0; i < tgvUzLz.size(); i++) {
            std::ostringstream l0, uz, amplitude;
            l0 << "l0_" << tgvUzLz.at(i);
            uz << "uz_" << tgvUzLz.at(i);
            amplitude << "Amplitude_" << tgvUzLz.at(i);
            tgvUzL0.push_back(StringUtil::toInt(input->getValue<std::string>(l0.str())));
            tgvUzUz.push_back(StringUtil::toDouble(input->getValue<std::string>(uz.str())));
            tgvUzAmp.push_back(StringUtil::toDouble(input->getValue<std::string>(amplitude.str())));
        }
        std::shared_ptr<TaylorGreenVortexUzLogFileDataImp> tgvUzLogFileData = TaylorGreenVortexUzLogFileDataImp::getNewInstance();
        tgvUzLogFileData->setL0(tgvUzL0);
        tgvUzLogFileData->setUz(tgvUzUz);
        tgvUzLogFileData->setAmplitude(tgvUzAmp);
        logFileData->setTaylorGreenVortexUzLogFileData(tgvUzLogFileData);
        simSigniture << logFileData->getKernel() << "TaylorGreenVortexUzViscosity" << logFileData->getViscosity() << "Uz" << tgvUzUz.at(0) << "Amp" << tgvUzAmp.at(0);
        logFileData->setBasicSimulation(TaylorGreenVortexUz);
    }
    char charsToRemove[] = ".-";
    std::string compatibleString = removeCharsFromString(simSigniture.str(), charsToRemove);
    logFileData->setSimulationSigniture(compatibleString);

    std::vector<int> simTime;
    std::vector<double> resultsCheckTime;
    std::vector<double> testTime;
    std::vector<int> analyticalVTKWritingTime;
    for (int i = 0; i < logFileData->getBasicGridLengths().size(); i++) {
        std::ostringstream simTimeOStringStream, resultsCheckTimeOStringStream, testTimeOStringStream, analyticalVTKWritingTimeOStringStream;
        simTimeOStringStream << "SimulationTime_" << logFileData->getBasicGridLengths().at(i);
        resultsCheckTimeOStringStream << "ResultsCheckTime_" << logFileData->getBasicGridLengths().at(i);
        testTimeOStringStream << "TestTime_" << logFileData->getBasicGridLengths().at(i);
        analyticalVTKWritingTimeOStringStream << "AnalyticalVTKFileWritingTime_" << logFileData->getBasicGridLengths().at(i);
        std::string simTimeString = input->getValue<std::string>(simTimeOStringStream.str());
        std::string resultCheckTimeString = input->getValue<std::string>(resultsCheckTimeOStringStream.str());
        std::string testTimeString = input->getValue<std::string>(testTimeOStringStream.str());
        std::string analyticalVTKWritingTimeString = input->getValue<std::string>(analyticalVTKWritingTimeOStringStream.str());
        simTimeString.erase(simTimeString.end() - 3, simTimeString.end());
        resultCheckTimeString.erase(resultCheckTimeString.end() - 3, resultCheckTimeString.end());
        testTimeString.erase(testTimeString.end() - 3, testTimeString.end());
        analyticalVTKWritingTimeString.erase(analyticalVTKWritingTimeString.end() - 3, analyticalVTKWritingTimeString.end());
        simTime.push_back(StringUtil::toInt(simTimeString));
        resultsCheckTime.push_back(StringUtil::toDouble(resultCheckTimeString));
        testTime.push_back(StringUtil::toDouble(testTimeString));
        analyticalVTKWritingTime.push_back(StringUtil::toInt(analyticalVTKWritingTimeString));
    }

    logFileData->setVTKFileWriting(StringUtil::toBool(input->getValue<std::string>("VTKFileWriting")));
    logFileData->setSimTime(simTime);
    logFileData->setResultCheckTime(resultsCheckTime);
    logFileData->setTestTime(testTime);
    logFileData->setAnalyticalVTKWritingTime(analyticalVTKWritingTime);
    
    logFileData->setPhiTestRun(StringUtil::toBool(input->getValue<std::string>("PhiTest")));
    logFileData->setNyTestRun(StringUtil::toBool(input->getValue<std::string>("NyTest")));
    logFileData->setL2NormTestRun(StringUtil::toBool(input->getValue<std::string>("L2NormTest")));
    logFileData->setL2NormTestBetweenKernelRun(StringUtil::toBool(input->getValue<std::string>("L2NormTestBetweenKernel")));

    if (logFileData->getPhiTestRun()) {
        std::vector<std::string> failPhi = StringUtil::toStringVector(input->getValue<std::string>("FailTests_Phi_PhiTest"));
        std::vector<std::string> failOOA = StringUtil::toStringVector(input->getValue<std::string>("FailTests_OOA_PhiTest"));

        std::vector<std::string> dataToCalc = StringUtil::toStringVector(input->getValue<std::string>("DataToCalc_PhiTest"));
        std::vector<std::shared_ptr<PhiLogFileData> > aPhiLogGroup;
        for (int i = 0; i < dataToCalc.size(); i++) {
            std::shared_ptr<PhiLogFileDataImp> phiLog = PhiLogFileDataImp::getNewInstance();
            phiLog->setBasicGridLengths(logFileData->getBasicGridLengths());
            phiLog->setDataToCalc(dataToCalc.at(i));
            phiLog->setStartTimeStepCalculation(StringUtil::toInt(input->getValue<std::string>("StartTimeStepCalculation_PhiTest")));
            phiLog->setEndTimeStepCalculation(StringUtil::toInt(input->getValue<std::string>("EndTimeStepCalculation_PhiTest")));

            std::vector<double> phiDiff;
            std::vector<std::vector<double> > orderOfAccuracy;
            for (int j = 0; j < logFileData->getBasicGridLengths().size(); j++) {
                std::ostringstream phiBasicString, phiString, phiDiffString;
                phiBasicString << logFileData->getBasicGridLengths().at(j) << "_" << dataToCalc.at(i);
                bool failData = false;
                for (int k = 0; k < failPhi.size(); k++) {
                    if (phiBasicString.str() == failPhi.at(k))
                        failData = true;
                }
                if (!failData) {
                    phiDiffString << "PhiDiff_" << logFileData->getBasicGridLengths().at(j) << "_" << dataToCalc.at(i);
                    phiDiff.push_back(StringUtil::toDouble(input->getValue<std::string>(phiDiffString.str())));
                }

                for (int k = j + 1; k < logFileData->getBasicGridLengths().size(); k++) {
                    std::vector<double> aOrderOfAccuracyGroup;
                    std::ostringstream phiDiffOOA, phiDiffBasicOOA;
                    phiDiffBasicOOA << logFileData->getBasicGridLengths().at(j) << "_" << logFileData->getBasicGridLengths().at(k) << "_" << dataToCalc.at(i);
                    bool failData = false;
                    for (int k = 0; k < failOOA.size(); k++) {
                        if (phiDiffBasicOOA.str() == failOOA.at(k))
                            failData = true;
                    }
                    if (!failData) {
                        phiDiffOOA << "OrderOfAccuracy_PhiDiff_" << phiDiffBasicOOA.str();
                        aOrderOfAccuracyGroup.push_back(logFileData->getBasicGridLengths().at(j));
                        aOrderOfAccuracyGroup.push_back(logFileData->getBasicGridLengths().at(k));
                        aOrderOfAccuracyGroup.push_back(StringUtil::toDouble(input->getValue<std::string>(phiDiffOOA.str())));
                    }
                    if (aOrderOfAccuracyGroup.size() > 0)
                        orderOfAccuracy.push_back(aOrderOfAccuracyGroup);

                }


            }
            if (phiDiff.size() > 0) {
                phiLog->setPhiDiff(phiDiff);
            }
            if (orderOfAccuracy.size() > 0)
                phiLog->setOrderOfAccuracy(orderOfAccuracy);
            if (phiDiff.size() > 0 || orderOfAccuracy.size() > 0)
                aPhiLogGroup.push_back(phiLog);
        }
        if (aPhiLogGroup.size() > 0)
            logFileData->setPhiLogFileData(aPhiLogGroup);
        else
            logFileData->setPhiTestRun(false);
    }


    if (logFileData->getNyTestRun()) {
        std::vector<std::string> failNy = StringUtil::toStringVector(input->getValue<std::string>("FailTests_Ny_NyTest"));
        std::vector<std::string> failOOA = StringUtil::toStringVector(input->getValue<std::string>("FailTests_OOA_NyTest"));

        std::vector<std::string> dataToCalc = StringUtil::toStringVector(input->getValue<std::string>("DataToCalc_NyTest"));
        std::vector<std::shared_ptr<NyLogFileData> > aNyLogGroup;
        for (int i = 0; i < dataToCalc.size(); i++) {
            std::shared_ptr<NyLogFileDataImp> nyLog = NyLogFileDataImp::getNewInstance();
            nyLog->setBasicGridLengths(logFileData->getBasicGridLengths());
            nyLog->setDataToCalc(dataToCalc.at(i));
            nyLog->setStartTimeStepCalculation(StringUtil::toInt(input->getValue<std::string>("StartTimeStepCalculation_NyTest")));
            nyLog->setEndTimeStepCalculation(StringUtil::toInt(input->getValue<std::string>("EndTimeStepCalculation_NyTest")));

            std::vector<double> ny, nyDiff;
            std::vector<std::vector<double> > orderOfAccuracy;
            for (int j = 0; j < logFileData->getBasicGridLengths().size(); j++) {
                std::ostringstream nyBasicString, nyString, nyDiffString;
                nyBasicString << logFileData->getBasicGridLengths().at(j) << "_" << dataToCalc.at(i);
                bool failData = false;
                for (int k = 0; k < failNy.size(); k++) {
                    if (nyBasicString.str() == failNy.at(k))
                        failData = true;
                }
                if (!failData) {
                    nyString << "Ny_" << nyBasicString.str();
                    ny.push_back(StringUtil::toDouble(input->getValue<std::string>(nyString.str())));
                    nyDiffString << "NyDiff_" << logFileData->getBasicGridLengths().at(j) << "_" << dataToCalc.at(i);
                    nyDiff.push_back(StringUtil::toDouble(input->getValue<std::string>(nyDiffString.str())));
                }            

                
                for (int k = j + 1; k < logFileData->getBasicGridLengths().size(); k++) {
                    std::vector<double> aOrderOfAccuracyGroup;
                    std::ostringstream nyDiffOOA, nyDiffBasicOOA;
                    nyDiffBasicOOA << logFileData->getBasicGridLengths().at(j) << "_" << logFileData->getBasicGridLengths().at(k) << "_" << dataToCalc.at(i);
                    bool failData = false;
                    for (int k = 0; k < failOOA.size(); k++) {
                        if (nyDiffBasicOOA.str() == failOOA.at(k))
                            failData = true;
                    }
                    if (!failData) {
                        nyDiffOOA << "OrderOfAccuracy_NyDiff_" << nyDiffBasicOOA.str();
                        aOrderOfAccuracyGroup.push_back(logFileData->getBasicGridLengths().at(j));
                        aOrderOfAccuracyGroup.push_back(logFileData->getBasicGridLengths().at(k));
                        aOrderOfAccuracyGroup.push_back(StringUtil::toDouble(input->getValue<std::string>(nyDiffOOA.str())));
                    }
                    if (aOrderOfAccuracyGroup.size() > 0)
                        orderOfAccuracy.push_back(aOrderOfAccuracyGroup);
                        
                }
                

            }
            if (ny.size() > 0) {
                nyLog->setNy(ny);
                nyLog->setNyDiff(nyDiff);
            }
            if (orderOfAccuracy.size() > 0)
                nyLog->setOrderOfAccuracy(orderOfAccuracy);
            if (ny.size() > 0 || orderOfAccuracy.size() > 0)
                aNyLogGroup.push_back(nyLog);
        }
        if (aNyLogGroup.size() > 0)
            logFileData->setNyLogFileData(aNyLogGroup);
        else
            logFileData->setNyTestRun(false);
    }

    if (logFileData->getL2NormTestRun()) {
        std::vector<std::shared_ptr<L2NormLogFileData> > l2NormGroup;
        std::vector<std::string> dataToCalcL2Norm = StringUtil::toStringVector(input->getValue<std::string>("DataToCalc_L2Norm"));
        std::vector<std::string> normData = StringUtil::toStringVector(input->getValue<std::string>("NormalizeData_L2Norm"));
        std::vector<std::string> failL2Norm = StringUtil::toStringVector(input->getValue<std::string>("FailTests_L2Norm"));
        for (int i = 0; i < dataToCalcL2Norm.size(); i++) {
            for (int k = 0; k < normData.size(); k++) {
                std::shared_ptr<L2NormLogFileDataImp> aL2Norm = L2NormLogFileDataImp::getNewInstance();
                aL2Norm->setDataToCalc(dataToCalcL2Norm.at(i));
                aL2Norm->setNormalizeData(normData.at(k));
                aL2Norm->setBasicGridLengths(logFileData->getBasicGridLengths());
                aL2Norm->setBasicTimeStep(StringUtil::toInt(input->getValue<std::string>("BasicTimeStep_L2Norm")));
                aL2Norm->setDivergentTimeStep(StringUtil::toInt(input->getValue<std::string>("DivergentTimeStep_L2Norm")));

                std::vector<double>  l2NormBasicTimeStep;
                std::vector<double>  l2NormDivergentTimeStep;
                std::vector<double>  l2NormDiff;
                for (int j = 0; j < logFileData->getBasicGridLengths().size(); j++) {
                    std::ostringstream basicTimeStep, divergentTimeStep, diff;
                    std::ostringstream basicString;
                    basicString << "L" << logFileData->getBasicGridLengths().at(j) << "_" << dataToCalcL2Norm.at(i) << "_" << normData.at(k);
                    bool fail = false;
                    for (int l = 0; l < failL2Norm.size(); l++)
                        if (basicString.str() == failL2Norm.at(l))
                            fail = true;
                    if (!fail) {
                        basicTimeStep << "L2Norm_BasicTimeStep_" << basicString.str();
                        divergentTimeStep << "L2Norm_DivergentTimeStep_" << basicString.str();
                        diff << "L2Norm_Diff_" << basicString.str();
                        l2NormBasicTimeStep.push_back(StringUtil::toDouble(input->getValue<std::string>(basicTimeStep.str())));
                        l2NormDivergentTimeStep.push_back(StringUtil::toDouble(input->getValue<std::string>(divergentTimeStep.str())));
                        l2NormDiff.push_back(StringUtil::toDouble(input->getValue<std::string>(diff.str())));
                    }
                }
                if (l2NormBasicTimeStep.size() > 0) {
                    if (l2NormBasicTimeStep.size() != logFileData->getBasicGridLengths().size() || l2NormDivergentTimeStep.size() != logFileData->getBasicGridLengths().size() || l2NormDiff.size() != logFileData->getBasicGridLengths().size()) {
                        std::vector<double> lengths;
                        std::vector<std::string> basicStrings;
                        for (int j = 0; j < logFileData->getBasicGridLengths().size(); j++) {
                            std::ostringstream basicString;
                            basicString << "L" << logFileData->getBasicGridLengths().at(j) << "_" << dataToCalcL2Norm.at(i) << "_" << normData.at(k);
                            basicStrings.push_back(basicString.str());
                            lengths.push_back(logFileData->getBasicGridLengths().at(j));
                        }
                        std::vector<double> failLengths;
                        for (int j = 0; j < basicStrings.size(); j++) {
                            bool lengthIsInFail = false;
                            for (int l = 0; l < failL2Norm.size(); l++) {
                                if (basicStrings.at(j) == failL2Norm.at(l))
                                    lengthIsInFail = true;
                            }
                            if (lengthIsInFail)
                                failLengths.push_back(lengths.at(j));
                        }
                        for (int j = 0; j < failLengths.size(); j++) {
                            for (int l = 0; l < lengths.size(); l++) {
                                if (checkEqualDouble(failLengths.at(j), lengths.at(l))) {
                                    std::vector<double>::iterator itBasic = l2NormBasicTimeStep.begin() + l;
                                    l2NormBasicTimeStep.insert(itBasic, 0.0);
                                    std::vector<double>::iterator itDiv = l2NormDivergentTimeStep.begin() + l;
                                    l2NormDivergentTimeStep.insert(itDiv, 0.0);
                                    std::vector<double>::iterator itDiff = l2NormDiff.begin() + l;
                                    l2NormDiff.insert(itDiff, 0.0);
                                }
                            }
                        }
                    }
                    aL2Norm->setL2NormForBasicTimeStep(l2NormBasicTimeStep);
                    aL2Norm->setL2NormForDivergentTimeStep(l2NormDivergentTimeStep);
                    aL2Norm->setL2NormDiff(l2NormDiff);
                    l2NormGroup.push_back(aL2Norm);
                }
            }
        }
        if (l2NormGroup.size() > 0)
            logFileData->setL2NormLogFileData(l2NormGroup);
        else
            logFileData->setL2NormTestRun(false);
    }

    if (logFileData->getL2NormTestBetweenKernelRun()) {
        std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > l2NormBetweenKernelsData;
        std::vector<std::string> dataToCalc = StringUtil::toStringVector(input->getValue<std::string>("DataToCalculate_L2Norm_BK"));
        std::vector<int> timeSteps = StringUtil::toIntVector(input->getValue<std::string>("TimeSteps_L2Norm_BK"));
        std::vector<std::string> normalizeData = StringUtil::toStringVector(input->getValue<std::string>("NormalizeWith_L2Norm_BK"));
        std::vector<std::string> failL2Norm = StringUtil::toStringVector(input->getValue<std::string>("FailTests_L2Norm_BK"));


        for (int i = 0; i < dataToCalc.size(); i++) {
            for (int j = 0; j < timeSteps.size(); j++) {
                for (int k = 0; k < normalizeData.size(); k++) {
                    std::vector<double> l2NormBasicKernel;
                    std::vector<double> l2NormDivergentKernel;
                    std::vector<double> l2NormBetweenKernels;
                    std::shared_ptr<L2NormBetweenKernelsLogFileDataImp> aL2NormLogFileData = L2NormBetweenKernelsLogFileDataImp::getNewInstance();
                    aL2NormLogFileData->setBasicKernel(input->getValue<std::string>("BasicKernel_L2Norm_BK"));
                    aL2NormLogFileData->setDivergentKernel(logFileData->getKernel());
                    aL2NormLogFileData->setDataToCalculate(dataToCalc.at(i));
                    aL2NormLogFileData->setTimeStep(timeSteps.at(j));
                    aL2NormLogFileData->setNormalizeData(normalizeData.at(k));
                    aL2NormLogFileData->setBasicGridLengths(logFileData->getBasicGridLengths());

                    for (int l = 0; l < logFileData->getBasicGridLengths().size(); l++) {
                        std::ostringstream basicKernel, divergentKernel, diff;
                        std::ostringstream basicString;
                        basicString << "L" << logFileData->getBasicGridLengths().at(l) << "_" << dataToCalc.at(i) << "_TimeStep_" << timeSteps.at(j) << "_" << normalizeData.at(k);

                        std::string myString = basicString.str();
                        bool fail = false;
                        for (int m = 0; m < failL2Norm.size(); m++) {
                            if (basicString.str() == failL2Norm.at(m))
                                fail = true;
                        }
                        if (!fail) {
                            basicKernel << "L2Norm_BasicKernel_" << basicString.str();
                            divergentKernel << "L2Norm_DivergentKernel_" << basicString.str();
                            diff << "L2Norm_Between_Kernels_" << basicString.str();
                            l2NormBasicKernel.push_back(StringUtil::toDouble(input->getValue<std::string>(basicKernel.str())));
                            l2NormDivergentKernel.push_back(StringUtil::toDouble(input->getValue<std::string>(divergentKernel.str())));
                            l2NormBetweenKernels.push_back(StringUtil::toDouble(input->getValue<std::string>(diff.str())));
                        }                        
                    }
                    if (l2NormBasicKernel.size() > 0) {
                        if (l2NormBasicKernel.size() != logFileData->getBasicGridLengths().size() || l2NormDivergentKernel.size() != logFileData->getBasicGridLengths().size() || l2NormBetweenKernels.size() != logFileData->getBasicGridLengths().size()) {
                            std::vector<double> lengths;
                            std::vector<std::string> basicStrings;
                            for (int l = 0; l < logFileData->getBasicGridLengths().size(); l++) {
                                std::ostringstream basicString;
                                basicString << "L" << logFileData->getBasicGridLengths().at(l) << "_" << dataToCalc.at(i) << "_TimeStep_" << timeSteps.at(j) << "_" << normalizeData.at(k);
                                basicStrings.push_back(basicString.str());
                                lengths.push_back(logFileData->getBasicGridLengths().at(l));
                            }
                            std::vector<double> failLengths;
                            for (int m = 0; m < basicStrings.size(); m++) {
                                bool lengthIsInFail = false;
                                for (int l = 0; l < failL2Norm.size(); l++) {
                                    if (basicStrings.at(m) == failL2Norm.at(l))
                                        lengthIsInFail = true;
                                }
                                if (lengthIsInFail)
                                    failLengths.push_back(lengths.at(m));
                            }
                            for (int m = 0; m < failLengths.size(); m++) {
                                for (int l = 0; l < lengths.size(); l++) {
                                    if (checkEqualDouble(failLengths.at(m), lengths.at(l))) {
                                        std::vector<double>::iterator itBasic = l2NormBasicKernel.begin() + l;
                                        l2NormBasicKernel.insert(itBasic, 0.0);
                                        std::vector<double>::iterator itDiv = l2NormDivergentKernel.begin() + l;
                                        l2NormDivergentKernel.insert(itDiv, 0.0);
                                        std::vector<double>::iterator itDiff = l2NormBetweenKernels.begin() + l;
                                        l2NormBetweenKernels.insert(itDiff, 0.0);
                                    }

                                }
                            }

                        }
                        aL2NormLogFileData->setL2NormForBasicKernel(l2NormBasicKernel);
                        aL2NormLogFileData->setL2NormForDivergentKernel(l2NormDivergentKernel);
                        aL2NormLogFileData->setL2NormBetweenKernels(l2NormBetweenKernels);
                        l2NormBetweenKernelsData.push_back(aL2NormLogFileData);
                    }
                }
            }
        }
        if (l2NormBetweenKernelsData.size() > 0)
            logFileData->setL2NormBetweenKernelsLogFileData(l2NormBetweenKernelsData);
        else
            logFileData->setL2NormTestBetweenKernelRun(false);
    }

    return logFileData;
}

std::vector<std::shared_ptr<LogFileData> > LogFileReader::readLogFilesInDirectoryToLogFileData(std::string directory)
{
    std::vector<std::shared_ptr<LogFileData> > logFileData;

    std::cout << "seaching for LogFiles in: " << directory << std::endl;
    std::vector<std::string> filePaths = getAllFilesInDir(directory, ".txt");
    std::cout << filePaths.size() << " LogFiles found." << std::endl;
    std::cout << "reading LogFiles.." << std::endl;
    for (int i = 0; i < filePaths.size(); i++) {
        logFileData.push_back(readLogFileToLogFileData(filePaths.at(i)));
    }
    return logFileData;
}

LogFileReader::LogFileReader()
{
}

std::vector<std::string> LogFileReader::getAllFilesInDir(const std::string &dirPath, const std::string &fileExtension)
{
    std::vector<std::string> listOfFiles;
    std::filesystem::path myPath = dirPath;
    if (std::filesystem::exists(myPath) && std::filesystem::is_directory(myPath))
    {
        for (auto& item : std::filesystem::recursive_directory_iterator(myPath))
        {
            if (std::filesystem::is_regular_file(item.path()) && item.path().extension() == fileExtension)
                listOfFiles.push_back(item.path().string());
        }
    }
    return listOfFiles;
}

std::string LogFileReader::removeCharsFromString(std::string str, char * charsToRemove)
{
    for (unsigned int i = 0; i < std::strlen(charsToRemove); ++i)
        str.erase(remove(str.begin(), str.end(), charsToRemove[i]), str.end());
    return str;
}

bool LogFileReader::checkEqualDouble(double one, double two)
{
    const FloatingPoint<double> lhs(one), rhs(two);

    if (lhs.AlmostEquals(rhs))
        return true;
    return false;
}

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
//! \file Parameter.h
//! \ingroup Parameter
//! \author Martin Schoenherr
//=======================================================================================
#include "Parameter.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#pragma clang diagnostic ignored "-Wunused-but-set-parameter"
#endif
#include <curand_kernel.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include "StringUtilities/StringUtil.h"

#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>
#include "Cuda/CudaStreamManager.h"

Parameter::Parameter() : Parameter(1, 0, {}) {}

Parameter::Parameter(const vf::basics::ConfigurationFile* configData) : Parameter(1, 0, configData) {}

Parameter::Parameter(int numberOfProcesses, int myId) : Parameter(numberOfProcesses, myId, {}) {}

Parameter::Parameter(int numberOfProcesses, int myId, std::optional<const vf::basics::ConfigurationFile*> configData)
{
    this->numprocs = numberOfProcesses;
    this->myProcessId = myId;

    this->setQuadricLimiters(0.01, 0.01, 0.01);
    this->setForcing(0.0, 0.0, 0.0);

    if(configData)
        readConfigData(**configData);

    initGridPaths();
    initGridBasePoints();
    initDefaultLBMkernelAllLevels();

    this->cudaStreamManager = std::make_unique<CudaStreamManager>();
}

Parameter::~Parameter() = default;

void Parameter::readConfigData(const vf::basics::ConfigurationFile &configData)
{
    if (configData.contains("NumberOfDevices"))
        this->setMaxDev(configData.getValue<int>("NumberOfDevices"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Devices"))
        this->setDevices(configData.getVector<uint>("Devices"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Path"))
        this->setOutputPath(configData.getValue<std::string>("Path"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Prefix"))
        this->setOutputPrefix(configData.getValue<std::string>("Prefix"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("WriteGrid"))
        this->setPrintFiles(configData.getValue<bool>("WriteGrid"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("GeometryValues"))
        this->setUseGeometryValues(configData.getValue<bool>("GeometryValues"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("calc2ndOrderMoments"))
        this->setCalc2ndOrderMoments(configData.getValue<bool>("calc2ndOrderMoments"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("calc3rdOrderMoments"))
        this->setCalc3rdOrderMoments(configData.getValue<bool>("calc3rdOrderMoments"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("calcHigherOrderMoments"))
        this->setCalcHighOrderMoments(configData.getValue<bool>("calcHigherOrderMoments"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("calcMean"))
        this->setCalcMean(configData.getValue<bool>("calcMean"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("calcCp"))
        this->calcCp = configData.getValue<bool>("calcCp");
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("calcDrafLift"))
        this->calcDragLift = configData.getValue<bool>("calcDrafLift");
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("writeVeloASCIIfiles"))
        this->writeVeloASCII = configData.getValue<bool>("writeVeloASCIIfiles");
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("calcPlaneConc"))
        this->calcPlaneConc = configData.getValue<bool>("calcPlaneConc");
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("UseConcFile"))
        this->setConcFile(configData.getValue<bool>("UseConcFile"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("UseMeasurePoints"))
        this->setUseMeasurePoints(configData.getValue<bool>("UseMeasurePoints"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("UseInitNeq"))
        this->setUseInitNeq(configData.getValue<bool>("UseInitNeq"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("D3Qxx"))
        this->setD3Qxx(configData.getValue<int>("D3Qxx"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("TimeEnd"))
        this->setTimestepEnd(configData.getValue<int>("TimeEnd"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("TimeOut"))
        this->setTimestepOut(configData.getValue<int>("TimeOut"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("TimeStartOut"))
        this->setTimestepStartOut(configData.getValue<int>("TimeStartOut"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("TimeStartCalcMean"))
        this->setTimeCalcMedStart(configData.getValue<int>("TimeStartCalcMean"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("TimeEndCalcMean"))
        this->setTimeCalcMedEnd(configData.getValue<int>("TimeEndCalcMean"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("PressInID"))
        this->setPressInID(configData.getValue<int>("PressInID"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("PressOutID"))
        this->setPressOutID(configData.getValue<int>("PressOutID"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("PressInZ"))
        this->setPressInZ(configData.getValue<int>("PressInZ"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("PressOutZ"))
        this->setPressOutZ(configData.getValue<int>("PressOutZ"));

    //////////////////////////////////////////////////////////////////////////
    // second component
    if (configData.contains("DiffOn"))
        this->setDiffOn(configData.getValue<bool>("DiffOn"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Diffusivity"))
        this->setDiffusivity(configData.getValue<real>("Diffusivity"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Temp"))
        this->setTemperatureInit(configData.getValue<real>("Temp"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("TempBC"))
        this->setTemperatureBC(configData.getValue<real>("TempBC"));

    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Viscosity_LB"))
        this->setViscosityLB(configData.getValue<real>("Viscosity_LB"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Velocity_LB"))
        this->setVelocityLB(configData.getValue<real>("Velocity_LB"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Viscosity_Ratio_World_to_LB"))
        this->setViscosityRatio(configData.getValue<real>("Viscosity_Ratio_World_to_LB"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Velocity_Ratio_World_to_LB"))
        this->setVelocityRatio(configData.getValue<real>("Velocity_Ratio_World_to_LB"));
    // //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Density_Ratio_World_to_LB"))
        this->setDensityRatio(configData.getValue<real>("Density_Ratio_World_to_LB"));

    if (configData.contains("Delta_Press"))
        this->setPressRatio(configData.getValue<real>("Delta_Press"));

    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("SliceRealX"))
        this->setRealX(configData.getValue<real>("SliceRealX"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("SliceRealY"))
        this->setRealY(configData.getValue<real>("SliceRealY"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("FactorPressBC"))
        this->setFactorPressBC(configData.getValue<real>("FactorPressBC"));

    //////////////////////////////////////////////////////////////////////////
    // CUDA streams and optimized communication
    if (this->getNumprocs() > 1) {
        if (configData.contains("useStreams")) {
            if (configData.getValue<bool>("useStreams"))
                this->setUseStreams(true);
        }

        if (configData.contains("useReducedCommunicationInInterpolation")) {
            this->useReducedCommunicationAfterFtoC =
                configData.getValue<bool>("useReducedCommunicationInInterpolation");
        }
    }
    //////////////////////////////////////////////////////////////////////////

    // read Geometry (STL)
    if (configData.contains("ReadGeometry"))
        this->setReadGeo(configData.getValue<bool>("ReadGeometry"));

    if (configData.contains("GeometryC"))
        this->setGeometryFileC(configData.getValue<std::string>("GeometryC"));
    else if (this->getReadGeo())
        throw std::runtime_error("readGeo is true, GeometryC has to be set as well!");

    if (configData.contains("GeometryM"))
        this->setGeometryFileM(configData.getValue<std::string>("GeometryM"));
    else if (this->getReadGeo())
        throw std::runtime_error("readGeo is true, GeometryM has to be set as well!");

    if (configData.contains("GeometryF"))
        this->setGeometryFileF(configData.getValue<std::string>("GeometryF"));
    else if (this->getReadGeo())
        throw std::runtime_error("readGeo is true, GeometryF has to be set as well!");

    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("measureClockCycle"))
        this->setclockCycleForMP(configData.getValue<real>("measureClockCycle"));

    if (configData.contains("measureTimestep"))
        this->settimestepForMP(configData.getValue<uint>("measureTimestep"));

    //////////////////////////////////////////////////////////////////////////

    if (configData.contains("GridPath"))
        this->setGridPath(configData.getValue<std::string>("GridPath"));

    // Forcing
    real forcingX = 0.0;
    real forcingY = 0.0;
    real forcingZ = 0.0;

    if (configData.contains("ForcingX"))
        forcingX = configData.getValue<real>("ForcingX");
    if (configData.contains("ForcingY"))
        forcingY = configData.getValue<real>("ForcingY");
    if (configData.contains("ForcingZ"))
        forcingZ = configData.getValue<real>("ForcingZ");

    this->setForcing(forcingX, forcingY, forcingZ);
    //////////////////////////////////////////////////////////////////////////
    // quadricLimiters
    real quadricLimiterP = (real)0.01;
    real quadricLimiterM = (real)0.01;
    real quadricLimiterD = (real)0.01;

    if (configData.contains("QuadricLimiterP"))
        quadricLimiterP = configData.getValue<real>("QuadricLimiterP");
    if (configData.contains("QuadricLimiterM"))
        quadricLimiterM = configData.getValue<real>("QuadricLimiterM");
    if (configData.contains("QuadricLimiterD"))
        quadricLimiterD = configData.getValue<real>("QuadricLimiterD");

    this->setQuadricLimiters(quadricLimiterP, quadricLimiterM, quadricLimiterD);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Restart
    if (configData.contains("TimeDoCheckPoint"))
        this->setTimeDoCheckPoint(configData.getValue<uint>("TimeDoCheckPoint"));

    if (configData.contains("TimeDoRestart"))
        this->setTimeDoRestart(configData.getValue<uint>("TimeDoRestart"));

    if (configData.contains("DoCheckPoint"))
        this->setDoCheckPoint(configData.getValue<bool>("DoCheckPoint"));

    if (configData.contains("DoRestart"))
        this->setDoRestart(configData.getValue<bool>("DoRestart"));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (configData.contains("NOGL"))
        setMaxLevel(configData.getValue<int>("NOGL"));

    if (configData.contains("GridX"))
        this->setGridX(configData.getVector<int>("GridX"));

    if (configData.contains("GridY"))
        this->setGridY(configData.getVector<int>("GridY"));

    if (configData.contains("GridZ"))
        this->setGridZ(configData.getVector<int>("GridZ"));

    if (configData.contains("DistX"))
        this->setDistX(configData.getVector<int>("DistX"));

    if (configData.contains("DistY"))
        this->setDistY(configData.getVector<int>("DistY"));

    if (configData.contains("DistZ"))
        this->setDistZ(configData.getVector<int>("DistZ"));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Kernel
    if (configData.contains("MainKernelName"))
        this->configureMainKernel(configData.getValue<std::string>("MainKernelName"));

    if (configData.contains("MultiKernelOn"))
        this->setMultiKernelOn(configData.getValue<bool>("MultiKernelOn"));

    if (configData.contains("MultiKernelLevel"))
        this->setMultiKernelLevel(configData.getVector<int>("MultiKernelLevel"));

    if (configData.contains("MultiKernelName"))
        this->setMultiKernel(configData.getVector<std::string>("MultiKernelName"));
}

void Parameter::initGridPaths(){
    std::string gridPath = this->getGridPath();

    // add missing slash to gridPath
    if (gridPath.back() != '/') {
        gridPath += "/";
        this->gridPath = gridPath;
    }

    // for multi-gpu add process id (if not already there)
    if (this->getNumprocs() > 1) {
        gridPath += StringUtil::toString(this->getMyProcessID()) + "/";
        this->gridPath = gridPath;
    }

    //////////////////////////////////////////////////////////////////////////

    this->setgeoVec(gridPath + "geoVec.dat");
    this->setcoordX(gridPath + "coordX.dat");
    this->setcoordY(gridPath + "coordY.dat");
    this->setcoordZ(gridPath + "coordZ.dat");
    this->setneighborX(gridPath + "neighborX.dat");
    this->setneighborY(gridPath + "neighborY.dat");
    this->setneighborZ(gridPath + "neighborZ.dat");
    this->setneighborWSB(gridPath + "neighborWSB.dat");
    this->setscaleCFC(gridPath + "scaleCFC.dat");
    this->setscaleCFF(gridPath + "scaleCFF.dat");
    this->setscaleFCC(gridPath + "scaleFCC.dat");
    this->setscaleFCF(gridPath + "scaleFCF.dat");
    this->setscaleOffsetCF(gridPath + "offsetVecCF.dat");
    this->setscaleOffsetFC(gridPath + "offsetVecFC.dat");
    this->setgeomBoundaryBcQs(gridPath + "geomBoundaryQs.dat");
    this->setgeomBoundaryBcValues(gridPath + "geomBoundaryValues.dat");
    this->setinletBcQs(gridPath + "inletBoundaryQs.dat");
    this->setinletBcValues(gridPath + "inletBoundaryValues.dat");
    this->setoutletBcQs(gridPath + "outletBoundaryQs.dat");
    this->setoutletBcValues(gridPath + "outletBoundaryValues.dat");
    this->settopBcQs(gridPath + "topBoundaryQs.dat");
    this->settopBcValues(gridPath + "topBoundaryValues.dat");
    this->setbottomBcQs(gridPath + "bottomBoundaryQs.dat");
    this->setbottomBcValues(gridPath + "bottomBoundaryValues.dat");
    this->setfrontBcQs(gridPath + "frontBoundaryQs.dat");
    this->setfrontBcValues(gridPath + "frontBoundaryValues.dat");
    this->setbackBcQs(gridPath + "backBoundaryQs.dat");
    this->setbackBcValues(gridPath + "backBoundaryValues.dat");
    this->setnumberNodes(gridPath + "numberNodes.dat");
    this->setLBMvsSI(gridPath + "LBMvsSI.dat");
    this->setmeasurePoints(gridPath + "measurePoints.dat");
    this->setcpTop(gridPath + "cpTop.dat");
    this->setcpBottom(gridPath + "cpBottom.dat");
    this->setcpBottom2(gridPath + "cpBottom2.dat");
    this->setConcentration(gridPath + "conc.dat");

    //////////////////////////////////////////////////////////////////////////
    // Normals - Geometry
    this->setgeomBoundaryNormalX(gridPath + "geomBoundaryNormalX.dat");
    this->setgeomBoundaryNormalY(gridPath + "geomBoundaryNormalY.dat");
    this->setgeomBoundaryNormalZ(gridPath + "geomBoundaryNormalZ.dat");
    // Normals - Inlet
    this->setInflowBoundaryNormalX(gridPath + "inletBoundaryNormalX.dat");
    this->setInflowBoundaryNormalY(gridPath + "inletBoundaryNormalY.dat");
    this->setInflowBoundaryNormalZ(gridPath + "inletBoundaryNormalZ.dat");
    // Normals - Outlet
    this->setOutflowBoundaryNormalX(gridPath + "outletBoundaryNormalX.dat");
    this->setOutflowBoundaryNormalY(gridPath + "outletBoundaryNormalY.dat");
    this->setOutflowBoundaryNormalZ(gridPath + "outletBoundaryNormalZ.dat");
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // for Multi GPU
    if (this->getNumprocs() > 1) {

        // 3D domain decomposition
        std::vector<std::string> sendProcNeighborsX, sendProcNeighborsY, sendProcNeighborsZ;
        std::vector<std::string> recvProcNeighborsX, recvProcNeighborsY, recvProcNeighborsZ;
        for (int i = 0; i < this->getNumprocs(); i++) {
            sendProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xs.dat");
            sendProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Ys.dat");
            sendProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zs.dat");
            recvProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xr.dat");
            recvProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Yr.dat");
            recvProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zr.dat");
        }
        this->setPossNeighborFilesX(sendProcNeighborsX, "send");
        this->setPossNeighborFilesY(sendProcNeighborsY, "send");
        this->setPossNeighborFilesZ(sendProcNeighborsZ, "send");
        this->setPossNeighborFilesX(recvProcNeighborsX, "recv");
        this->setPossNeighborFilesY(recvProcNeighborsY, "recv");
        this->setPossNeighborFilesZ(recvProcNeighborsZ, "recv");

    //////////////////////////////////////////////////////////////////////////
    }
}

void Parameter::initGridBasePoints()
{
    if (this->getGridX().empty())
        this->setGridX(std::vector<int>(this->getMaxLevel() + 1, 32));
    if (this->getGridY().empty())
        this->setGridY(std::vector<int>(this->getMaxLevel() + 1, 32));
    if (this->getGridZ().empty())
        this->setGridZ(std::vector<int>(this->getMaxLevel() + 1, 32));

    if (this->getDistX().empty())
        this->setDistX(std::vector<int>(this->getMaxLevel() + 1, 32));
    if (this->getDistY().empty())
        this->setDistY(std::vector<int>(this->getMaxLevel() + 1, 32));
    if (this->getDistZ().empty())
        this->setDistZ(std::vector<int>(this->getMaxLevel() + 1, 32));
}

void Parameter::initDefaultLBMkernelAllLevels(){
    if (this->getMultiKernelOn() && this->getMultiKernelLevel().empty()) {
        std::vector<int> tmp;
        for (int i = 0; i < this->getMaxLevel() + 1; i++) {
            tmp.push_back(i);
        }
        this->setMultiKernelLevel(tmp);
    }

    if (this->getMultiKernelOn() && this->getMultiKernel().empty()) {
        std::vector<std::string> tmp;
        for (int i = 0; i < this->getMaxLevel() + 1; i++) {
            tmp.push_back("CumulantK17Comp");
        }
        this->setMultiKernel(tmp);
    }
}

void Parameter::initLBMSimulationParameter()
{
    // host
    for (int i = coarse; i <= fine; i++) {
        parH[i]                   = std::make_shared<LBMSimulationParameter>();
        parH[i]->numberofthreads  = 64; // 128;
        parH[i]->gridNX           = getGridX().at(i);
        parH[i]->gridNY           = getGridY().at(i);
        parH[i]->gridNZ           = getGridZ().at(i);
        parH[i]->viscosity        = this->vis * pow((real)2.0, i);
        parH[i]->diffusivity      = this->Diffusivity * pow((real)2.0, i);
        parH[i]->omega            = (real)1.0 / (real(3.0) * parH[i]->viscosity + real(0.5)); // omega :-) not s9 = -1.0f/(3.0f*parH[i]->vis+0.5f);//
        parH[i]->nx               = parH[i]->gridNX + 2 * STARTOFFX;
        parH[i]->ny               = parH[i]->gridNY + 2 * STARTOFFY;
        parH[i]->nz               = parH[i]->gridNZ + 2 * STARTOFFZ;
        parH[i]->size_Mat         = parH[i]->nx * parH[i]->ny * parH[i]->nz;
        parH[i]->sizePlaneXY      = parH[i]->nx * parH[i]->ny;
        parH[i]->sizePlaneYZ      = parH[i]->ny * parH[i]->nz;
        parH[i]->sizePlaneXZ      = parH[i]->nx * parH[i]->nz;
//        parH[i]->mem_size_real    = sizeof(real) * parH[i]->size_Mat;         //DEPRECATED: related to full matrix
//        parH[i]->mem_size_int     = sizeof(unsigned int) * parH[i]->size_Mat; //DEPRECATED: related to full matrix
//        parH[i]->mem_size_bool    = sizeof(bool) * parH[i]->size_Mat;         //DEPRECATED: related to full matrix
//        parH[i]->mem_size_real_yz = sizeof(real) * parH[i]->ny * parH[i]->nz; //DEPRECATED: related to full matrix
        parH[i]->isEvenTimestep        = true;
        parH[i]->startz           = parH[i]->gridNZ * this->myProcessId;
        parH[i]->endz             = parH[i]->gridNZ * this->myProcessId + parH[i]->gridNZ;
        parH[i]->Lx               = ((real)1.0 * parH[i]->gridNX - (real)1.0) / (pow((real)2.0, i));
        parH[i]->Ly               = ((real)1.0 * parH[i]->gridNY - (real)1.0) / (pow((real)2.0, i));
        parH[i]->Lz               = ((real)1.0 * parH[i]->gridNZ - (real)1.0) / (pow((real)2.0, i));
        parH[i]->dx               = (real)1.0 / pow((real)2.0, i);
        parH[i]->XdistKn          = getDistX().at(i);
        parH[i]->YdistKn          = getDistY().at(i);
        parH[i]->ZdistKn          = getDistZ().at(i);
        if (i == coarse) {
            parH[i]->distX  = (real)getDistX().at(i);
            parH[i]->distY  = (real)getDistY().at(i);
            parH[i]->distZ  = (real)getDistZ().at(i);
            parH[i]->mTtoWx = (real)1.0;
            parH[i]->mTtoWy = (real)1.0;
            parH[i]->mTtoWz = (real)1.0;
            parH[i]->cTtoWx = (real)0.0;
            parH[i]->cTtoWy = (real)0.0;
            parH[i]->cTtoWz = (real)0.0;
            ////MGs Trafo///////////////////////////////////////////////////////////////
            // parH[i]->cStartx               = (real)parH[i]->XdistKn;
            // parH[i]->cStarty               = (real)parH[i]->XdistKn;
            // parH[i]->cStartz               = (real)parH[i]->XdistKn;
            ////////////////////////////////////////////////////////////////////////////
        } else {
            // Geller
            parH[i]->distX = ((real)getDistX().at(i) + (real)0.25) * parH[i - 1]->dx;
            parH[i]->distY = ((real)getDistY().at(i) + (real)0.25) * parH[i - 1]->dx;
            parH[i]->distZ = ((real)getDistZ().at(i) + (real)0.25) * parH[i - 1]->dx;
            // parH[i]->distX                 = ((real)getDistX().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distX;
            // parH[i]->distY                 = ((real)getDistY().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distY;
            // parH[i]->distZ                 = ((real)getDistZ().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distZ;
            parH[i]->mTtoWx = (real)pow(0.5f, i);
            parH[i]->mTtoWy = (real)pow(0.5f, i);
            parH[i]->mTtoWz = (real)pow(0.5f, i);
            parH[i]->cTtoWx = (real)(STARTOFFX / 2.f + (parH[i]->gridNX + 1.f) / 4.f); // funzt nur fuer zwei level
            parH[i]->cTtoWy = (real)(STARTOFFY / 2.f + (parH[i]->gridNY + 1.f) / 4.f); // funzt nur fuer zwei level
            parH[i]->cTtoWz = (real)(STARTOFFZ / 2.f + (parH[i]->gridNZ + 1.f) / 4.f); // funzt nur fuer zwei level
            ////MGs Trafo///////////////////////////////////////////////////////////////
            // parH[i]->cStartx               = (real)parH[i]->XdistKn;
            // parH[i]->cStarty               = (real)parH[i]->XdistKn;
            // parH[i]->cStartz               = (real)parH[i]->XdistKn;
            ////////////////////////////////////////////////////////////////////////////
        }
    }

    // device
    for (int i = coarse; i <= fine; i++) {
        parD[i]                   = std::make_shared<LBMSimulationParameter>();
        parD[i]->numberofthreads  = parH[i]->numberofthreads;
        parD[i]->gridNX           = parH[i]->gridNX;
        parD[i]->gridNY           = parH[i]->gridNY;
        parD[i]->gridNZ           = parH[i]->gridNZ;
        parD[i]->viscosity        = parH[i]->viscosity;
        parD[i]->diffusivity      = parH[i]->diffusivity;
        parD[i]->omega            = parH[i]->omega;
        parD[i]->nx               = parH[i]->nx;
        parD[i]->ny               = parH[i]->ny;
        parD[i]->nz               = parH[i]->nz;
        parD[i]->size_Mat         = parH[i]->size_Mat;
        parD[i]->sizePlaneXY      = parH[i]->sizePlaneXY;
        parD[i]->sizePlaneYZ      = parH[i]->sizePlaneYZ;
        parD[i]->sizePlaneXZ      = parH[i]->sizePlaneXZ;
        //parD[i]->mem_size_real    = sizeof(real) * parD[i]->size_Mat;          //DEPRECATED: related to full matrix
        //parD[i]->mem_size_int     = sizeof(unsigned int) * parD[i]->size_Mat;  //DEPRECATED: related to full matrix
        //parD[i]->mem_size_bool    = sizeof(bool) * parD[i]->size_Mat;          //DEPRECATED: related to full matrix
        //parD[i]->mem_size_real_yz = sizeof(real) * parD[i]->ny * parD[i]->nz;  //DEPRECATED: related to full matrix
        parD[i]->isEvenTimestep        = parH[i]->isEvenTimestep;
        parD[i]->startz           = parH[i]->startz;
        parD[i]->endz             = parH[i]->endz;
        parD[i]->Lx               = parH[i]->Lx;
        parD[i]->Ly               = parH[i]->Ly;
        parD[i]->Lz               = parH[i]->Lz;
        parD[i]->dx               = parH[i]->dx;
        parD[i]->XdistKn          = parH[i]->XdistKn;
        parD[i]->YdistKn          = parH[i]->YdistKn;
        parD[i]->ZdistKn          = parH[i]->ZdistKn;
        parD[i]->distX            = parH[i]->distX;
        parD[i]->distY            = parH[i]->distY;
        parD[i]->distZ            = parH[i]->distZ;
    }

    checkParameterValidityCumulantK17();
}

void Parameter::checkParameterValidityCumulantK17() const
{
    if (this->mainKernel != vf::collisionKernel::compressible::K17CompressibleNavierStokes)
        return;

    const real viscosity = this->parH[maxlevel]->viscosity;
    const real viscosityLimit = 1.0 / 42.0;
    if (viscosity > viscosityLimit) {
        VF_LOG_WARNING("The viscosity (in LB units) at level {} is {:1.3g}. It is recommended to keep it smaller than {:1.3g} "
                       "for the CumulantK17 collision kernel.",
                       maxlevel, viscosity, viscosityLimit);
    }

    const real velocity = this->u0;
    const real velocityLimit = 0.1;
    if (velocity > velocityLimit) {
        VF_LOG_WARNING("The velocity (in LB units) is {:1.4g}. It is recommended to keep it smaller than {:1.4g} for the "
                       "CumulantK17 collision kernel.",
                       velocity, velocityLimit);
    }
}

void Parameter::copyMeasurePointsArrayToVector(int lev)
{
    int valuesPerClockCycle = (int)(getclockCycleForMP() / getTimestepForMP());
    for (int i = 0; i < (int)parH[lev]->MP.size(); i++) {
        for (int j = 0; j < valuesPerClockCycle; j++) {
            int index = i * valuesPerClockCycle + j;
            parH[lev]->MP[i].Vx.push_back(parH[lev]->VxMP[index]);
            parH[lev]->MP[i].Vy.push_back(parH[lev]->VyMP[index]);
            parH[lev]->MP[i].Vz.push_back(parH[lev]->VzMP[index]);
            parH[lev]->MP[i].Rho.push_back(parH[lev]->RhoMP[index]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// set-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::setForcing(real forcingX, real forcingY, real forcingZ)
{
    this->hostForcing[0] = forcingX;
    this->hostForcing[1] = forcingY;
    this->hostForcing[2] = forcingZ;
}
void Parameter::setQuadricLimiters(real quadricLimiterP, real quadricLimiterM, real quadricLimiterD)
{
    this->hostQuadricLimiters[0] = quadricLimiterP;
    this->hostQuadricLimiters[1] = quadricLimiterM;
    this->hostQuadricLimiters[2] = quadricLimiterD;
}
void Parameter::setPhi(real inPhi)
{
    Phi = inPhi;
}
void Parameter::setAngularVelocity(real inAngVel)
{
    angularVelocity = inAngVel;
}
void Parameter::setStepEnsight(unsigned int step)
{
    this->stepEnsight = step;
}
void Parameter::setOutputCount(unsigned int outputCount)
{
    this->outputCount = outputCount;
}
void Parameter::setStartTurn(unsigned int inStartTurn)
{
    startTurn = inStartTurn;
}
void Parameter::setDiffOn(bool isDiff)
{
    diffOn = isDiff;
}
void Parameter::setCompOn(bool isComp)
{
    compOn = isComp;
}
void Parameter::setD3Qxx(int d3qxx)
{
    this->D3Qxx = d3qxx;
}
void Parameter::setMaxLevel(int numberOfLevels)
{
    this->maxlevel = numberOfLevels - 1;
    this->fine = this->maxlevel;
    parH.resize(this->maxlevel + 1);
    parD.resize(this->maxlevel + 1);
}
void Parameter::setTimestepEnd(unsigned int tend)
{
    this->tend = tend;
}
void Parameter::setTimestepOut(unsigned int tout)
{
    this->tout = tout;
}
void Parameter::setTimestepStartOut(unsigned int tStartOut)
{
    this->tStartOut = tStartOut;
}
void Parameter::setTimestepOfCoarseLevel(unsigned int timestep)
{
    this->timestep = timestep;
}
void Parameter::setCalcTurbulenceIntensity(bool calcVelocityAndFluctuations)
{
    this->calcVelocityAndFluctuations = calcVelocityAndFluctuations;
}
void Parameter::setCalcMean(bool calcMean)
{
    this->calcMean = calcMean;
}
void Parameter::setCalcDragLift(bool calcDragLift)
{
    this->calcDragLift = calcDragLift;
}
void Parameter::setCalcCp(bool calcCp)
{
    this->calcCp = calcCp;
}
void Parameter::setWriteVeloASCIIfiles(bool writeVeloASCII)
{
    this->writeVeloASCII = writeVeloASCII;
}
void Parameter::setCalcPlaneConc(bool calcPlaneConc)
{
    this->calcPlaneConc = calcPlaneConc;
}
void Parameter::setTimeCalcMedStart(int CalcMedStart)
{
    this->tCalcMedStart = CalcMedStart;
}
void Parameter::setTimeCalcMedEnd(int CalcMedEnd)
{
    this->tCalcMedEnd = CalcMedEnd;
}
void Parameter::setOutputPath(std::string oPath)
{
    // add missing slash to outputPath
    if (oPath.back() != '/')
        oPath += "/";

    this->oPath = oPath;
    this->setPathAndFilename(this->getOutputPath() + this->getOutputPrefix());
}
void Parameter::setOutputPrefix(std::string oPrefix)
{
    this->oPrefix = oPrefix;
    this->setPathAndFilename(this->getOutputPath() + this->getOutputPrefix());
}
void Parameter::setPathAndFilename(std::string fname)
{
    this->fname = fname;
}
void Parameter::setGridPath(std::string gridPath)
{
    this->gridPath = gridPath;
    this->initGridPaths();
}
void Parameter::setPrintFiles(bool printfiles)
{
    this->printFiles = printfiles;
}
void Parameter::setReadGeo(bool readGeo)
{
    this->readGeo = readGeo;
}
void Parameter::setDiffusivity(real Diffusivity)
{
    this->Diffusivity = Diffusivity;
}
void Parameter::setTemperatureInit(real Temp)
{
    this->Temp = Temp;
}
void Parameter::setTemperatureBC(real TempBC)
{
    this->TempBC = TempBC;
}
void Parameter::setViscosityLB(real Viscosity)
{
    this->vis = Viscosity;
}
void Parameter::setVelocityLB(real Velocity)
{
    this->u0 = Velocity;
}
void Parameter::setViscosityRatio(real ViscosityRatio)
{
    this->vis_ratio = ViscosityRatio;
}
void Parameter::setVelocityRatio(real VelocityRatio)
{
    this->u0_ratio = VelocityRatio;
}
void Parameter::setDensityRatio(real DensityRatio)
{
    this->delta_rho = DensityRatio;
}
void Parameter::setPressRatio(real PressRatio)
{
    this->delta_press = PressRatio;
}
real Parameter::getViscosityRatio()
{
    return this->vis_ratio;
}
real Parameter::getVelocityRatio()
{
    return this->u0_ratio;
}
real Parameter::getDensityRatio()
{
    return this->delta_rho;
}
real Parameter::getPressureRatio()
{
    return this->delta_press;
}
real Parameter::getTimeRatio()
{
    return this->getViscosityRatio() * pow(this->getVelocityRatio(), -2);
}
real Parameter::getLengthRatio()
{
    return this->getViscosityRatio() / this->getVelocityRatio();
}
real Parameter::getForceRatio()
{
    return (this->getDensityRatio()+1.0) * this->getVelocityRatio()/this->getTimeRatio();
}
real Parameter::getScaledViscosityRatio(int level)
{
    return this->getViscosityRatio()/(level+1);
}
real Parameter::getScaledVelocityRatio(int level)
{
    return this->getVelocityRatio();
}
real Parameter::getScaledDensityRatio(int level)
{
    return this->getDensityRatio();
}
real Parameter::getScaledPressureRatio(int level)
{
    return this->getPressureRatio();
}
real Parameter::getScaledTimeRatio(int level)
{
    return this->getTimeRatio()/(level+1);
}
real Parameter::getScaledLengthRatio(int level)
{
    return this->getLengthRatio()/(level+1);
}
real Parameter::getScaledForceRatio(int level)
{
    return this->getForceRatio()*(level+1);
}
real Parameter::getScaledStressRatio(int level)
{
    return this->getVelocityRatio()*this->getVelocityRatio();
}
void Parameter::setRealX(real RealX)
{
    this->RealX = RealX;
}
void Parameter::setRealY(real RealY)
{
    this->RealY = RealY;
}
void Parameter::setPressInID(unsigned int PressInID)
{
    this->PressInID = PressInID;
}
void Parameter::setPressOutID(unsigned int PressOutID)
{
    this->PressOutID = PressOutID;
}
void Parameter::setPressInZ(unsigned int PressInZ)
{
    this->PressInZ = PressInZ;
}
void Parameter::setPressOutZ(unsigned int PressOutZ)
{
    this->PressOutZ = PressOutZ;
}
void Parameter::setOutflowPressureCorrectionFactor(real pressBCrhoCorrectionFactor)
{
    this->outflowPressureCorrectionFactor = pressBCrhoCorrectionFactor;
}
void Parameter::setMaxDev(int maxdev)
{
    this->maxdev = maxdev;
}
void Parameter::setMyID(int myid)
{
    this->myProcessId = myid;
}
void Parameter::setNumprocs(int numprocs)
{
    this->numprocs = numprocs;
}
void Parameter::setDevices(std::vector<uint> devices)
{
    this->devices = devices;
}
void Parameter::setGeometryFileC(std::string GeometryFileC)
{
    this->geometryFileC = GeometryFileC;
}
void Parameter::setGeometryFileM(std::string GeometryFileM)
{
    this->geometryFileM = GeometryFileM;
}
void Parameter::setGeometryFileF(std::string GeometryFileF)
{
    this->geometryFileF = GeometryFileF;
}
void Parameter::setRe(real Re)
{
    this->Re = Re;
}
void Parameter::setFactorPressBC(real factorPressBC)
{
    this->factorPressBC = factorPressBC;
}
void Parameter::setIsGeo(bool isGeo)
{
    this->isGeo = isGeo;
}
void Parameter::setIsGeoNormal(bool isGeoNormal)
{
    this->isGeoNormal = isGeoNormal;
}
void Parameter::setIsInflowNormal(bool isInflowNormal)
{
    this->isInflowNormal = isInflowNormal;
}
void Parameter::setIsOutflowNormal(bool isOutflowNormal)
{
    this->isOutflowNormal = isOutflowNormal;
}
void Parameter::setIsCp(bool isCp)
{
    this->isCp = isCp;
}
void Parameter::setConcFile(bool concFile)
{
    this->isConc = concFile;
}
void Parameter::setUseMeasurePoints(bool useMeasurePoints)
{
    this->isMeasurePoints = useMeasurePoints;
}
void Parameter::setUseInitNeq(bool useInitNeq)
{
    this->isInitNeq = useInitNeq;
}
void Parameter::setUseTurbulentViscosity(bool useTurbulentViscosity)
{
    this->isTurbulentViscosity = useTurbulentViscosity;
}
void Parameter::setTurbulenceModel(vf::lbm::TurbulenceModel turbulenceModel)
{
    this->turbulenceModel = turbulenceModel;
}
void Parameter::setSGSConstant(real SGSConstant)
{
    this->SGSConstant = SGSConstant;
}
void Parameter::setHasWallModelMonitor(bool hasWallModelMonitor)
{
    this->hasWallModelMonitor = hasWallModelMonitor;
}

void Parameter::setIsBodyForce(bool isBodyForce)
{
    this->isBodyForce = isBodyForce;
}

void Parameter::setGridX(std::vector<int> GridX)
{
    this->GridX = GridX;
}
void Parameter::setGridY(std::vector<int> GridY)
{
    this->GridY = GridY;
}
void Parameter::setGridZ(std::vector<int> GridZ)
{
    this->GridZ = GridZ;
}
void Parameter::setDistX(std::vector<int> DistX)
{
    this->DistX = DistX;
}
void Parameter::setDistY(std::vector<int> DistY)
{
    this->DistY = DistY;
}
void Parameter::setDistZ(std::vector<int> DistZ)
{
    this->DistZ = DistZ;
}
void Parameter::setScaleLBMtoSI(std::vector<real> scaleLBMtoSI)
{
    this->scaleLBMtoSI = scaleLBMtoSI;
}
void Parameter::setTranslateLBMtoSI(std::vector<real> translateLBMtoSI)
{
    this->translateLBMtoSI = translateLBMtoSI;
}
void Parameter::setMinCoordX(std::vector<real> MinCoordX)
{
    this->minCoordX = MinCoordX;
}
void Parameter::setMinCoordY(std::vector<real> MinCoordY)
{
    this->minCoordY = MinCoordY;
}
void Parameter::setMinCoordZ(std::vector<real> MinCoordZ)
{
    this->minCoordZ = MinCoordZ;
}
void Parameter::setMaxCoordX(std::vector<real> MaxCoordX)
{
    this->maxCoordX = MaxCoordX;
}
void Parameter::setMaxCoordY(std::vector<real> MaxCoordY)
{
    this->maxCoordY = MaxCoordY;
}
void Parameter::setMaxCoordZ(std::vector<real> MaxCoordZ)
{
    this->maxCoordZ = MaxCoordZ;
}
void Parameter::setTempH(TempforBoundaryConditions *TempH)
{
    this->TempH = TempH;
}
void Parameter::setTempD(TempforBoundaryConditions *TempD)
{
    this->TempD = TempD;
}
void Parameter::setTempVelH(TempVelforBoundaryConditions *TempVelH)
{
    this->TempVelH = TempVelH;
}
void Parameter::setTempVelD(TempVelforBoundaryConditions *TempVelD)
{
    this->TempVelD = TempVelD;
}
void Parameter::setTempPressH(TempPressforBoundaryConditions *TempPressH)
{
    this->TempPressH = TempPressH;
}
void Parameter::setTempPressD(TempPressforBoundaryConditions *TempPressD)
{
    this->TempPressD = TempPressD;
}
// void Parameter::setQinflowH(QforBoundaryConditions* QinflowH)
//{
//   this->QinflowH = QinflowH;
//}
// void Parameter::setQinflowD(QforBoundaryConditions* QinflowD)
//{
//   this->QinflowD = QinflowD;
//}
// void Parameter::setQoutflowH(QforBoundaryConditions* QoutflowH)
//{
//   this->QoutflowH = QoutflowH;
//}
// void Parameter::setQoutflowD(QforBoundaryConditions* QoutflowD)
//{
//   this->QoutflowD = QoutflowD;
//}
void Parameter::setkFull(std::string kFull)
{
    this->kFull = kFull;
}
void Parameter::setgeoFull(std::string geoFull)
{
    this->geoFull = geoFull;
}
void Parameter::setgeoVec(std::string geoVec)
{
    this->geoVec = geoVec;
}
void Parameter::setcoordX(std::string coordX)
{
    this->coordX = coordX;
}
void Parameter::setcoordY(std::string coordY)
{
    this->coordY = coordY;
}
void Parameter::setcoordZ(std::string coordZ)
{
    this->coordZ = coordZ;
}
void Parameter::setneighborX(std::string neighborX)
{
    this->neighborX = neighborX;
}
void Parameter::setneighborY(std::string neighborY)
{
    this->neighborY = neighborY;
}
void Parameter::setneighborZ(std::string neighborZ)
{
    this->neighborZ = neighborZ;
}
void Parameter::setneighborWSB(std::string neighborWSB)
{
    this->neighborWSB = neighborWSB;
}
void Parameter::setscaleCFC(std::string scaleCFC)
{
    this->scaleCFC = scaleCFC;
}
void Parameter::setscaleCFF(std::string scaleCFF)
{
    this->scaleCFF = scaleCFF;
}
void Parameter::setscaleFCC(std::string scaleFCC)
{
    this->scaleFCC = scaleFCC;
}
void Parameter::setscaleFCF(std::string scaleFCF)
{
    this->scaleFCF = scaleFCF;
}
void Parameter::setscaleOffsetCF(std::string scaleOffsetCF)
{
    this->scaleOffsetCF = scaleOffsetCF;
}
void Parameter::setscaleOffsetFC(std::string scaleOffsetFC)
{
    this->scaleOffsetFC = scaleOffsetFC;
}
void Parameter::setgeomBoundaryBcQs(std::string geomBoundaryBcQs)
{
    this->geomBoundaryBcQs = geomBoundaryBcQs;
}
void Parameter::setgeomBoundaryBcValues(std::string geomBoundaryBcValues)
{
    this->geomBoundaryBcValues = geomBoundaryBcValues;
}
void Parameter::setnoSlipBcPos(std::string noSlipBcPos)
{
    this->noSlipBcPos = noSlipBcPos;
}
void Parameter::setnoSlipBcQs(std::string noSlipBcQs)
{
    this->noSlipBcQs = noSlipBcQs;
}
void Parameter::setnoSlipBcValue(std::string noSlipBcValue)
{
    this->noSlipBcValue = noSlipBcValue;
}
void Parameter::setnoSlipBcValues(std::string noSlipBcValues)
{
    this->noSlipBcValues = noSlipBcValues;
}
void Parameter::setslipBcPos(std::string slipBcPos)
{
    this->slipBcPos = slipBcPos;
}
void Parameter::setslipBcQs(std::string slipBcQs)
{
    this->slipBcQs = slipBcQs;
}
void Parameter::setslipBcValue(std::string slipBcValue)
{
    this->slipBcValue = slipBcValue;
}
void Parameter::setpressBcPos(std::string pressBcPos)
{
    this->pressBcPos = pressBcPos;
}
void Parameter::setpressBcQs(std::string pressBcQs)
{
    this->pressBcQs = pressBcQs;
}
void Parameter::setpressBcValue(std::string pressBcValue)
{
    this->pressBcValue = pressBcValue;
}
void Parameter::setpressBcValues(std::string pressBcValues)
{
    this->pressBcValues = pressBcValues;
}
void Parameter::setvelBcQs(std::string velBcQs)
{
    this->velBcQs = velBcQs;
}
void Parameter::setvelBcValues(std::string velBcValues)
{
    this->velBcValues = velBcValues;
}
void Parameter::setinletBcQs(std::string inletBcQs)
{
    this->inletBcQs = inletBcQs;
}
void Parameter::setinletBcValues(std::string inletBcValues)
{
    this->inletBcValues = inletBcValues;
}
void Parameter::setoutletBcQs(std::string outletBcQs)
{
    this->outletBcQs = outletBcQs;
}
void Parameter::setoutletBcValues(std::string outletBcValues)
{
    this->outletBcValues = outletBcValues;
}
void Parameter::settopBcQs(std::string topBcQs)
{
    this->topBcQs = topBcQs;
}
void Parameter::settopBcValues(std::string topBcValues)
{
    this->topBcValues = topBcValues;
}
void Parameter::setbottomBcQs(std::string bottomBcQs)
{
    this->bottomBcQs = bottomBcQs;
}
void Parameter::setbottomBcValues(std::string bottomBcValues)
{
    this->bottomBcValues = bottomBcValues;
}
void Parameter::setfrontBcQs(std::string frontBcQs)
{
    this->frontBcQs = frontBcQs;
}
void Parameter::setfrontBcValues(std::string frontBcValues)
{
    this->frontBcValues = frontBcValues;
}
void Parameter::setbackBcQs(std::string backBcQs)
{
    this->backBcQs = backBcQs;
}
void Parameter::setbackBcValues(std::string backBcValues)
{
    this->backBcValues = backBcValues;
}
void Parameter::setwallBcQs(std::string wallBcQs)
{
    this->wallBcQs = wallBcQs;
}
void Parameter::setwallBcValues(std::string wallBcValues)
{
    this->wallBcValues = wallBcValues;
}
void Parameter::setperiodicBcQs(std::string periodicBcQs)
{
    this->periodicBcQs = periodicBcQs;
}
void Parameter::setperiodicBcValues(std::string periodicBcValues)
{
    this->periodicBcValues = periodicBcValues;
}
void Parameter::setmeasurePoints(std::string measurePoints)
{
    this->measurePoints = measurePoints;
}
void Parameter::setnumberNodes(std::string numberNodes)
{
    this->numberNodes = numberNodes;
}
void Parameter::setLBMvsSI(std::string LBMvsSI)
{
    this->LBMvsSI = LBMvsSI;
}
void Parameter::setcpTop(std::string cpTop)
{
    this->cpTop = cpTop;
}
void Parameter::setcpBottom(std::string cpBottom)
{
    this->cpBottom = cpBottom;
}
void Parameter::setcpBottom2(std::string cpBottom2)
{
    this->cpBottom2 = cpBottom2;
}
void Parameter::setConcentration(std::string concFile)
{
    this->concentration = concFile;
}
void Parameter::setclockCycleForMP(real clockCycleForMP)
{
    this->clockCycleForMP = clockCycleForMP;
}
void Parameter::setTimeDoCheckPoint(unsigned int tDoCheckPoint)
{
    this->tDoCheckPoint = tDoCheckPoint;
}
void Parameter::setTimeDoRestart(unsigned int tDoRestart)
{
    this->tDoRestart = tDoRestart;
}
void Parameter::setDoCheckPoint(bool doCheckPoint)
{
    this->doCheckPoint = doCheckPoint;
}
void Parameter::setDoRestart(bool doRestart)
{
    this->doRestart = doRestart;
}
void Parameter::settimestepForMP(unsigned int timestepForMP)
{
    this->timeStepForMP = timestepForMP;
}
void Parameter::setObj(std::string str, bool isObj)
{
    if (str == "geo") {
        this->setIsGeo(isObj);
    } else if (str == "cp") {
        this->setIsCp(isObj);
    } else if (str == "geoNormal") {
        this->setIsGeoNormal(isObj);
    } else if (str == "inflowNormal") {
        this->setIsInflowNormal(isObj);
    } else if (str == "outflowNormal") {
        this->setIsOutflowNormal(isObj);
    }
}
void Parameter::setUseGeometryValues(bool useGeometryValues)
{
    this->GeometryValues = useGeometryValues;
}
void Parameter::setCalc2ndOrderMoments(bool is2ndOrderMoments)
{
    this->is2ndOrderMoments = is2ndOrderMoments;
}
void Parameter::setCalc3rdOrderMoments(bool is3rdOrderMoments)
{
    this->is3rdOrderMoments = is3rdOrderMoments;
}
void Parameter::setCalcHighOrderMoments(bool isHighOrderMoments)
{
    this->isHighOrderMoments = isHighOrderMoments;
}
void Parameter::setMemsizeGPU(double admem, bool reset)
{
    if (reset == true) {
        this->memsizeGPU = 0.;
    } else {
        this->memsizeGPU += admem;
    }
}
// 1D domain decomposition
void Parameter::setPossNeighborFiles(std::vector<std::string> possNeighborFiles, std::string sor)
{
    if (sor == "send") {
        this->possNeighborFilesSend = possNeighborFiles;
    } else if (sor == "recv") {
        this->possNeighborFilesRecv = possNeighborFiles;
    }
}
void Parameter::setNumberOfProcessNeighbors(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
    if (sor == "send") {
        parH[level]->sendProcessNeighbor.resize(numberOfProcessNeighbors);
        parD[level]->sendProcessNeighbor.resize(numberOfProcessNeighbors);
    } else if (sor == "recv") {
        parH[level]->recvProcessNeighbor.resize(numberOfProcessNeighbors);
        parD[level]->recvProcessNeighbor.resize(numberOfProcessNeighbors);
    }
}
void Parameter::setIsNeighbor(bool isNeigbor)
{
    this->isNeigbor = isNeigbor;
}
// 3D domain decomposition
void Parameter::setPossNeighborFilesX(std::vector<std::string> possNeighborFiles, std::string sor)
{
    if (sor == "send") {
        this->possNeighborFilesSendX = possNeighborFiles;
    } else if (sor == "recv") {
        this->possNeighborFilesRecvX = possNeighborFiles;
    }
}
void Parameter::setPossNeighborFilesY(std::vector<std::string> possNeighborFiles, std::string sor)
{
    if (sor == "send") {
        this->possNeighborFilesSendY = possNeighborFiles;
    } else if (sor == "recv") {
        this->possNeighborFilesRecvY = possNeighborFiles;
    }
}
void Parameter::setPossNeighborFilesZ(std::vector<std::string> possNeighborFiles, std::string sor)
{
    if (sor == "send") {
        this->possNeighborFilesSendZ = possNeighborFiles;
    } else if (sor == "recv") {
        this->possNeighborFilesRecvZ = possNeighborFiles;
    }
}
void Parameter::setNumberOfProcessNeighborsX(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
    if (sor == "send") {
        parH[level]->sendProcessNeighborX.resize(numberOfProcessNeighbors);
        parD[level]->sendProcessNeighborX.resize(numberOfProcessNeighbors);
        //////////////////////////////////////////////////////////////////////////
        if (getDiffOn() == true) {
            parH[level]->sendProcessNeighborADX.resize(numberOfProcessNeighbors);
            parD[level]->sendProcessNeighborADX.resize(numberOfProcessNeighbors);
        }
        //////////////////////////////////////////////////////////////////////////
    } else if (sor == "recv") {
        parH[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
        parD[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
        //////////////////////////////////////////////////////////////////////////
        if (getDiffOn() == true) {
            parH[level]->recvProcessNeighborADX.resize(numberOfProcessNeighbors);
            parD[level]->recvProcessNeighborADX.resize(numberOfProcessNeighbors);
        }
        //////////////////////////////////////////////////////////////////////////
    }
}
void Parameter::setNumberOfProcessNeighborsY(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
    if (sor == "send") {
        parH[level]->sendProcessNeighborY.resize(numberOfProcessNeighbors);
        parD[level]->sendProcessNeighborY.resize(numberOfProcessNeighbors);
        //////////////////////////////////////////////////////////////////////////
        if (getDiffOn() == true) {
            parH[level]->sendProcessNeighborADY.resize(numberOfProcessNeighbors);
            parD[level]->sendProcessNeighborADY.resize(numberOfProcessNeighbors);
        }
        //////////////////////////////////////////////////////////////////////////
    } else if (sor == "recv") {
        parH[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
        parD[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
        //////////////////////////////////////////////////////////////////////////
        if (getDiffOn() == true) {
            parH[level]->recvProcessNeighborADY.resize(numberOfProcessNeighbors);
            parD[level]->recvProcessNeighborADY.resize(numberOfProcessNeighbors);
        }
        //////////////////////////////////////////////////////////////////////////
    }
}
void Parameter::setNumberOfProcessNeighborsZ(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
    if (sor == "send") {
        parH[level]->sendProcessNeighborZ.resize(numberOfProcessNeighbors);
        parD[level]->sendProcessNeighborZ.resize(numberOfProcessNeighbors);
        //////////////////////////////////////////////////////////////////////////
        if (getDiffOn() == true) {
            parH[level]->sendProcessNeighborADZ.resize(numberOfProcessNeighbors);
            parD[level]->sendProcessNeighborADZ.resize(numberOfProcessNeighbors);
        }
        //////////////////////////////////////////////////////////////////////////
    } else if (sor == "recv") {
        parH[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
        parD[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
        //////////////////////////////////////////////////////////////////////////
        if (getDiffOn() == true) {
            parH[level]->recvProcessNeighborADZ.resize(numberOfProcessNeighbors);
            parD[level]->recvProcessNeighborADZ.resize(numberOfProcessNeighbors);
        }
        //////////////////////////////////////////////////////////////////////////
    }
}
void Parameter::setIsNeighborX(bool isNeigbor)
{
    this->isNeigborX = isNeigbor;
}
void Parameter::setIsNeighborY(bool isNeigbor)
{
    this->isNeigborY = isNeigbor;
}
void Parameter::setIsNeighborZ(bool isNeigbor)
{
    this->isNeigborZ = isNeigbor;
}
void Parameter::setSendProcessNeighborsAfterFtoCX(int numberOfNodes, int level, int arrayIndex)
{
    this->getParH(level)->sendProcessNeighborsAfterFtoCX[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCX[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParH(level)->sendProcessNeighborsAfterFtoCX[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCX[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParH(level)->sendProcessNeighborsAfterFtoCX[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCX[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
}
void Parameter::setSendProcessNeighborsAfterFtoCY(int numberOfNodes, int level, int arrayIndex)
{
    this->getParH(level)->sendProcessNeighborsAfterFtoCY[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCY[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParH(level)->sendProcessNeighborsAfterFtoCY[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCY[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParH(level)->sendProcessNeighborsAfterFtoCY[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCY[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
}
void Parameter::setSendProcessNeighborsAfterFtoCZ(int numberOfNodes, int level, int arrayIndex)
{
    this->getParH(level)->sendProcessNeighborsAfterFtoCZ[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCZ[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParH(level)->sendProcessNeighborsAfterFtoCZ[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCZ[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParH(level)->sendProcessNeighborsAfterFtoCZ[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
    this->getParD(level)->sendProcessNeighborsAfterFtoCZ[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
}
void Parameter::setRecvProcessNeighborsAfterFtoCX(int numberOfNodes, int level, int arrayIndex)
{
    this->getParH(level)->recvProcessNeighborsAfterFtoCX[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCX[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParH(level)->recvProcessNeighborsAfterFtoCX[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCX[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParH(level)->recvProcessNeighborsAfterFtoCX[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCX[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
}
void Parameter::setRecvProcessNeighborsAfterFtoCY(int numberOfNodes, int level, int arrayIndex)
{
    this->getParH(level)->recvProcessNeighborsAfterFtoCY[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCY[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParH(level)->recvProcessNeighborsAfterFtoCY[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCY[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParH(level)->recvProcessNeighborsAfterFtoCY[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCY[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
}
void Parameter::setRecvProcessNeighborsAfterFtoCZ(int numberOfNodes, int level, int arrayIndex)
{
    this->getParH(level)->recvProcessNeighborsAfterFtoCZ[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCZ[arrayIndex].numberOfNodes = numberOfNodes;
    this->getParH(level)->recvProcessNeighborsAfterFtoCZ[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCZ[arrayIndex].memsizeFs     = sizeof(real) * numberOfNodes;
    this->getParH(level)->recvProcessNeighborsAfterFtoCZ[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
    this->getParD(level)->recvProcessNeighborsAfterFtoCZ[arrayIndex].numberOfFs    = this->D3Qxx * numberOfNodes;
}
void Parameter::setgeomBoundaryNormalX(std::string geomNormalX)
{
    this->geomNormalX = geomNormalX;
}
void Parameter::setgeomBoundaryNormalY(std::string geomNormalY)
{
    this->geomNormalY = geomNormalY;
}
void Parameter::setgeomBoundaryNormalZ(std::string geomNormalZ)
{
    this->geomNormalZ = geomNormalZ;
}
void Parameter::setInflowBoundaryNormalX(std::string inflowNormalX)
{
    this->inflowNormalX = inflowNormalX;
}
void Parameter::setInflowBoundaryNormalY(std::string inflowNormalY)
{
    this->inflowNormalY = inflowNormalY;
}
void Parameter::setInflowBoundaryNormalZ(std::string inflowNormalZ)
{
    this->inflowNormalZ = inflowNormalZ;
}
void Parameter::setOutflowBoundaryNormalX(std::string outflowNormalX)
{
    this->outflowNormalX = outflowNormalX;
}
void Parameter::setOutflowBoundaryNormalY(std::string outflowNormalY)
{
    this->outflowNormalY = outflowNormalY;
}
void Parameter::setOutflowBoundaryNormalZ(std::string outflowNormalZ)
{
    this->outflowNormalZ = outflowNormalZ;
}
void Parameter::configureMainKernel(std::string kernel)
{
    this->mainKernel = kernel;
    if (kernel == vf::collisionKernel::compressible::K17CompressibleNavierStokes)
        this->kernelNeedsFluidNodeIndicesToRun = true;
}
void Parameter::setMultiKernelOn(bool isOn)
{
    this->multiKernelOn = isOn;
}
void Parameter::setMultiKernelLevel(std::vector<int> kernelLevel)
{
    this->multiKernelLevel = kernelLevel;
}
void Parameter::setMultiKernel(std::vector<std::string> kernel)
{
    this->multiKernel = kernel;
}
void Parameter::setADKernel(std::string adKernel)
{
    this->adKernel = adKernel;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// add-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::addActuator(SPtr<PreCollisionInteractor> actuator)
{
    actuators.push_back(actuator);
}
void Parameter::addProbe(SPtr<PreCollisionInteractor> probe)
{
    probes.push_back(probe);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// get-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double *Parameter::getForcesDouble()
{
    return this->hostForcing;
}
real *Parameter::getForcesHost()
{
    return this->forcingH;
}
real *Parameter::getForcesDev()
{
    return this->forcingD;
}
double *Parameter::getQuadricLimitersDouble()
{
    return this->hostQuadricLimiters;
}
real *Parameter::getQuadricLimitersHost()
{
    return this->quadricLimitersH;
}
real *Parameter::getQuadricLimitersDev()
{
    return this->quadricLimitersD;
}
real Parameter::getPhi()
{
    return Phi;
}
real Parameter::getAngularVelocity()
{
    return angularVelocity;
}
unsigned int Parameter::getStepEnsight()
{
    return this->stepEnsight;
}
unsigned int Parameter::getOutputCount()
{
    return this->outputCount;
}
unsigned int Parameter::getStartTurn()
{
    return startTurn;
}
std::shared_ptr<LBMSimulationParameter> Parameter::getParD(int level)
{
    return parD[level];
}
std::shared_ptr<LBMSimulationParameter> Parameter::getParH(int level)
{
    return parH[level];
}

LBMSimulationParameter& Parameter::getParDeviceAsReference(int level) const
{
    return *parD[level];
}

LBMSimulationParameter& Parameter::getParHostAsReference(int level) const
{
    return *parH[level];
}

const std::vector<std::shared_ptr<LBMSimulationParameter>> &Parameter::getParHallLevels()
{
    return parH;
}
const std::vector<std::shared_ptr<LBMSimulationParameter>> &Parameter::getParDallLevels()
{
    return parD;
}

unsigned int Parameter::getSizeMat(int level)
{
    return parH[level]->size_Mat;
}
//unsigned int Parameter::getMemSizereal(int level)      //DEPRECATED: related to full matrix
//{
//    return parH[level]->mem_size_real;
//}
//unsigned int Parameter::getMemSizeInt(int level)     //DEPRECATED: related to full matrix
//{
//    return parH[level]->mem_size_int;
//}
//unsigned int Parameter::getMemSizeBool(int level)    //DEPRECATED: related to full matrix
//{
//    return parH[level]->mem_size_bool;
//}
//unsigned int Parameter::getMemSizerealYZ(int level)  //DEPRECATED: related to full matrix
//{
//    return parH[level]->mem_size_real_yz;
//}
int Parameter::getFine() const
{
    return fine;
}
int Parameter::getCoarse() const
{
    return coarse;
}
bool Parameter::getEvenOrOdd(int level)
{
    return parD[level]->isEvenTimestep;
}
bool Parameter::getDiffOn()
{
    return diffOn;
}
bool Parameter::getCompOn()
{
    return compOn;
}
int Parameter::getFactorNZ()
{
    return factor_gridNZ;
}
int Parameter::getD3Qxx()
{
    return this->D3Qxx;
}
int Parameter::getMaxLevel() const
{
    return this->maxlevel;
}
unsigned int Parameter::getTimestepStart()
{
    if (getDoRestart()) {
        return getTimeDoRestart() + 1;
    } else {
        return 1;
    }
}
unsigned int Parameter::getTimestepInit()
{
    if (getDoRestart()) {
        return getTimeDoRestart();
    } else {
        return 0;
    }
}
unsigned int Parameter::getTimestepEnd()
{
    return this->tend;
}
unsigned int Parameter::getTimestepOut()
{
    return this->tout;
}
unsigned int Parameter::getTimestepStartOut()
{
    return this->tStartOut;
}
bool Parameter::getCalcMean()
{
    return this->calcMean;
}
bool Parameter::getCalcDragLift()
{
    return this->calcDragLift;
}
bool Parameter::getCalcCp()
{
    return this->calcCp;
}
bool Parameter::getWriteVeloASCIIfiles()
{
    return this->writeVeloASCII;
}
bool Parameter::getCalcPlaneConc()
{
    return this->calcPlaneConc;
}
int Parameter::getTimeCalcMedStart()
{
    return this->tCalcMedStart;
}
int Parameter::getTimeCalcMedEnd()
{
    return this->tCalcMedEnd;
}
std::string Parameter::getOutputPath()
{
    return this->oPath;
}
std::string Parameter::getOutputPrefix()
{
    return this->oPrefix;
}
std::string Parameter::getFName() const
{
    return this->fname;
}
std::string Parameter::getGridPath()
{
    return this->gridPath;
}
bool Parameter::getPrintFiles()
{
    return this->printFiles;
}
bool Parameter::getReadGeo()
{
    return this->readGeo;
}
bool Parameter::getCalcTurbulenceIntensity()
{
    return this->calcVelocityAndFluctuations;
}
real Parameter::getDiffusivity()
{
    return this->Diffusivity;
}
real Parameter::getTemperatureInit()
{
    return this->Temp;
}
real Parameter::getTemperatureBC()
{
    return this->TempBC;
}
real Parameter::getViscosity()
{
    return this->vis;
}
real Parameter::getVelocity()
{
    return this->u0;
}
real Parameter::getRealX()
{
    return this->RealX;
}
real Parameter::getRealY()
{
    return this->RealY;
}
unsigned int Parameter::getPressInID()
{
    return this->PressInID;
}
unsigned int Parameter::getPressOutID()
{
    return this->PressOutID;
}
unsigned int Parameter::getPressInZ()
{
    return this->PressInZ;
}
unsigned int Parameter::getPressOutZ()
{
    return this->PressOutZ;
}
real Parameter::getOutflowPressureCorrectionFactor()
{
    return this->outflowPressureCorrectionFactor;
}
int Parameter::getMaxDev()
{
    return this->maxdev;
}
int Parameter::getMyProcessID() const
{
    return this->myProcessId;
}
int Parameter::getNumprocs() const
{
    return this->numprocs;
}
std::vector<uint> Parameter::getDevices()
{
    return this->devices;
}
std::string Parameter::getGeometryFileC()
{
    return this->geometryFileC;
}
std::string Parameter::getGeometryFileM()
{
    return this->geometryFileM;
}
std::string Parameter::getGeometryFileF()
{
    return this->geometryFileF;
}
real Parameter::getRe()
{
    return this->Re;
}
real Parameter::getFactorPressBC()
{
    return this->factorPressBC;
}
std::vector<int> Parameter::getGridX()
{
    return this->GridX;
}
std::vector<int> Parameter::getGridY()
{
    return this->GridY;
}
std::vector<int> Parameter::getGridZ()
{
    return this->GridZ;
}
std::vector<int> Parameter::getDistX()
{
    return this->DistX;
}
std::vector<int> Parameter::getDistY()
{
    return this->DistY;
}
std::vector<int> Parameter::getDistZ()
{
    return this->DistZ;
}
std::vector<real> Parameter::getScaleLBMtoSI()
{
    return this->scaleLBMtoSI;
}
std::vector<real> Parameter::getTranslateLBMtoSI()
{
    return this->translateLBMtoSI;
}
std::vector<real> Parameter::getMinCoordX()
{
    return this->minCoordX;
}
std::vector<real> Parameter::getMinCoordY()
{
    return this->minCoordY;
}
std::vector<real> Parameter::getMinCoordZ()
{
    return this->minCoordZ;
}
std::vector<real> Parameter::getMaxCoordX()
{
    return this->maxCoordX;
}
std::vector<real> Parameter::getMaxCoordY()
{
    return this->maxCoordY;
}
std::vector<real> Parameter::getMaxCoordZ()
{
    return this->maxCoordZ;
}
TempforBoundaryConditions *Parameter::getTempH()
{
    return this->TempH;
}
TempforBoundaryConditions *Parameter::getTempD()
{
    return this->TempD;
}
TempVelforBoundaryConditions *Parameter::getTempVelH()
{
    return this->TempVelH;
}
TempVelforBoundaryConditions *Parameter::getTempVelD()
{
    return this->TempVelD;
}
TempPressforBoundaryConditions *Parameter::getTempPressH()
{
    return this->TempPressH;
}
TempPressforBoundaryConditions *Parameter::getTempPressD()
{
    return this->TempPressD;
}
// QforBoundaryConditions* Parameter::getQinflowH()
//{
//   return this->QinflowH;
//}
// QforBoundaryConditions* Parameter::getQinflowD()
//{
//   return this->QinflowD;
//}
// QforBoundaryConditions* Parameter::getQoutflowH()
//{
//   return this->QoutflowH;
//}
// QforBoundaryConditions* Parameter::getQoutflowD()
//{
//   return this->QoutflowD;
//}
std::string Parameter::getkFull()
{
    return this->kFull;
}
std::string Parameter::getgeoFull()
{
    return this->geoFull;
}
std::string Parameter::getgeoVec()
{
    return this->geoVec;
}
std::string Parameter::getcoordX()
{
    return this->coordX;
}
std::string Parameter::getcoordY()
{
    return this->coordY;
}
std::string Parameter::getcoordZ()
{
    return this->coordZ;
}
std::string Parameter::getneighborX()
{
    return this->neighborX;
}
std::string Parameter::getneighborY()
{
    return this->neighborY;
}
std::string Parameter::getneighborZ()
{
    return this->neighborZ;
}
std::string Parameter::getneighborWSB()
{
    return this->neighborWSB;
}
std::string Parameter::getscaleCFC()
{
    return this->scaleCFC;
}
std::string Parameter::getscaleCFF()
{
    return this->scaleCFF;
}
std::string Parameter::getscaleFCC()
{
    return this->scaleFCC;
}
std::string Parameter::getscaleFCF()
{
    return this->scaleFCF;
}
std::string Parameter::getscaleOffsetCF()
{
    return this->scaleOffsetCF;
}
std::string Parameter::getscaleOffsetFC()
{
    return this->scaleOffsetFC;
}
std::string Parameter::getgeomBoundaryBcQs()
{
    return this->geomBoundaryBcQs;
}
std::string Parameter::getgeomBoundaryBcValues()
{
    return this->geomBoundaryBcValues;
}
std::string Parameter::getnoSlipBcPos()
{
    return this->noSlipBcPos;
}
std::string Parameter::getnoSlipBcQs()
{
    return this->noSlipBcQs;
}
std::string Parameter::getnoSlipBcValue()
{
    return this->noSlipBcValue;
}
std::string Parameter::getnoSlipBcValues()
{
    return this->noSlipBcValues;
}
std::string Parameter::getslipBcPos()
{
    return this->slipBcPos;
}
std::string Parameter::getslipBcQs()
{
    return this->slipBcQs;
}
std::string Parameter::getslipBcValue()
{
    return this->slipBcValue;
}
std::string Parameter::getpressBcPos()
{
    return this->pressBcPos;
}
std::string Parameter::getpressBcQs()
{
    return this->pressBcQs;
}
std::string Parameter::getpressBcValue()
{
    return this->pressBcValue;
}
std::string Parameter::getpressBcValues()
{
    return this->pressBcValues;
}
std::string Parameter::getvelBcQs()
{
    return this->velBcQs;
}
std::string Parameter::getvelBcValues()
{
    return this->velBcValues;
}
std::string Parameter::getinletBcQs()
{
    return this->inletBcQs;
}
std::string Parameter::getinletBcValues()
{
    return this->inletBcValues;
}
std::string Parameter::getoutletBcQs()
{
    return this->outletBcQs;
}
std::string Parameter::getoutletBcValues()
{
    return this->outletBcValues;
}
std::string Parameter::gettopBcQs()
{
    return this->topBcQs;
}
std::string Parameter::gettopBcValues()
{
    return this->topBcValues;
}
std::string Parameter::getbottomBcQs()
{
    return this->bottomBcQs;
}
std::string Parameter::getbottomBcValues()
{
    return this->bottomBcValues;
}
std::string Parameter::getfrontBcQs()
{
    return this->frontBcQs;
}
std::string Parameter::getfrontBcValues()
{
    return this->frontBcValues;
}
std::string Parameter::getbackBcQs()
{
    return this->backBcQs;
}
std::string Parameter::getbackBcValues()
{
    return this->backBcValues;
}
std::string Parameter::getwallBcQs()
{
    return this->wallBcQs;
}
std::string Parameter::getwallBcValues()
{
    return this->wallBcValues;
}
std::string Parameter::getperiodicBcQs()
{
    return this->periodicBcQs;
}
std::string Parameter::getperiodicBcValues()
{
    return this->periodicBcValues;
}
std::string Parameter::getmeasurePoints()
{
    return this->measurePoints;
}
std::string Parameter::getLBMvsSI()
{
    return this->LBMvsSI;
}
std::string Parameter::getnumberNodes()
{
    return this->numberNodes;
}
std::string Parameter::getcpTop()
{
    return this->cpTop;
}
std::string Parameter::getcpBottom()
{
    return this->cpBottom;
}
std::string Parameter::getcpBottom2()
{
    return this->cpBottom2;
}
std::string Parameter::getConcentration()
{
    return this->concentration;
}
real Parameter::getclockCycleForMP()
{
    return this->clockCycleForMP;
}
unsigned int Parameter::getTimeDoCheckPoint()
{
    return this->tDoCheckPoint;
}
unsigned int Parameter::getTimeDoRestart()
{
    return this->tDoRestart;
}

//=======================================================================================
//! \brief Get current (sub)time step of a given level.
//! \param level 
//! \param t current time step (of level 0)
//! \param isPostCollision whether getTimeStep is called post- (before swap) or pre- (after swap) collision
//!
unsigned int Parameter::getTimeStep(int level, unsigned int t, bool isPostCollision)
{
    if(level>this->getMaxLevel()) throw std::runtime_error("Parameter::getTimeStep: level>this->getMaxLevel()!");
    unsigned int tLevel = t;                                                                  
    if(level>0)
    {
        for(int i=1; i<level; i++){ tLevel = 1 + 2*(tLevel-1) + !this->getEvenOrOdd(i); }     
        bool addOne = isPostCollision? !this->getEvenOrOdd(level): this->getEvenOrOdd(level); 
        tLevel = 1 + 2*(tLevel-1) + addOne;
    }
    return tLevel;
}

bool Parameter::getDoCheckPoint()
{
    return this->doCheckPoint;
}
bool Parameter::getDoRestart()
{
    return this->doRestart;
}
bool Parameter::getIsGeo()
{
    return this->isGeo;
}
bool Parameter::getIsGeoNormal()
{
    return this->isGeoNormal;
}
bool Parameter::getIsInflowNormal()
{
    return this->isInflowNormal;
}
bool Parameter::getIsOutflowNormal()
{
    return this->isOutflowNormal;
}
bool Parameter::getIsCp()
{
    return this->isCp;
}
bool Parameter::getConcFile()
{
    return this->isConc;
}
bool Parameter::getUseMeasurePoints()
{
    return this->isMeasurePoints;
}
vf::lbm::TurbulenceModel Parameter::getTurbulenceModel()
{
    return this->turbulenceModel;
}
bool Parameter::getUseTurbulentViscosity()
{
    return this->isTurbulentViscosity;
}
real Parameter::getSGSConstant()
{
    return this->SGSConstant;
}
bool Parameter::getHasWallModelMonitor()
{
    return this->hasWallModelMonitor;
}
std::vector<SPtr<PreCollisionInteractor>> Parameter::getActuators()
{
    return actuators;
}
std::vector<SPtr<PreCollisionInteractor>> Parameter::getProbes()
{
    return probes;
}
bool Parameter::getUseInitNeq()
{
    return this->isInitNeq;
}

bool Parameter::getIsBodyForce()
{
    return this->isBodyForce;
}

bool Parameter::getIsGeometryValues()
{
    return this->GeometryValues;
}
bool Parameter::getCalc2ndOrderMoments()
{
    return this->is2ndOrderMoments;
}
bool Parameter::getCalc3rdOrderMoments()
{
    return this->is3rdOrderMoments;
}
bool Parameter::getCalcHighOrderMoments()
{
    return this->isHighOrderMoments;
}
bool Parameter::overWritingRestart(uint t)
{
    return t == getTimeDoRestart();
}
unsigned int Parameter::getTimestepForMP()
{
    return this->timeStepForMP;
}
unsigned int Parameter::getTimestepOfCoarseLevel()
{
    return this->timestep;
}
double Parameter::getMemsizeGPU()
{
    return this->memsizeGPU;
}
// 1D domain decomposition
std::vector<std::string> Parameter::getPossNeighborFiles(std::string sor)
{
    if (sor == "send") {
        return this->possNeighborFilesSend;
    } else if (sor == "recv") {
        return this->possNeighborFilesRecv;
    }
    throw std::runtime_error("Parameter string invalid.");
}
unsigned int Parameter::getNumberOfProcessNeighbors(int level, std::string sor)
{
    if (sor == "send") {
        return (unsigned int)parH[level]->sendProcessNeighbor.size();
    } else if (sor == "recv") {
        return (unsigned int)parH[level]->recvProcessNeighbor.size();
    }
    throw std::runtime_error("Parameter string invalid.");
}
bool Parameter::getIsNeighbor()
{
    return this->isNeigbor;
}
// 3D domain decomposition
std::vector<std::string> Parameter::getPossNeighborFilesX(std::string sor)
{
    if (sor == "send") {
        return this->possNeighborFilesSendX;
    } else if (sor == "recv") {
        return this->possNeighborFilesRecvX;
    }
    throw std::runtime_error("Parameter string invalid.");
}
std::vector<std::string> Parameter::getPossNeighborFilesY(std::string sor)
{
    if (sor == "send") {
        return this->possNeighborFilesSendY;
    } else if (sor == "recv") {
        return this->possNeighborFilesRecvY;
    }
    throw std::runtime_error("Parameter string invalid.");
}
std::vector<std::string> Parameter::getPossNeighborFilesZ(std::string sor)
{
    if (sor == "send") {
        return this->possNeighborFilesSendZ;
    } else if (sor == "recv") {
        return this->possNeighborFilesRecvZ;
    }
    throw std::runtime_error("Parameter string invalid.");
}
unsigned int Parameter::getNumberOfProcessNeighborsX(int level, std::string sor)
{
    if (sor == "send") {
        return (unsigned int)parH[level]->sendProcessNeighborX.size();
    } else if (sor == "recv") {
        return (unsigned int)parH[level]->recvProcessNeighborX.size();
    }
    throw std::runtime_error("getNumberOfProcessNeighborsX: Parameter string invalid.");
}
unsigned int Parameter::getNumberOfProcessNeighborsY(int level, std::string sor)
{
    if (sor == "send") {
        return (unsigned int)parH[level]->sendProcessNeighborY.size();
    } else if (sor == "recv") {
        return (unsigned int)parH[level]->recvProcessNeighborY.size();
    }
    throw std::runtime_error("getNumberOfProcessNeighborsY: Parameter string invalid.");
}
unsigned int Parameter::getNumberOfProcessNeighborsZ(int level, std::string sor)
{
    if (sor == "send") {
        return (unsigned int)parH[level]->sendProcessNeighborZ.size();
    } else if (sor == "recv") {
        return (unsigned int)parH[level]->recvProcessNeighborZ.size();
    }
    throw std::runtime_error("getNumberOfProcessNeighborsZ: Parameter string invalid.");
}

bool Parameter::getIsNeighborX()
{
    return this->isNeigborX;
}
bool Parameter::getIsNeighborY()
{
    return this->isNeigborY;
}
bool Parameter::getIsNeighborZ()
{
    return this->isNeigborZ;
}
std::string Parameter::getgeomBoundaryNormalX()
{
    return this->geomNormalX;
}
std::string Parameter::getgeomBoundaryNormalY()
{
    return this->geomNormalY;
}
std::string Parameter::getgeomBoundaryNormalZ()
{
    return this->geomNormalZ;
}
std::string Parameter::getInflowBoundaryNormalX()
{
    return this->inflowNormalX;
}
std::string Parameter::getInflowBoundaryNormalY()
{
    return this->inflowNormalY;
}
std::string Parameter::getInflowBoundaryNormalZ()
{
    return this->inflowNormalZ;
}
std::string Parameter::getOutflowBoundaryNormalX()
{
    return this->outflowNormalX;
}
std::string Parameter::getOutflowBoundaryNormalY()
{
    return this->outflowNormalY;
}
std::string Parameter::getOutflowBoundaryNormalZ()
{
    return this->outflowNormalZ;
}

std::string Parameter::getMainKernel()
{
    return mainKernel;
}
bool Parameter::getMultiKernelOn()
{
    return multiKernelOn;
}
std::vector<int> Parameter::getMultiKernelLevel()
{
    return multiKernelLevel;
}
std::vector<std::string> Parameter::getMultiKernel()
{
    return multiKernel;
}
std::string Parameter::getADKernel()
{
    return adKernel;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// initial condition fluid
void Parameter::setInitialCondition(
    std::function<void(real, real, real, real &, real &, real &, real &)> initialCondition)
{
    this->initialCondition = initialCondition;
}

std::function<void(real, real, real, real &, real &, real &, real &)> &Parameter::getInitialCondition()
{
    return this->initialCondition;
}

// initial condition concentration
void Parameter::setInitialConditionAD(std::function<void(real, real, real, real&)> initialConditionAD)
{
    this->initialConditionAD = initialConditionAD;
}

std::function<void(real, real, real, real&)>& Parameter::getInitialConditionAD()
{
    return this->initialConditionAD;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

real Parameter::TrafoXtoWorld(int CoordX, int level)
{
    return (parH[level]->mTtoWx * CoordX + parH[level]->cTtoWx);
}
real Parameter::TrafoYtoWorld(int CoordY, int level)
{
    return (parH[level]->mTtoWy * CoordY + parH[level]->cTtoWy);
}
real Parameter::TrafoZtoWorld(int CoordZ, int level)
{
    return (parH[level]->mTtoWz * CoordZ + parH[level]->cTtoWz);
}
real Parameter::TrafoXtoMGsWorld(int CoordX, int level)
{
    real temp = 0;
    for (int i = 0; i <= level; i++) {
        temp += (parH[i]->XdistKn + 0.25f) * 2.f * parH[i]->dx;
    }
    temp += (real)((CoordX)*parH[level]->dx);
    return temp;
}
real Parameter::TrafoYtoMGsWorld(int CoordY, int level)
{
    real temp = 0;
    for (int i = 0; i <= level; i++) {
        temp += (parH[i]->YdistKn + 0.25f) * 2.f * parH[i]->dx;
    }
    temp += (real)((CoordY)*parH[level]->dx);
    return temp;
}
real Parameter::TrafoZtoMGsWorld(int CoordZ, int level)
{
    real temp = 0;
    for (int i = 0; i <= level; i++) {
        temp += (parH[i]->ZdistKn + 0.25f) * 2.f * parH[i]->dx;
    }
    temp += (real)((CoordZ)*parH[level]->dx);
    return temp;
}

void Parameter::setUseStreams(bool useStreams)
{
    if (useStreams) {
        if (this->getNumprocs() != 1) {
            this->useStreams = useStreams;
            return; 
        } else {
            VF_LOG_INFO( "Can't use streams with only one process!");
        }
    }
    this->useStreams = false;
}

bool Parameter::getUseStreams()
{
    return this->useStreams;
}

std::unique_ptr<CudaStreamManager> &Parameter::getStreamManager()
{
    return this->cudaStreamManager;
}

bool Parameter::getKernelNeedsFluidNodeIndicesToRun()
{
    return this->kernelNeedsFluidNodeIndicesToRun;
}

void Parameter::setKernelNeedsFluidNodeIndicesToRun(bool  kernelNeedsFluidNodeIndicesToRun){
    this->kernelNeedsFluidNodeIndicesToRun = kernelNeedsFluidNodeIndicesToRun;
}

void Parameter::initProcessNeighborsAfterFtoCX(int level)
{
    this->getParH(level)->sendProcessNeighborsAfterFtoCX.resize(this->getParH(level)->sendProcessNeighborX.size());
    this->getParH(level)->recvProcessNeighborsAfterFtoCX.resize(this->getParH(level)->recvProcessNeighborX.size());
    this->getParD(level)->sendProcessNeighborsAfterFtoCX.resize(
        this->getParH(level)->sendProcessNeighborsAfterFtoCX.size());
    this->getParD(level)->recvProcessNeighborsAfterFtoCX.resize(
        this->getParH(level)->recvProcessNeighborsAfterFtoCX.size());
}

void Parameter::initProcessNeighborsAfterFtoCY(int level)
{
    this->getParH(level)->sendProcessNeighborsAfterFtoCY.resize(this->getParH(level)->sendProcessNeighborY.size());
    this->getParH(level)->recvProcessNeighborsAfterFtoCY.resize(this->getParH(level)->recvProcessNeighborY.size());
    this->getParD(level)->sendProcessNeighborsAfterFtoCY.resize(
        this->getParH(level)->sendProcessNeighborsAfterFtoCY.size());
    this->getParD(level)->recvProcessNeighborsAfterFtoCY.resize(
        this->getParH(level)->recvProcessNeighborsAfterFtoCY.size());
}

void Parameter::initProcessNeighborsAfterFtoCZ(int level)
{
    this->getParH(level)->sendProcessNeighborsAfterFtoCZ.resize(this->getParH(level)->sendProcessNeighborZ.size());
    this->getParH(level)->recvProcessNeighborsAfterFtoCZ.resize(this->getParH(level)->recvProcessNeighborZ.size());
    this->getParD(level)->sendProcessNeighborsAfterFtoCZ.resize(
        this->getParH(level)->sendProcessNeighborsAfterFtoCZ.size());
    this->getParD(level)->recvProcessNeighborsAfterFtoCZ.resize(
        this->getParH(level)->recvProcessNeighborsAfterFtoCZ.size());
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

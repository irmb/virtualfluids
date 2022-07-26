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
#ifndef GPU_PARAMETER_H
#define GPU_PARAMETER_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lbm/constants/D3Q27.h"
#include "LBM/LB.h"
#include "PreCollisionInteractor/PreCollisionInteractor.h"

#include "VirtualFluids_GPU_export.h"

struct curandStateXORWOW;
using curandState = struct curandStateXORWOW;
namespace vf:: basics
{
class ConfigurationFile;
}
class CudaStreamManager;

//! \struct LBMSimulationParameter
//! \brief struct holds and manages the LB-parameter of the simulation
//! \brief For this purpose it holds structures and pointer for host and device data, respectively.
struct LBMSimulationParameter {
    //////////////////////////////////////////////////////////////////////////
    //! \brief decides if the simulation timestep is even or odd
    //! \brief this information is important for the esoteric twist
    bool isEvenTimestep;
    //////////////////////////////////////////////////////////////////////////
    //! \brief stores the number of threads per GPU block
    uint numberofthreads;

    // distributions///////////
    // Distributions19 d0;
    Distributions27 d0;  // DEPRECATED: distribution functions for full matrix (not sparse)
    //! \brief store all distribution functions for the D3Q27
    Distributions27 distributions;

    // distributions F3////////
    Distributions6 g6;

    // advection diffusion //////////////////
    //! \brief store all distribution functions for the D3Q7 advection diffusion field
    Distributions7 distributionsAD7;
    //! \brief store all distribution functions for the D3Q27 advection diffusion field
    Distributions27 distributionsAD27;
    //! \brief stores a field of concentration values
    real *Conc, *Conc_Full;
    //! \brief stores the diffusivity
    real diffusivity;
    //! \brief stores the value for omega (for the diffusivity)
    real omegaDiffusivity;
    // BC NoSlip
    TempforBoundaryConditions Temp;
    // BC Velocity
    TempVelforBoundaryConditions TempVel;
    // BC Pressure
    TempPressforBoundaryConditions TempPress;
    // Plane Conc
    real *ConcPlaneIn, *ConcPlaneOut1, *ConcPlaneOut2;
    std::vector<double> PlaneConcVectorIn, PlaneConcVectorOut1, PlaneConcVectorOut2;

    // trafo///////////////////
    real mTtoWx, mTtoWy, mTtoWz;
    real cTtoWx, cTtoWy, cTtoWz;

    // MGstrafo////////////////
    real cStartx, cStarty, cStartz;
    real cFx, cFy, cFz;

    // typeOfGridNode (formerly known as "geo") /////////////////////
    int *geo; // DEPRECATED: typeOfGridNode for full matrix (not sparse)
    //! \brief stores the type for every lattice node (f.e. fluid node)
    unsigned int *typeOfGridNode;

    // k///////////////////////
    unsigned int *k; // DEPRECATED: index for full matrix

    // neighbor///////////////////////////////////////////////////////////////
    //! \brief store the neighbors in +X, +Y, +Z, and in diagonal negative direction
    //! \brief this information is important because we use an indirect addressing scheme
    uint *neighborX, *neighborY, *neighborZ, *neighborInverse;

    // coordinates////////////////////////////////////////////////////////////
    //! \brief store the coordinates for every lattice node 
    real *coordinateX, *coordinateY, *coordinateZ;

    // body forces////////////
    real *forceX_SP, *forceY_SP, *forceZ_SP;

    // vel parab///////////////
    real *vParab;

    // turbulent viscosity ///
    real *turbViscosity;
    real *gSij, *gSDij, *gDxvx, *gDyvx, *gDzvx, *gDxvy, *gDyvy, *gDzvy, *gDxvz, *gDyvz, *gDzvz; // DebugInformation

    // turbulence intensity //
    real *vx_mean, *vy_mean, *vz_mean;       // means
    real *vxx, *vyy, *vzz, *vxy, *vxz, *vyz; // fluctuations
    std::vector<real> turbulenceIntensity;

    // macroscopic values//////
    // real *vx, *vy, *vz, *rho;  // DEPRECATED: macroscopic values for full matrix
    //! \brief store the macroscopic values (velocity, density, pressure)
    //! \brief for every lattice node
    real *velocityX, *velocityY, *velocityZ, *rho, *pressure;
    //! \brief stores the value for omega
    real omega;
    //! \brief stores the value for viscosity (on level 0)
    real vis;

    // derivations for iso test
    real *dxxUx, *dyyUy, *dzzUz;

    // median-macro-values/////
    real *vx_SP_Med, *vy_SP_Med, *vz_SP_Med, *rho_SP_Med, *press_SP_Med;
    real *vx_SP_Med_Out, *vy_SP_Med_Out, *vz_SP_Med_Out, *rho_SP_Med_Out, *press_SP_Med_Out;
    // Advection-Diffusion
    real *Conc_Med, *Conc_Med_Out;

    // grid////////////////////
    unsigned int nx, ny, nz;
    unsigned int gridNX, gridNY, gridNZ;

    // size of matrix//////////
    unsigned int size_Mat;
    unsigned int sizePlaneXY, sizePlaneYZ, sizePlaneXZ;

    // size of sparse matrix//////////
    //! \brief stores the number of nodes (based on indirect addressing scheme)
    unsigned int numberOfNodes;
    unsigned int size_Array_SP;

    // size of Plane btw. 2 GPUs//////
    unsigned int sizePlaneSB, sizePlaneRB, startB, endB;
    unsigned int sizePlaneST, sizePlaneRT, startT, endT;
    bool isSetSendB, isSetRecvB, isSetSendT, isSetRecvT;
    int *SendT, *SendB, *RecvT, *RecvB;

    // size of Plane for PressMess
    unsigned int sizePlanePress, startP;
    unsigned int sizePlanePressIN, startPIN;
    unsigned int sizePlanePressOUT, startPOUT;
    bool isSetPress;

    // memsizeSP/////////////////
    //! \brief stores the size of the memory consumption for real/int values of the arrays (e.g. coordinates, velocity)
    unsigned int mem_size_real_SP;
    unsigned int mem_size_int_SP;

    // memsize/////////////////
    unsigned int mem_size_real;
    unsigned int mem_size_int;
    unsigned int mem_size_bool;
    unsigned int mem_size_real_yz;

    // print///////////////////
    unsigned int startz, endz;
    real Lx, Ly, Lz, dx;
    real distX, distY, distZ;

    // interface////////////////
    bool need_interface[6];
    unsigned int XdistKn, YdistKn, ZdistKn;
    InterpolationCellCF intCF;
    InterpolationCellFC intFC;
    unsigned int K_CF;
    unsigned int K_FC;
    unsigned int mem_size_kCF;
    unsigned int mem_size_kFC;

    InterpolationCellFC intFCBorder;
    InterpolationCellFC intFCBulk;
    InterpolationCellCF intCFBorder;
    InterpolationCellCF intCFBulk;

    // offset//////////////////
    OffsetCF offCF;
    OffsetCF offCFBulk;
    OffsetFC offFC;
    unsigned int mem_size_kCF_off;
    unsigned int mem_size_kFC_off;

    // BC's////////////////////
    //! \brief stores the boundary condition data
    QforBoundaryConditions noSlipBC, velocityBC, outflowBC, slipBC, stressBC, pressureBC;
    //! \brief number of lattice nodes for the boundary conditions
    unsigned int numberOfNoSlipBCnodesRead, numberOfVeloBCnodesRead, numberOfOutflowBCnodesRead, numberOfSlipBCnodesRead, numberOfStressBCnodesRead, numberOfPressureBCnodesRead;

    QforBoundaryConditions QpressX0, QpressX1, QpressY0, QpressY1, QpressZ0, QpressZ1; // DEPRECATED
    QforBoundaryConditions propellerBC;
    QforBoundaryConditions geometryBC;
    QforBoundaryConditions geometryBCnormalX, geometryBCnormalY, geometryBCnormalZ;
    QforBoundaryConditions inflowBCnormalX, inflowBCnormalY, inflowBCnormalZ;
    QforBoundaryConditions outflowBCnormalX, outflowBCnormalY, outflowBCnormalZ;
    QforBoundaryConditions QInlet, QOutlet, QPeriodic; // DEPRECATED
    unsigned int kInletQread, kOutletQread;  // DEPRECATED

    WallModelParameters wallModel;

    // testRoundoffError
    Distributions27 kDistTestRE;

    //////////////////////////////////////////////////////////////////////////
    // velocities to fit the force
    real *VxForce, *VyForce, *VzForce;
    //////////////////////////////////////////////////////////////////////////
    //! \brief sets the forcing uniform on every fluid node in all three space dimensions
    real *forcing;

    // Measure Points/////////
    std::vector<MeasurePoints> MP;
    unsigned int *kMP;
    real *VxMP;
    real *VyMP;
    real *VzMP;
    real *RhoMP;
    unsigned int memSizerealkMP, memSizeIntkMP, numberOfPointskMP;
    unsigned int numberOfValuesMP;

    // Drag Lift//////////////
    double *DragPreX, *DragPostX;
    double *DragPreY, *DragPostY;
    double *DragPreZ, *DragPostZ;
    std::vector<double> DragXvector;
    std::vector<double> DragYvector;
    std::vector<double> DragZvector;

    // 2ndMoments////////////
    real *kxyFromfcNEQ, *kyzFromfcNEQ, *kxzFromfcNEQ, *kxxMyyFromfcNEQ, *kxxMzzFromfcNEQ;

    // 3rdMoments////////////
    real *CUMbbb, *CUMabc, *CUMbac, *CUMbca, *CUMcba, *CUMacb, *CUMcab;

    // HigherMoments/////////
    real *CUMcbb, *CUMbcb, *CUMbbc, *CUMcca, *CUMcac, *CUMacc, *CUMbcc, *CUMcbc, *CUMccb, *CUMccc;

    // CpTop/////////////////
    int *cpTopIndex;
    double *cpPressTop;
    unsigned int numberOfPointsCpTop;
    std::vector<std::vector<double>> cpTop;
    std::vector<double> pressMirror;
    std::vector<bool> isOutsideInterface;
    unsigned int numberOfPointsPressWindow;

    // CpBottom/////////////
    int *cpBottomIndex;
    double *cpPressBottom;
    unsigned int numberOfPointsCpBottom;
    std::vector<std::vector<double>> cpBottom;

    // CpBottom2////////////
    int *cpBottom2Index;
    double *cpPressBottom2;
    unsigned int numberOfPointsCpBottom2;
    std::vector<std::vector<double>> cpBottom2;

    // Concentration////////
    int *concIndex;
    real *concentration;
    unsigned int numberOfPointsConc;

    // street X and Y velocity fractions///////
    real *streetFractionXvelocity;
    real *streetFractionYvelocity;
    int *naschVelocity;
    uint numberOfStreetNodes;

    // deltaPhi
    real deltaPhi;

    ////////////////////////////////////////////////////////////////////////////
    // particles
    PathLineParticles plp;
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // 1D domain decomposition
    std::vector<ProcessNeighbor27> sendProcessNeighbor;
    std::vector<ProcessNeighbor27> recvProcessNeighbor;
    ///////////////////////////////////////////////////////
    // 3D domain decomposition
    std::vector<ProcessNeighbor27> sendProcessNeighborX;
    std::vector<ProcessNeighbor27> sendProcessNeighborY;
    std::vector<ProcessNeighbor27> sendProcessNeighborZ;
    std::vector<ProcessNeighbor27> recvProcessNeighborX;
    std::vector<ProcessNeighbor27> recvProcessNeighborY;
    std::vector<ProcessNeighbor27> recvProcessNeighborZ;

    std::vector<ProcessNeighbor27> sendProcessNeighborsAfterFtoCX;
    std::vector<ProcessNeighbor27> sendProcessNeighborsAfterFtoCY;
    std::vector<ProcessNeighbor27> sendProcessNeighborsAfterFtoCZ;
    std::vector<ProcessNeighbor27> recvProcessNeighborsAfterFtoCX;
    std::vector<ProcessNeighbor27> recvProcessNeighborsAfterFtoCY;
    std::vector<ProcessNeighbor27> recvProcessNeighborsAfterFtoCZ;
    ///////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    std::vector<ProcessNeighbor27> sendProcessNeighborADX;
    std::vector<ProcessNeighbor27> sendProcessNeighborADY;
    std::vector<ProcessNeighbor27> sendProcessNeighborADZ;
    std::vector<ProcessNeighbor27> recvProcessNeighborADX;
    std::vector<ProcessNeighbor27> recvProcessNeighborADY;
    std::vector<ProcessNeighbor27> recvProcessNeighborADZ;
    ///////////////////////////////////////////////////////
    // 3D domain decomposition F3
    std::vector<ProcessNeighborF3> sendProcessNeighborF3X;
    std::vector<ProcessNeighborF3> sendProcessNeighborF3Y;
    std::vector<ProcessNeighborF3> sendProcessNeighborF3Z;
    std::vector<ProcessNeighborF3> recvProcessNeighborF3X;
    std::vector<ProcessNeighborF3> recvProcessNeighborF3Y;
    std::vector<ProcessNeighborF3> recvProcessNeighborF3Z;
    ////////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition: position (index in array) of corner nodes in ProcessNeighbor27
    struct EdgeNodePositions {
        int indexOfProcessNeighborRecv;
        int indexInRecvBuffer;
        int indexOfProcessNeighborSend;
        int indexInSendBuffer;
        EdgeNodePositions(int indexOfProcessNeighborRecv, int indexInRecvBuffer, int indexOfProcessNeighborSend,
                          int indexInSendBuffer)
            : indexOfProcessNeighborRecv(indexOfProcessNeighborRecv), indexInRecvBuffer(indexInRecvBuffer),
              indexOfProcessNeighborSend(indexOfProcessNeighborSend), indexInSendBuffer(indexInSendBuffer)
        {
        }
    };
    std::vector<EdgeNodePositions> edgeNodesXtoY;
    std::vector<EdgeNodePositions> edgeNodesXtoZ;
    std::vector<EdgeNodePositions> edgeNodesYtoZ;

    ///////////////////////////////////////////////////////
    uint *fluidNodeIndices;
    uint numberOfFluidNodes;
    uint *fluidNodeIndicesBorder;
    uint numberOfFluidNodesBorder;
};

//! \brief Class for LBM-parameter management
class VIRTUALFLUIDS_GPU_EXPORT Parameter
{
public:
    explicit Parameter(const vf::basics::ConfigurationFile &configData, const int numberOfProcesses = 1, const int myId = 0);
    Parameter(const int numberOfProcesses = 1, const int myId = 0);
    ~Parameter();
    void initLBMSimulationParameter();

    //! \brief Pointer to instance of LBMSimulationParameter - stored on Host System
    std::shared_ptr<LBMSimulationParameter> getParH(int level);
    //! \brief Pointer to instance of LBMSimulationParameter - stored on Device (GPU)
    std::shared_ptr<LBMSimulationParameter> getParD(int level);

    void copyMeasurePointsArrayToVector(int lev);

    //////////////////////////////////////////////////////////////////////////
    // setter
    void setForcing(real forcingX, real forcingY, real forcingZ);
    void setQuadricLimiters(real quadricLimiterP, real quadricLimiterM, real quadricLimiterD);
    void setPhi(real inPhi);
    void setAngularVelocity(real inAngVel);
    void setStepEnsight(unsigned int step);
    void setOutputCount(unsigned int outputCount);
    void setlimitOfNodesForVTK(unsigned int limitOfNodesForVTK);
    void setStartTurn(unsigned int inStartTurn);
    void setDiffOn(bool isDiff);
    void setCompOn(bool isComp);
    void setDiffMod(int DiffMod);
    void setDiffusivity(real Diffusivity);
    void setD3Qxx(int d3qxx);
    void setMaxLevel(int maxlevel);
    void setParticleBasicLevel(int pbl);
    void setParticleInitLevel(int pil);
    void setNumberOfParticles(int nop);
    void setCalcParticles(bool calcParticles);
    void setStartXHotWall(real startXHotWall);
    void setEndXHotWall(real endXHotWall);
    void setTimestepEnd(unsigned int tend);
    void setTimestepOut(unsigned int tout);
    void setTimestepStartOut(unsigned int tStartOut);
    void setTimestepOfCoarseLevel(unsigned int timestep);
    void setCalcTurbulenceIntensity(bool calcVelocityAndFluctuations);
    void setCalcMedian(bool calcMedian);
    void setCalcDragLift(bool calcDragLift);
    void setCalcCp(bool calcCp);
    void setWriteVeloASCIIfiles(bool writeVeloASCII);
    void setCalcPlaneConc(bool calcPlaneConc);
    void setTimeCalcMedStart(int CalcMedStart);
    void setTimeCalcMedEnd(int CalcMedEnd);
    void setMaxDev(int maxdev);
    void setMyID(int myid);
    void setNumprocs(int numprocs);
    void setPressInID(unsigned int PressInID);
    void setPressOutID(unsigned int PressOutID);
    void setPressInZ(unsigned int PressInZ);
    void setPressOutZ(unsigned int PressOutZ);
    void settimestepForMP(unsigned int timestepForMP);
    void setOutputPath(std::string oPath);
    void setOutputPrefix(std::string oPrefix);
    void setPathAndFilename(std::string fname);
    void setGridPath(std::string gridPath);
    void setGeometryFileC(std::string GeometryFileC);
    void setGeometryFileM(std::string GeometryFileM);
    void setGeometryFileF(std::string GeometryFileF);
    void setkFull(std::string kFull);
    void setgeoFull(std::string geoFull);
    void setgeoVec(std::string geoVec);
    void setcoordX(std::string coordX);
    void setcoordY(std::string coordY);
    void setcoordZ(std::string coordZ);
    void setneighborX(std::string neighborX);
    void setneighborY(std::string neighborY);
    void setneighborZ(std::string neighborZ);
    void setneighborWSB(std::string neighborWSB);
    void setscaleCFC(std::string scaleCFC);
    void setscaleCFF(std::string scaleCFF);
    void setscaleFCC(std::string scaleFCC);
    void setscaleFCF(std::string scaleFCF);
    void setscaleOffsetCF(std::string scaleOffsetCF);
    void setscaleOffsetFC(std::string scaleOffsetFC);
    void setgeomBoundaryBcQs(std::string geomBoundaryBcQs);
    void setgeomBoundaryBcValues(std::string geomBoundaryBcValues);
    void setnoSlipBcPos(std::string noSlipBcPos);
    void setnoSlipBcQs(std::string noSlipBcQs);
    void setnoSlipBcValue(std::string noSlipBcValue);
    void setnoSlipBcValues(std::string noSlipBcValues);
    void setslipBcPos(std::string slipBcPos);
    void setslipBcQs(std::string slipBcQs);
    void setslipBcValue(std::string slipBcValue);
    void setpressBcPos(std::string pressBcPos);
    void setpressBcQs(std::string pressBcQs);
    void setpressBcValue(std::string pressBcValue);
    void setpressBcValues(std::string pressBcValues);
    void setvelBcQs(std::string velBcQs);
    void setvelBcValues(std::string velBcValues);
    void setinletBcQs(std::string inletBcQs);
    void setinletBcValues(std::string inletBcValues);
    void setoutletBcQs(std::string outletBcQs);
    void setoutletBcValues(std::string outletBcValues);
    void settopBcQs(std::string topBcQs);
    void settopBcValues(std::string topBcValues);
    void setbottomBcQs(std::string bottomBcQs);
    void setbottomBcValues(std::string bottomBcValues);
    void setfrontBcQs(std::string frontBcQs);
    void setfrontBcValues(std::string frontBcValues);
    void setbackBcQs(std::string backBcQs);
    void setbackBcValues(std::string backBcValues);
    void setwallBcQs(std::string wallBcQs);
    void setwallBcValues(std::string wallBcValues);
    void setperiodicBcQs(std::string periodicBcQs);
    void setperiodicBcValues(std::string periodicBcValues);
    void setpropellerCylinder(std::string propellerCylinder);
    void setpropellerValues(std::string propellerValues);
    void setpropellerQs(std::string propellerQs);
    void setmeasurePoints(std::string measurePoints);
    void setnumberNodes(std::string numberNodes);
    void setLBMvsSI(std::string LBMvsSI);
    void setcpTop(std::string cpTop);
    void setcpBottom(std::string cpBottom);
    void setcpBottom2(std::string cpBottom2);
    void setConcentration(std::string concFile);
    void setStreetVelocity(std::string streetVelocity);
    void setPrintFiles(bool printfiles);
    void setReadGeo(bool readGeo);
    void setTemperatureInit(real Temp);
    void setTemperatureBC(real TempBC);
    void setViscosityLB(real Viscosity);
    void setVelocityLB(real Velocity);
    void setViscosityRatio(real ViscosityRatio);
    void setVelocityRatio(real VelocityRatio);
    void setDensityRatio(real DensityRatio);
    void setPressRatio(real PressRatio);
    void setRealX(real RealX);
    void setRealY(real RealY);
    void setRe(real Re);
    void setFactorPressBC(real factorPressBC);
    void setIsGeo(bool isGeo);
    void setIsGeoNormal(bool isGeoNormal);
    void setIsInflowNormal(bool isInflowNormal);
    void setIsOutflowNormal(bool isOutflowNormal);
    void setIsProp(bool isProp);
    void setIsCp(bool isCp);
    void setConcFile(bool concFile);
    void setStreetVelocityFile(bool streetVelocityFile);
    void setUseMeasurePoints(bool useMeasurePoints);
    void setUseWale(bool useWale);
    void setUseTurbulentViscosity(bool useTurbulentViscosity);
    void setUseAMD(bool useAMD);
    void setSGSConstant(real SGSConstant);
    void setHasWallModelMonitor(bool hasWallModelMonitor);
    void setUseInitNeq(bool useInitNeq);
    void setSimulatePorousMedia(bool simulatePorousMedia);
    void setIsF3(bool isF3);
    void setIsBodyForce(bool isBodyForce);
    void setclockCycleForMP(real clockCycleForMP);
    void setDevices(std::vector<uint> devices);
    void setGridX(std::vector<int> GridX);
    void setGridY(std::vector<int> GridY);
    void setGridZ(std::vector<int> GridZ);
    void setDistX(std::vector<int> DistX);
    void setDistY(std::vector<int> DistY);
    void setDistZ(std::vector<int> DistZ);
    void setScaleLBMtoSI(std::vector<real> scaleLBMtoSI);
    void setTranslateLBMtoSI(std::vector<real> translateLBMtoSI);
    void setMinCoordX(std::vector<real> MinCoordX);
    void setMinCoordY(std::vector<real> MinCoordY);
    void setMinCoordZ(std::vector<real> MinCoordZ);
    void setMaxCoordX(std::vector<real> MaxCoordX);
    void setMaxCoordY(std::vector<real> MaxCoordY);
    void setMaxCoordZ(std::vector<real> MaxCoordZ);
    void setTempH(TempforBoundaryConditions *TempH);
    void setTempD(TempforBoundaryConditions *TempD);
    void setTempVelH(TempVelforBoundaryConditions *TempVelH);
    void setTempVelD(TempVelforBoundaryConditions *TempVelD);
    void setTempPressH(TempPressforBoundaryConditions *TempPressH);
    void setTempPressD(TempPressforBoundaryConditions *TempPressD);
    void setTimeDoCheckPoint(unsigned int tDoCheckPoint);
    void setTimeDoRestart(unsigned int tDoRestart);
    void setDoCheckPoint(bool doCheckPoint);
    void setDoRestart(bool doRestart);
    void setObj(std::string str, bool isObj);
    void setUseGeometryValues(bool GeometryValues);
    void setCalc2ndOrderMoments(bool is2ndOrderMoments);
    void setCalc3rdOrderMoments(bool is3rdOrderMoments);
    void setCalcHighOrderMoments(bool isHighOrderMoments);
    void setMemsizeGPU(double admem, bool reset);
    // 1D domain decomposition
    void setPossNeighborFiles(std::vector<std::string> possNeighborFiles, std::string sor);
    void setNumberOfProcessNeighbors(unsigned int numberOfProcessNeighbors, int level, std::string sor);
    void setIsNeighbor(bool isNeighbor);
    // 3D domain decomposition
    void setPossNeighborFilesX(std::vector<std::string> possNeighborFiles, std::string sor);
    void setPossNeighborFilesY(std::vector<std::string> possNeighborFiles, std::string sor);
    void setPossNeighborFilesZ(std::vector<std::string> possNeighborFiles, std::string sor);
    void setNumberOfProcessNeighborsX(unsigned int numberOfProcessNeighbors, int level, std::string sor);
    void setNumberOfProcessNeighborsY(unsigned int numberOfProcessNeighbors, int level, std::string sor);
    void setNumberOfProcessNeighborsZ(unsigned int numberOfProcessNeighbors, int level, std::string sor);
    void setIsNeighborX(bool isNeighbor);
    void setIsNeighborY(bool isNeighbor);
    void setIsNeighborZ(bool isNeighbor);
    void setSendProcessNeighborsAfterFtoCX(int numberOfNodes, int level, int arrayIndex);
    void setSendProcessNeighborsAfterFtoCY(int numberOfNodes, int level, int arrayIndex);
    void setSendProcessNeighborsAfterFtoCZ(int numberOfNodes, int level, int arrayIndex);
    void setRecvProcessNeighborsAfterFtoCX(int numberOfNodes, int level, int arrayIndex);
    void setRecvProcessNeighborsAfterFtoCY(int numberOfNodes, int level, int arrayIndex);
    void setRecvProcessNeighborsAfterFtoCZ(int numberOfNodes, int level, int arrayIndex);
    // void setQinflowH(QforBoundaryConditions* QinflowH);
    // void setQinflowD(QforBoundaryConditions* QinflowD);
    // void setQoutflowH(QforBoundaryConditions* QoutflowH);
    // void setQoutflowD(QforBoundaryConditions* QoutflowD);
    // Normals
    void setgeomBoundaryNormalX(std::string geomNormalX);
    void setgeomBoundaryNormalY(std::string geomNormalY);
    void setgeomBoundaryNormalZ(std::string geomNormalZ);
    void setInflowBoundaryNormalX(std::string inflowNormalX);
    void setInflowBoundaryNormalY(std::string inflowNormalY);
    void setInflowBoundaryNormalZ(std::string inflowNormalZ);
    void setOutflowBoundaryNormalX(std::string outflowNormalX);
    void setOutflowBoundaryNormalY(std::string outflowNormalY);
    void setOutflowBoundaryNormalZ(std::string outflowNormalZ);
    // Kernel
    void setMainKernel(std::string kernel);
    void setMultiKernelOn(bool isOn);
    void setMultiKernelLevel(std::vector<int> kernelLevel);
    void setMultiKernel(std::vector<std::string> kernel);

    void setADKernel(std::string adKernel);

    // adder
    void addActuator(SPtr<PreCollisionInteractor> actuator);
    void addProbe(SPtr<PreCollisionInteractor> probes);

    // getter
    double *getForcesDouble();
    real *getForcesHost();
    real *getForcesDev();
    double *getQuadricLimitersDouble();
    real *getQuadricLimitersHost();
    real *getQuadricLimitersDev();
    real getPhi();
    real getAngularVelocity();
    real getStartXHotWall();
    real getEndXHotWall();
    unsigned int getStepEnsight();
    unsigned int getOutputCount();
    unsigned int getlimitOfNodesForVTK();
    unsigned int getStartTurn();
    bool getEvenOrOdd(int level);
    bool getDiffOn();
    bool getCompOn();
    bool getPrintFiles();
    bool getReadGeo();
    bool getCalcTurbulenceIntensity();
    bool getCalcMedian();
    bool getCalcDragLift();
    bool getCalcCp();
    bool getCalcParticles();
    bool getWriteVeloASCIIfiles();
    bool getCalcPlaneConc();
    //! \returns index of finest level
    int getFine();
    //! \returns index of coarsest level
    int getCoarse();
    int getParticleBasicLevel();
    int getParticleInitLevel();
    int getNumberOfParticles();
    int getDiffMod();
    int getFactorNZ();
    int getD3Qxx();
    //! \returns the maximum level of grid refinement
    int getMaxLevel();
    int getTimeCalcMedStart();
    int getTimeCalcMedEnd();
    int getMaxDev();
    //! \returns the ID of the current MPI process
    int getMyProcessID();
    int getNumprocs();
    std::string getOutputPath();
    std::string getOutputPrefix();
    std::string getFName();
    std::string getGridPath();
    std::string getGeometryFileC();
    std::string getGeometryFileM();
    std::string getGeometryFileF();
    std::string getkFull();
    std::string getgeoFull();
    std::string getgeoVec();
    std::string getcoordX();
    std::string getcoordY();
    std::string getcoordZ();
    std::string getneighborX();
    std::string getneighborY();
    std::string getneighborZ();
    std::string getneighborWSB();
    std::string getscaleCFC();
    std::string getscaleCFF();
    std::string getscaleFCC();
    std::string getscaleFCF();
    std::string getscaleOffsetCF();
    std::string getscaleOffsetFC();
    std::string getgeomBoundaryBcQs();
    std::string getgeomBoundaryBcValues();
    std::string getnoSlipBcPos();
    std::string getnoSlipBcQs();
    std::string getnoSlipBcValue();
    std::string getnoSlipBcValues();
    std::string getslipBcPos();
    std::string getslipBcQs();
    std::string getslipBcValue();
    std::string getpressBcPos();
    std::string getpressBcQs();
    std::string getpressBcValue();
    std::string getpressBcValues();
    std::string getvelBcQs();
    std::string getvelBcValues();
    std::string getinletBcQs();
    std::string getinletBcValues();
    std::string getoutletBcQs();
    std::string getoutletBcValues();
    std::string gettopBcQs();
    std::string gettopBcValues();
    std::string getbottomBcQs();
    std::string getbottomBcValues();
    std::string getfrontBcQs();
    std::string getfrontBcValues();
    std::string getbackBcQs();
    std::string getbackBcValues();
    std::string getwallBcQs();
    std::string getwallBcValues();
    std::string getperiodicBcQs();
    std::string getperiodicBcValues();
    std::string getpropellerQs();
    std::string getpropellerCylinder();
    std::string getpropellerValues();
    std::string getmeasurePoints();
    std::string getnumberNodes();
    std::string getLBMvsSI();
    std::string getcpTop();
    std::string getcpBottom();
    std::string getcpBottom2();
    std::string getConcentration();
    std::string getStreetVelocityFilePath();
    unsigned int getPressInID();
    unsigned int getPressOutID();
    unsigned int getPressInZ();
    unsigned int getPressOutZ();
    unsigned int getMemSizereal(int level);
    unsigned int getMemSizeInt(int level);
    unsigned int getMemSizeBool(int level);
    unsigned int getMemSizerealYZ(int level);
    unsigned int getSizeMat(int level);
    unsigned int getTimestepStart();
    unsigned int getTimestepInit();
    unsigned int getTimestepEnd();
    unsigned int getTimestepOut();
    unsigned int getTimestepStartOut();
    unsigned int getTimestepForMP();
    unsigned int getTimestepOfCoarseLevel();
    real getDiffusivity();
    real getTemperatureInit();
    real getTemperatureBC();
    real getViscosity();
    real getVelocity();
    //! \returns the viscosity ratio in SI/LB units
    real getViscosityRatio();
    //! \returns the velocity ratio in SI/LB units
    real getVelocityRatio();
    //! \returns the density ratio in SI/LB units
    real getDensityRatio();
    //! \returns the pressure ratio in SI/LB units
    real getPressureRatio();
    //! \returns the time ratio in SI/LB units
    real getTimeRatio();
    //! \returns the length ratio in SI/LB units
    real getLengthRatio();
    //! \returns the force ratio in SI/LB units
    real getForceRatio();
    real getRealX();
    real getRealY();
    real getRe();
    real getFactorPressBC();
    real getclockCycleForMP();
    std::vector<uint> getDevices();
    std::vector<int> getGridX();
    std::vector<int> getGridY();
    std::vector<int> getGridZ();
    std::vector<int> getDistX();
    std::vector<int> getDistY();
    std::vector<int> getDistZ();
    std::vector<real> getScaleLBMtoSI();
    std::vector<real> getTranslateLBMtoSI();
    std::vector<real> getMinCoordX();
    std::vector<real> getMinCoordY();
    std::vector<real> getMinCoordZ();
    std::vector<real> getMaxCoordX();
    std::vector<real> getMaxCoordY();
    std::vector<real> getMaxCoordZ();
    TempforBoundaryConditions *getTempH();
    TempforBoundaryConditions *getTempD();
    TempVelforBoundaryConditions *getTempVelH();
    TempVelforBoundaryConditions *getTempVelD();
    TempPressforBoundaryConditions *getTempPressH();
    TempPressforBoundaryConditions *getTempPressD();
    std::vector<SPtr<PreCollisionInteractor>> getActuators();
    //! \returns the probes, e.g. point or plane probe
    std::vector<SPtr<PreCollisionInteractor>> getProbes();
    unsigned int getTimeDoCheckPoint();
    unsigned int getTimeDoRestart();
    bool getDoCheckPoint();
    bool getDoRestart();
    bool overWritingRestart(unsigned int t);
    bool getIsGeo();
    bool getIsGeoNormal();
    bool getIsInflowNormal();
    bool getIsOutflowNormal();
    bool getIsProp();
    bool getIsCp();
    bool getIsGeometryValues();
    bool getCalc2ndOrderMoments();
    bool getCalc3rdOrderMoments();
    bool getCalcHighOrderMoments();
    bool getConcFile();
    bool isStreetVelocityFile();
    bool getUseMeasurePoints();
    bool getUseWale();
    bool getUseTurbulentViscosity();
    bool getUseAMD();
    real getSGSConstant();
    bool getHasWallModelMonitor();
    bool getUseInitNeq();
    bool getSimulatePorousMedia();
    bool getIsF3();
    bool getIsBodyForce();
    double getMemsizeGPU();
    // 1D domain decomposition
    std::vector<std::string> getPossNeighborFiles(std::string sor);
    unsigned int getNumberOfProcessNeighbors(int level, std::string sor);
    bool getIsNeighbor();
    // 3D domain decomposition
    std::vector<std::string> getPossNeighborFilesX(std::string sor);
    std::vector<std::string> getPossNeighborFilesY(std::string sor);
    std::vector<std::string> getPossNeighborFilesZ(std::string sor);
    unsigned int getNumberOfProcessNeighborsX(int level, std::string sor);
    unsigned int getNumberOfProcessNeighborsY(int level, std::string sor);
    unsigned int getNumberOfProcessNeighborsZ(int level, std::string sor);
    bool getIsNeighborX();
    bool getIsNeighborY();
    bool getIsNeighborZ();
    // Normals
    std::string getgeomBoundaryNormalX();
    std::string getgeomBoundaryNormalY();
    std::string getgeomBoundaryNormalZ();
    std::string getInflowBoundaryNormalX();
    std::string getInflowBoundaryNormalY();
    std::string getInflowBoundaryNormalZ();
    std::string getOutflowBoundaryNormalX();
    std::string getOutflowBoundaryNormalY();
    std::string getOutflowBoundaryNormalZ();
    // CUDA random number
    curandState *getRandomState();
    // Kernel
    std::string getMainKernel();
    bool getMultiKernelOn();
    std::vector<int> getMultiKernelLevel();
    std::vector<std::string> getMultiKernel();

    std::string getADKernel();

    // Forcing///////////////
    real *forcingH, *forcingD;
    double hostForcing[3];

    //////////////////////////////////////////////////////////////////////////
    // limiters
    real *quadricLimitersH, *quadricLimitersD;
    double hostQuadricLimiters[3];

    ////////////////////////////////////////////////////////////////////////////
    // initial condition fluid
    void setInitialCondition(std::function<void(real, real, real, real &, real &, real &, real &)> initialCondition);
    std::function<void(real, real, real, real &, real &, real &, real &)> &getInitialCondition();

    // initial condition concentration
    void setInitialConditionAD(std::function<void(real, real, real, real &)> initialConditionAD);
    std::function<void(real, real, real, real &)> &getInitialConditionAD();

    ////////////////////////////////////////////////////////////////////////////

    std::vector<std::shared_ptr<LBMSimulationParameter>> parH = std::vector<std::shared_ptr<LBMSimulationParameter>>(1);
    std::vector<std::shared_ptr<LBMSimulationParameter>> parD = std::vector<std::shared_ptr<LBMSimulationParameter>>(1);

    ////////////////////////////////////////////////////////////////////////////

private:
    void readConfigData(const vf::basics::ConfigurationFile &configData);
    void initGridPaths();
    void initGridBasePoints();
    void initDefaultLBMkernelAllLevels();

private:
    bool compOn{ false };
    bool diffOn{ false };
    bool isF3{ false };
    bool calcDragLift{ false };
    bool calcCp{ false };
    bool writeVeloASCII{ false };
    bool calcPlaneConc{ false };
    bool calcVelocityAndFluctuations{ false };
    bool isBodyForce{ false };
    int diffMod{ 27 };
    //! \property maximum level of grid refinement
    int maxlevel{ 0 };
    int coarse{ 0 };
    int fine{ 0 };
    int factor_gridNZ{ 2 };
    int D3Qxx{ 27 };
    InitCondition ic;
    double memsizeGPU;
    unsigned int limitOfNodesForVTK;
    unsigned int outputCount;
    unsigned int timestep;

    // Kernel
    std::string mainKernel{ "CumulantK17CompChim" };
    bool multiKernelOn{ false };
    std::vector<int> multiKernelLevel;
    std::vector<std::string> multiKernel;
    bool kernelNeedsFluidNodeIndicesToRun = false;
    std::string adKernel;

    //////////////////////////////////////////////////////////////////////////
    // particles
    int particleBasicLevel{ 0 };
    int particleInitLevel{ 0 };
    int numberOfParticles{ 0 };
    bool calcParticles{ false };
    real startXHotWall{ (real)0.0 };
    real endXHotWall{ (real)0.0 };
    //////////////////////////////////////////////////////////////////////////
    // CUDA random number generation
    curandState *devState;
    //////////////////////////////////////////////////////////////////////////

    // Temperature
    TempforBoundaryConditions *TempH, *TempD;
    // Temperature Velocity
    TempVelforBoundaryConditions *TempVelH, *TempVelD;
    // Temperature Pressure
    TempPressforBoundaryConditions *TempPressH, *TempPressD;

    // Drehung///////////////
    real Phi{ 0.0 };
    real angularVelocity;
    unsigned int startTurn;

    // PreCollisionInteractors //////////////
    std::vector<SPtr<PreCollisionInteractor>> actuators;
    std::vector<SPtr<PreCollisionInteractor>> probes;

    // Step of Ensight writing//
    unsigned int stepEnsight;

    real TrafoXtoWorld(int CoordX, int level);
    real TrafoYtoWorld(int CoordY, int level);
    real TrafoZtoWorld(int CoordZ, int level);

public:
    real TrafoXtoMGsWorld(int CoordX, int level);
    real TrafoYtoMGsWorld(int CoordY, int level);
    real TrafoZtoMGsWorld(int CoordZ, int level);

private:
    // Multi GPGPU///////////////
    // 1D domain decomposition
    std::vector<std::string> possNeighborFilesSend;
    std::vector<std::string> possNeighborFilesRecv;
    bool isNeigbor;
    // 3D domain decomposition
    std::vector<std::string> possNeighborFilesSendX, possNeighborFilesSendY, possNeighborFilesSendZ;
    std::vector<std::string> possNeighborFilesRecvX, possNeighborFilesRecvY, possNeighborFilesRecvZ;
    bool isNeigborX, isNeigborY, isNeigborZ;

    ////////////////////////////////////////////////////////////////////////////
    //! \brief initial condition fluid
    std::function<void(real, real, real, real &, real &, real &, real &)> initialCondition;
    //! \brief initial condition concentration
    std::function<void(real, real, real, real &)> initialConditionAD;

    ////////////////////////////////////////////////////////////////////////////
    // cuda streams

    //! determines whether streams and thus communication hiding should be used
    bool useStreams{ false };
    std::unique_ptr<CudaStreamManager> cudaStreamManager;

public:
    //! \brief sets whether streams and thus communication hiding should be used
    //! \details This function is only useful for simulations on multiple GPUs. If there is only one MPI process, the
    //! passed value is automatically overwritten with false.
    void setUseStreams(bool useStreams);
    bool getUseStreams();
    std::unique_ptr<CudaStreamManager> &getStreamManager();
    bool getKernelNeedsFluidNodeIndicesToRun();
    void setKernelNeedsFluidNodeIndicesToRun(bool  kernelNeedsFluidNodeIndicesToRun);

    void initProcessNeighborsAfterFtoCX(int level);
    void initProcessNeighborsAfterFtoCY(int level);
    void initProcessNeighborsAfterFtoCZ(int level);

    bool useReducedCommunicationAfterFtoC{ true };
};

#endif

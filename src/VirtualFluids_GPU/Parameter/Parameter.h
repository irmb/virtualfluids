#ifndef PARAMETER_H
#define PARAMETER_H

#include <vector>
#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include "Calculation/PorousMedia.h"
//#include "Output/LogWriter.hpp"
//#include "boost/serialization/serialization.hpp"
//#include "boost/serialization/vector.hpp"

#include <cuda_runtime.h>
#include <helper_cuda.h>
//random numbers
#include <curand.h>
#include <curand_kernel.h>
#include "core/PointerDefinitions.h"
#include "VirtualFluidsDefinitions.h"

//struct
struct ParameterStruct{
	bool evenOrOdd;
	unsigned int numberofthreads;

	//distributions///////////
	//Distributions19 d0;
	Distributions27 d0;
	Distributions27 d0SP;

	//distributions F3////////
	Distributions6 g6;

	//thermo//////////////////
	Distributions7 d7;
	Distributions27 d27;
	real *Conc, *Conc_Full;
	real diffusivity;
	//BC NoSlip
	TempforBoundaryConditions Temp;
	//BC Velocity
	TempVelforBoundaryConditions TempVel;
	//BC Pressure
	TempPressforBoundaryConditions TempPress;
	//Plane Conc
	real *ConcPlaneIn, *ConcPlaneOut1, *ConcPlaneOut2;
	std::vector<double> PlaneConcVectorIn, PlaneConcVectorOut1, PlaneConcVectorOut2;

	//trafo///////////////////
	real mTtoWx, mTtoWy, mTtoWz;
	real cTtoWx, cTtoWy, cTtoWz;

	//MGstrafo////////////////
	real cStartx, cStarty, cStartz;
	real cFx, cFy, cFz;

	//geo/////////////////////
	int *geo;
	unsigned int *geoSP;

	//k///////////////////////
	unsigned int *k;

	//neighbor////////////////
	//unsigned int *neighborX, *neighborY, *neighborZ;
	unsigned int *neighborX_SP, *neighborY_SP, *neighborZ_SP, *neighborWSB_SP;

	//coordinates////////////
	//unsigned int *coordX_SP, *coordY_SP, *coordZ_SP;
	real *coordX_SP, *coordY_SP, *coordZ_SP;

	//vel parab///////////////
	real *vParab;

	// turbulent viscosity ///
	real *turbViscosity;

	//macroscopic values//////
	real *vx,    *vy,    *vz,    *rho;
	real *vx_SP, *vy_SP, *vz_SP, *rho_SP, *press_SP;
	real vis, omega;

	//median-macro-values/////
	real *vx_SP_Med, *vy_SP_Med, *vz_SP_Med, *rho_SP_Med, *press_SP_Med;
	real *vx_SP_Med_Out, *vy_SP_Med_Out, *vz_SP_Med_Out, *rho_SP_Med_Out, *press_SP_Med_Out;

	//grid////////////////////
	unsigned int nx, ny, nz;
	unsigned int gridNX, gridNY, gridNZ;

	//size of matrix//////////
	unsigned int size_Mat;
	unsigned int sizePlaneXY, sizePlaneYZ, sizePlaneXZ;

	//size of sparse matrix//////////
	unsigned int size_Mat_SP;
	unsigned int size_Array_SP;

	//size of Plane btw. 2 GPUs//////
	unsigned int sizePlaneSB, sizePlaneRB, startB, endB;
	unsigned int sizePlaneST, sizePlaneRT, startT, endT;
	bool isSetSendB, isSetRecvB, isSetSendT, isSetRecvT;
	int *SendT, *SendB, *RecvT, *RecvB;

	//size of Plane for PressMess
	unsigned int sizePlanePress, startP;
	unsigned int sizePlanePressIN, startPIN;
	unsigned int sizePlanePressOUT, startPOUT;
	bool isSetPress;

	//memsizeSP/////////////////
	unsigned int mem_size_real_SP;
	unsigned int mem_size_int_SP;

	//memsize/////////////////
	unsigned int mem_size_real;
	unsigned int mem_size_int;
	unsigned int mem_size_bool;
	unsigned int mem_size_real_yz;

	//print///////////////////
	unsigned int startz, endz;
	real Lx,Ly,Lz,dx;
	real distX, distY, distZ;

	//interface////////////////
	bool need_interface[6];
	unsigned int XdistKn, YdistKn, ZdistKn;
	InterpolationCellCF intCF;
	InterpolationCellFC intFC;
	unsigned int K_CF;
	unsigned int K_FC;
	unsigned int mem_size_kCF;
	unsigned int mem_size_kFC;

	//offset//////////////////
	OffsetCF offCF;
	OffsetFC offFC;
	unsigned int mem_size_kCF_off;
	unsigned int mem_size_kFC_off;

	//BC's////////////////////
	QforBoundaryConditions  QWall,   Qinflow,      Qoutflow,      QSlip;
	unsigned int            kQ,      kInflowQ,     kOutflowQ,     kSlipQ;
	unsigned int            kQread,  kInflowQread, kOutflowQread, kSlipQread;

	QforBoundaryConditions  QpressX0,QpressX1,QpressY0,QpressY1,QpressZ0,QpressZ1;
	QforBoundaryConditions  QPropeller;
	QforBoundaryConditions  QPress;
	QforBoundaryConditions  QGeom;
	QforBoundaryConditions  QGeomNormalX,    QGeomNormalY,    QGeomNormalZ;
	QforBoundaryConditions  QInflowNormalX,  QInflowNormalY,  QInflowNormalZ;
	QforBoundaryConditions  QOutflowNormalX, QOutflowNormalY, QOutflowNormalZ;
	QforBoundaryConditions  QInlet, QOutlet, QPeriodic;
	unsigned int            kInletQread, kOutletQread;
	unsigned int            kPressQ, kPressQread;
	//testRoundoffError
	Distributions27         kDistTestRE;

	//////////////////////////////////////////////////////////////////////////
	//velocities to fit the force
	real *VxForce, *VyForce, *VzForce;
	//////////////////////////////////////////////////////////////////////////

	//Measure Points/////////
	std::vector<MeasurePoints> MP; 
	unsigned int* kMP;
	real* VxMP;
	real* VyMP;
	real* VzMP;
	real* RhoMP;
	unsigned int memSizerealkMP, memSizeIntkMP,    numberOfPointskMP;
	unsigned int numberOfValuesMP;

	//Drag Lift//////////////
	double *DragPreX, *DragPostX;
	double *DragPreY, *DragPostY;
	double *DragPreZ, *DragPostZ;
	std::vector<double> DragXvector;
	std::vector<double> DragYvector;
	std::vector<double> DragZvector;

	//2ndMoments////////////
	real *kxyFromfcNEQ, *kyzFromfcNEQ, *kxzFromfcNEQ, *kxxMyyFromfcNEQ, *kxxMzzFromfcNEQ;

	//3rdMoments////////////
	real *CUMbbb, *CUMabc, *CUMbac, *CUMbca, *CUMcba, *CUMacb, *CUMcab;
																				
	//HigherMoments/////////
	real *CUMcbb, *CUMbcb, *CUMbbc, *CUMcca, *CUMcac, *CUMacc, *CUMbcc, *CUMcbc, *CUMccb, *CUMccc;

	//CpTop/////////////////													
	int *cpTopIndex;															
	double *cpPressTop;															
	unsigned int numberOfPointsCpTop;							
	std::vector< std::vector< double > > cpTop;	
	std::vector< double > pressMirror;
	std::vector< bool > isOutsideInterface;
	unsigned int numberOfPointsPressWindow;

																
	//CpBottom/////////////										
	int *cpBottomIndex;											
	double *cpPressBottom;										
	unsigned int numberOfPointsCpBottom;						
	std::vector< std::vector< double > > cpBottom;

	//CpBottom2////////////
	int *cpBottom2Index;
	double *cpPressBottom2;
	unsigned int numberOfPointsCpBottom2;
	std::vector< std::vector< double > > cpBottom2;

	//Concentration////////
	int *concIndex;
	unsigned int numberOfPointsConc;

	//deltaPhi
	real deltaPhi;

	////////////////////////////////////////////////////////////////////////////
	//particles
	PathLineParticles plp;
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	//1D domain decomposition
	std::vector< ProcessNeighbor27 > sendProcessNeighbor;
	std::vector< ProcessNeighbor27 > recvProcessNeighbor;
	///////////////////////////////////////////////////////
	//3D domain decomposition
	std::vector< ProcessNeighbor27 > sendProcessNeighborX;
	std::vector< ProcessNeighbor27 > sendProcessNeighborY;
	std::vector< ProcessNeighbor27 > sendProcessNeighborZ;
	std::vector< ProcessNeighbor27 > recvProcessNeighborX;
	std::vector< ProcessNeighbor27 > recvProcessNeighborY;
	std::vector< ProcessNeighbor27 > recvProcessNeighborZ;
	///////////////////////////////////////////////////////
	//3D domain decomposition convection diffusion
	std::vector< ProcessNeighbor27 > sendProcessNeighborADX;
	std::vector< ProcessNeighbor27 > sendProcessNeighborADY;
	std::vector< ProcessNeighbor27 > sendProcessNeighborADZ;
	std::vector< ProcessNeighbor27 > recvProcessNeighborADX;
	std::vector< ProcessNeighbor27 > recvProcessNeighborADY;
	std::vector< ProcessNeighbor27 > recvProcessNeighborADZ;
	////////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////////
	////Restart
	//friend class boost::serialization::access;
	//template<class Archive>
	//void serialize(Archive & ar, const unsigned int version)
	//{
	// unsigned int i;
	// for (i=0; i<size_Mat_SP;i++)
	// {
	//  ar & d0SP.f[0][i];
	// }

	// for (i=0; i<size_Mat;i++)
	// {
	//  ar & k[i];
	// }
	//}
	////////////////////////////////////////////////////////////////////////////
};

class VF_PUBLIC Parameter
{
public:
	////////////////////////////////////////////////////////////////////////////
	////really ugly...should be in private...
	//Parameter();
	////////////////////////////////////////////////////////////////////////////
    static SPtr<Parameter> make();


	static Parameter* getInstanz();
	ParameterStruct* getParH(int level);
	ParameterStruct* getParD(int level);
	void initParameter();
	void fillSparse(int level);

	//measure points
	void copyMeasurePointsArrayToVector(int lev);

	//alloc
	void cudaAllocFull(int lev);
	void cudaFreeFull(int lev);

	void cudaAllocCoord(int lev);
	void cudaCopyCoord(int lev);
	void cudaFreeCoord(int lev);

	void cudaCopyPrint(int lev);
	void cudaCopyMedianPrint(int lev);

	void cudaAllocSP(int lev);
	void cudaCopySP(int lev);
	void cudaFreeSP(int lev);

	void cudaAllocF3SP(int lev);

	void cudaAllocNeighborWSB(int lev);
	void cudaCopyNeighborWSB(int lev);
	void cudaFreeNeighborWSB(int lev);

	void cudaAllocTurbulentViscosity(int lev);
	void cudaCopyTurbulentViscosityHD(int lev);
	void cudaCopyTurbulentViscosityDH(int lev);
	void cudaFreeTurbulentViscosity(int lev);

	void cudaAllocMedianSP(int lev);
	void cudaCopyMedianSP(int lev);
	void cudaFreeMedianSP(int lev);

	void cudaAllocMedianOut(int lev);
	void cudaFreeMedianOut(int lev);

	void cudaAllocInterfaceCF(int lev);
	void cudaCopyInterfaceCF(int lev);
	void cudaFreeInterfaceCF(int lev);
	void cudaAllocInterfaceFC(int lev);
	void cudaCopyInterfaceFC(int lev);
	void cudaFreeInterfaceFC(int lev);
	void cudaAllocInterfaceOffCF(int lev);
	void cudaCopyInterfaceOffCF(int lev);
	void cudaFreeInterfaceOffCF(int lev);
	void cudaAllocInterfaceOffFC(int lev);
	void cudaCopyInterfaceOffFC(int lev);
	void cudaFreeInterfaceOffFC(int lev);

	void cudaAllocVeloBC(int lev);
	void cudaCopyVeloBC(int lev);
	void cudaFreeVeloBC(int lev);
	void cudaAllocOutflowBC(int lev);
	void cudaCopyOutflowBC(int lev);
	void cudaFreeOutflowBC(int lev);
	void cudaAllocWallBC(int lev);
	void cudaCopyWallBC(int lev);
	void cudaFreeWallBC(int lev);
	void cudaAllocSlipBC(int lev);
	void cudaCopySlipBC(int lev);
	void cudaFreeSlipBC(int lev);

	void cudaAllocGeomValuesBC(int lev);
	void cudaCopyGeomValuesBC(int lev);
	void cudaFreeGeomValuesBC(int lev);
	void cudaAllocGeomBC(int lev);
	void cudaCopyGeomBC(int lev);
	void cudaFreeGeomBC(int lev);
	//Normals
	void cudaAllocGeomNormals(int lev);
	void cudaCopyGeomNormals(int lev);
	void cudaFreeGeomNormals(int lev);
	void cudaAllocInflowNormals(int lev);
	void cudaCopyInflowNormals(int lev);
	void cudaFreeInflowNormals(int lev);
	void cudaAllocOutflowNormals(int lev);
	void cudaCopyOutflowNormals(int lev);
	void cudaFreeOutflowNormals(int lev);

	void cudaAllocPress(int lev);
	void cudaCopyPress(int lev);
	void cudaFreePress(int lev);
	void cudaAllocTestRE(int lev, unsigned int size);
	void cudaCopyTestREtoDevice(int lev, unsigned int size);
	void cudaCopyTestREtoHost(int lev, unsigned int size);
	void cudaFreeTestRE(int lev);

	void cudaAllocCpTop(int lev);
	void cudaCopyCpTopInit(int lev);
	void cudaCopyCpTop(int lev);
	void cudaFreeCpTop(int lev);

	void cudaAllocCpBottom(int lev);
	void cudaCopyCpBottomInit(int lev);
	void cudaCopyCpBottom(int lev);
	void cudaFreeCpBottom(int lev);

	void cudaAllocCpBottom2(int lev);
	void cudaCopyCpBottom2Init(int lev);
	void cudaCopyCpBottom2(int lev);
	void cudaFreeCpBottom2(int lev);

	void cudaAllocConcFile(int lev);
	void cudaCopyConcFile(int lev);
	void cudaFreeConcFile(int lev);

	void cudaAllocInlet(int lev);
	void cudaCopyInlet(int lev);
	void cudaFreeInlet(int lev);
	void cudaAllocOutlet(int lev);
	void cudaCopyOutlet(int lev);
	void cudaFreeOutlet(int lev);


	void cudaAllocPressX0(int lev);
	void cudaCopyPressX0(int lev);
	void cudaFreePressX0(int lev);
	void cudaAllocPressX1(int lev);
	void cudaCopyPressX1(int lev);
	void cudaFreePressX1(int lev);

	void cudaAllocVeloPropeller(int lev);
	void cudaCopyVeloPropeller(int lev);
	void cudaFreeVeloPropeller(int lev);

	void cudaAllocMeasurePoints(int lev, int i);
	void cudaCopyMeasurePoints(int lev, int i);
	void cudaFreeMeasurePoints(int lev, int i);
	void cudaAllocMeasurePointsIndex(int lev);
	void cudaCopyMeasurePointsIndex(int lev);
	void cudaCopyMeasurePointsToHost(int lev);
	void cudaFreeMeasurePointsIndex(int lev);

	void cudaAllocFsForCheckPointAndRestart(int lev);
	void cudaCopyFsForRestart(int lev);
	void cudaCopyFsForCheckPoint(int lev);
	void cudaFreeFsForCheckPointAndRestart(int lev);

	void cudaAllocDragLift(int lev, int numofelem);
	void cudaCopyDragLift(int lev, int numofelem);
	void cudaFreeDragLift(int lev);

	void cudaAlloc2ndMoments(int lev, int numofelem);
	void cudaCopy2ndMoments(int lev, int numofelem);
	void cudaFree2ndMoments(int lev);

	void cudaAlloc3rdMoments(int lev, int numofelem);
	void cudaCopy3rdMoments(int lev, int numofelem);
	void cudaFree3rdMoments(int lev);

	void cudaAllocHigherMoments(int lev, int numofelem);
	void cudaCopyHigherMoments(int lev, int numofelem);
	void cudaFreeHigherMoments(int lev);

	void cudaAllocForceVelo(int lev, int numofelem);
	void cudaCopyForceVelo(int lev, int numofelem);
	void cudaFreeForceVelo(int lev);

	void cudaAllocForcing();
	void cudaCopyForcingToDevice();
	void cudaCopyForcingToHost();
	void cudaFreeForcing();

	//////////////////////////////////////////////////////////////////////////
	//Particles
	void cudaAllocParticles(int lev);
	void cudaCopyParticles(int lev);
	void cudaFreeParticles(int lev);
	//random value
	void cudaAllocRandomValues();
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//Porous Media
	void cudaAllocPorousMedia(PorousMedia* pm, int lev);
	void cudaCopyPorousMedia(PorousMedia* pm, int lev);
	void cudaFreePorousMedia(PorousMedia* pm, int lev);
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//Temperature
	void cudaAllocConc(int lev);
	void cudaCopyConcDH(int lev);
	void cudaCopyConcHD(int lev);
	void cudaFreeConc(int lev);
	//////////////////////////////////////////////////////////////////////////
	void cudaAllocTempFs(int lev);
	//////////////////////////////////////////////////////////////////////////
	void cudaAllocTempPressBC(int lev);
	void cudaCopyTempPressBCHD(int lev);
	void cudaFreeTempPressBC(int lev);
	//////////////////////////////////////////////////////////////////////////
	void cudaAllocTempVeloBC(int lev);
	void cudaCopyTempVeloBCHD(int lev);
	void cudaFreeTempVeloBC(int lev);
	//////////////////////////////////////////////////////////////////////////
	void cudaAllocTempNoSlipBC(int lev);
	void cudaCopyTempNoSlipBCHD(int lev);
	void cudaFreeTempNoSlipBC(int lev);
	//////////////////////////////////////////////////////////////////////////
	void cudaAllocPlaneConcIn(int lev, int numofelem);
	void cudaCopyPlaneConcIn(int lev, int numofelem);
	void cudaAllocPlaneConcOut1(int lev, int numofelem);
	void cudaCopyPlaneConcOut1(int lev, int numofelem);
	void cudaAllocPlaneConcOut2(int lev, int numofelem);
	void cudaCopyPlaneConcOut2(int lev, int numofelem);
	void cudaFreePlaneConc(int lev);
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//1D domain decomposition
	void cudaAllocProcessNeighbor(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighbor(int lev, unsigned int processNeighbor);
	//////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	void cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborXFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborXFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborX(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborY(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborYFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborYFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborYIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborY(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborZ(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborZFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborZFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborZIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborZ(int lev, unsigned int processNeighbor);
	//////////////////////////////////////////////////////////////////////////
	//3D domain decomposition convection diffusion
	void cudaAllocProcessNeighborADX(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADXFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADXFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADXIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborADX(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborADY(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADYFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADYFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADYIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborADY(int lev, unsigned int processNeighbor);
	//
	void cudaAllocProcessNeighborADZ(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADZFsHD(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADZFsDH(int lev, unsigned int processNeighbor);
	void cudaCopyProcessNeighborADZIndex(int lev, unsigned int processNeighbor);
	void cudaFreeProcessNeighborADZ(int lev, unsigned int processNeighbor);
	//////////////////////////////////////////////////////////////////////////





	//////////////////////////////////////////////////////////////////////////
	//setter
	void setForcing(real forcingX, real forcingY, real forcingZ);
	void setPhi(real inPhi);
	void setAngularVelocity(real inAngVel);
	void setStepEnsight(unsigned int step);
	void setOutputCount(unsigned int outputCount);
	void setlimitOfNodesForVTK(unsigned int limitOfNodesForVTK);
	void setStartTurn(unsigned int inStartTurn);
	void setSizeMatSparse(int level);
	void setDiffOn(bool isDiff);
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
	void setTEnd(unsigned int tend);
	void setTOut(unsigned int tout);
	void setTStartOut(unsigned int tStartOut);
	void setCalcMedian(bool calcMedian);
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
	void setFName(std::string fname);
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
	void setPrintFiles(bool printfiles);
	void setReadGeo(bool readGeo);
	void setTemperatureInit(real Temp);
	void setTemperatureBC(real TempBC);
	void setViscosity(real Viscosity);
	void setVelocity(real Velocity);
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
	void setUseMeasurePoints(bool useMeasurePoints);
	void setUseWale(bool useWale);
	void setSimulatePorousMedia(bool simulatePorousMedia);
	void setclockCycleForMP(real clockCycleForMP);
	void setDevices(std::vector<int> devices);
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
	void setNeedInterface(std::vector<bool> NeedInterface);
	void setTempH(TempforBoundaryConditions* TempH);
	void setTempD(TempforBoundaryConditions* TempD);
	void setTempVelH(TempVelforBoundaryConditions* TempVelH);
	void setTempVelD(TempVelforBoundaryConditions* TempVelD);
	void setTempPressH(TempPressforBoundaryConditions* TempPressH);
	void setTempPressD(TempPressforBoundaryConditions* TempPressD);
	void setTimeDoCheckPoint(unsigned int tDoCheckPoint);
	void setTimeDoRestart(unsigned int tDoRestart);   
	void setDoCheckPoint(bool doCheckPoint);
	void setDoRestart(bool doRestart);
	void setObj(std::string str, bool isObj);
	void setGeometryValues(bool GeometryValues);
	void setCalc2ndOrderMoments(bool is2ndOrderMoments);
	void setCalc3rdOrderMoments(bool is3rdOrderMoments);
	void setCalcHighOrderMoments(bool isHighOrderMoments);
	void setMemsizeGPU(double admem, bool reset);
	//1D domain decomposition
	void setPossNeighborFiles(std::vector<std::string> possNeighborFiles, std::string sor);
	void setNumberOfProcessNeighbors(unsigned int numberOfProcessNeighbors, int level, std::string sor);
	void setIsNeighbor(bool isNeighbor);
	//3D domain decomposition
	void setPossNeighborFilesX(std::vector<std::string> possNeighborFiles, std::string sor);
	void setPossNeighborFilesY(std::vector<std::string> possNeighborFiles, std::string sor);
	void setPossNeighborFilesZ(std::vector<std::string> possNeighborFiles, std::string sor);
	void setNumberOfProcessNeighborsX(unsigned int numberOfProcessNeighbors, int level, std::string sor);
	void setNumberOfProcessNeighborsY(unsigned int numberOfProcessNeighbors, int level, std::string sor);
	void setNumberOfProcessNeighborsZ(unsigned int numberOfProcessNeighbors, int level, std::string sor);
	void setIsNeighborX(bool isNeighbor);
	void setIsNeighborY(bool isNeighbor);
	void setIsNeighborZ(bool isNeighbor);
	//void setkInflowQ(unsigned int kInflowQ);
	//void setkOutflowQ(unsigned int kOutflowQ);
	//void setQinflowH(QforBoundaryConditions* QinflowH);
	//void setQinflowD(QforBoundaryConditions* QinflowD);
	//void setQoutflowH(QforBoundaryConditions* QoutflowH);
	//void setQoutflowD(QforBoundaryConditions* QoutflowD);
	//Normals
	void setgeomBoundaryNormalX(std::string geomNormalX);
	void setgeomBoundaryNormalY(std::string geomNormalY);
	void setgeomBoundaryNormalZ(std::string geomNormalZ);
	void setInflowBoundaryNormalX(std::string inflowNormalX);
	void setInflowBoundaryNormalY(std::string inflowNormalY);
	void setInflowBoundaryNormalZ(std::string inflowNormalZ);
	void setOutflowBoundaryNormalX(std::string outflowNormalX);
	void setOutflowBoundaryNormalY(std::string outflowNormalY);
	void setOutflowBoundaryNormalZ(std::string outflowNormalZ);

	//getter
	double* getForcesDouble();
	real* getForcesHost();
	real* getForcesDev();
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
	bool getPrintFiles();
	bool getReadGeo();
	bool getCalcMedian();
	bool getCalcParticle();
	int getFine();
	int getCoarse();
	int getParticleBasicLevel();
	int getParticleInitLevel();
	int getNumberOfParticles();
	int getDiffMod();
	int getFactorNZ();
	int getD3Qxx();
	int getMaxLevel();
	int getTimeCalcMedStart();
	int getTimeCalcMedEnd();
	int getMaxDev();
	int getMyID();
	int getNumprocs();
	std::string getOutputPath();
	std::string getOutputPrefix();
	std::string getFName();
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
	unsigned int getPressInID();
	unsigned int getPressOutID();
	unsigned int getPressInZ();
	unsigned int getPressOutZ();
	unsigned int getMemSizereal(int level);
	unsigned int getMemSizeInt(int level);
	unsigned int getMemSizeBool(int level);
	unsigned int getMemSizerealYZ(int level);
	unsigned int getSizeMat(int level);
	unsigned int getTStart();
	unsigned int getTInit();
	unsigned int getTEnd();
	unsigned int getTOut();
	unsigned int getTStartOut();
	unsigned int getTimestepForMP();
	real getDiffusivity();
	real getTemperatureInit();
	real getTemperatureBC();
	real getViscosity();
	real getVelocity();
	real getViscosityRatio();
	real getVelocityRatio();
	real getDensityRatio();
	real getPressRatio();
	real getRealX();
	real getRealY();
	real getRe();
	real getFactorPressBC();
	real getclockCycleForMP();
	std::vector<int> getDevices();
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
	std::vector<bool> getNeedInterface();
	TempforBoundaryConditions* getTempH();
	TempforBoundaryConditions* getTempD();
	TempVelforBoundaryConditions* getTempVelH();
	TempVelforBoundaryConditions* getTempVelD();
	TempPressforBoundaryConditions* getTempPressH();
	TempPressforBoundaryConditions* getTempPressD();
	unsigned int getTimeDoCheckPoint();
	unsigned int	getTimeDoRestart();   
	bool	getDoCheckPoint();
	bool	getDoRestart();
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
	bool getUseMeasurePoints();
	bool getUseWale();
	bool getSimulatePorousMedia();
	double getMemsizeGPU();
	//1D domain decomposition
	std::vector<std::string> getPossNeighborFiles(std::string sor);
	unsigned int getNumberOfProcessNeighbors(int level, std::string sor);
	bool getIsNeighbor();
	//3D domain decomposition
	std::vector<std::string> getPossNeighborFilesX(std::string sor);
	std::vector<std::string> getPossNeighborFilesY(std::string sor);
	std::vector<std::string> getPossNeighborFilesZ(std::string sor);
	unsigned int getNumberOfProcessNeighborsX(int level, std::string sor);
	unsigned int getNumberOfProcessNeighborsY(int level, std::string sor);
	unsigned int getNumberOfProcessNeighborsZ(int level, std::string sor);
	bool getIsNeighborX();
	bool getIsNeighborY();
	bool getIsNeighborZ();
	//Normals
	std::string getgeomBoundaryNormalX();
	std::string getgeomBoundaryNormalY();
	std::string getgeomBoundaryNormalZ();
	std::string getInflowBoundaryNormalX();
	std::string getInflowBoundaryNormalY();
	std::string getInflowBoundaryNormalZ();
	std::string getOutflowBoundaryNormalX();
	std::string getOutflowBoundaryNormalY();
	std::string getOutflowBoundaryNormalZ();
	//CUDA random number
	curandState* getRandomState();


    public:
        //Forcing///////////////
        real *forcingH, *forcingD;
        double hostForcing[3];

protected:
private:
	static Parameter* instanz;
	bool diffOn;
	int diffMod;
	int coarse, fine, maxlevel;
	int factor_gridNZ;
	int D3Qxx;
	InitCondition ic;
	double memsizeGPU;
	unsigned int limitOfNodesForVTK;
	unsigned int outputCount;

	//////////////////////////////////////////////////////////////////////////
	//particles
	int particleBasicLevel, particleInitLevel;
	int numberOfParticles;
	bool calcParticles;
	real stickToSolid;
	real startXHotWall, endXHotWall;
	//////////////////////////////////////////////////////////////////////////
	//CUDA random number generation
	curandState* devState;
	//////////////////////////////////////////////////////////////////////////

	//Temperature
	TempforBoundaryConditions *TempH, *TempD;
	//Temperature Velocity
	TempVelforBoundaryConditions *TempVelH, *TempVelD;
	//Temperature Pressure
	TempPressforBoundaryConditions *TempPressH, *TempPressD;

	//Drehung///////////////
	real Phi, angularVelocity;
	unsigned int startTurn;

	//Step of Ensight writing//
	unsigned int stepEnsight;

	std::vector<ParameterStruct*> parH;
	std::vector<ParameterStruct*> parD;
	//LogWriter output;

	Parameter();
	Parameter(const Parameter&);
	void initInterfaceParameter(int level);
	real TrafoXtoWorld(int CoordX, int level);
	real TrafoYtoWorld(int CoordY, int level);
	real TrafoZtoWorld(int CoordZ, int level);
public:
	real TrafoXtoMGsWorld(int CoordX, int level);
	real TrafoYtoMGsWorld(int CoordY, int level);
	real TrafoZtoMGsWorld(int CoordZ, int level);
private:
	//Multi GPGPU///////////////
	//1D domain decomposition
	std::vector<std::string> possNeighborFilesSend;
	std::vector<std::string> possNeighborFilesRecv;
	bool isNeigbor;
	//3D domain decomposition
	std::vector<std::string> possNeighborFilesSendX, possNeighborFilesSendY, possNeighborFilesSendZ;
	std::vector<std::string> possNeighborFilesRecvX, possNeighborFilesRecvY, possNeighborFilesRecvZ;
	bool isNeigborX, isNeigborY, isNeigborZ;
	////////////////////////////////////////////////////////////////////////////
	////Restart
	//friend class boost::serialization::access;
	//template<class Archive>
	//void serialize(Archive & ar, const unsigned int version)
	//{
	// unsigned int i;
	// int j;
	// for (j=coarse; j<=fine; j++)
	// {
	//  for (i=0; i<parH[j]->size_Mat_SP;i++)
	//  {
	//   ar & parH[j]->d0SP.f[0][i];
	//  }

	//  for (i=0; i<parH[j]->size_Mat;i++)
	//  {
	//   ar & parH[j]->k[j];
	//  }
	// }
	//}
	////////////////////////////////////////////////////////////////////////////
};

#endif


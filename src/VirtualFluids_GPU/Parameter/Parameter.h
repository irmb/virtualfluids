#ifndef PARAMETER_H
#define PARAMETER_H

#include <vector>
#include <string>
#include <memory>

#include "LBM/LB.h"
#include "LBM/D3Q27.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "VirtualFluids_GPU_EXPORT.h"

struct ParameterStruct
{
	bool evenOrOdd;
	unsigned int numberofthreads;

	//distributions///////////
	Distributions27 d0;
	Distributions27 d0SP;

	//geo/////////////////////
	int *geo;
	unsigned int *geoSP;

	//k///////////////////////
	unsigned int *k;

	//neighbor////////////////
	unsigned int *neighborX_SP, *neighborY_SP, *neighborZ_SP;

	//coordinates////////////
	doubflo *coordX_SP, *coordY_SP, *coordZ_SP;

	//vel parab///////////////
	doubflo *vParab;

	//macroscopic values//////
	doubflo *vx,    *vy,    *vz,    *rho;
	doubflo *vx_SP, *vy_SP, *vz_SP, *rho_SP, *press_SP;
	doubflo vis, omega;

	//grid////////////////////
	unsigned int nx, ny, nz;
	unsigned int gridNX, gridNY, gridNZ;

	//size of matrix//////////
	unsigned int size_Mat;
	unsigned int sizePlaneXY, sizePlaneYZ, sizePlaneXZ;

	//size of sparse matrix//////////
	unsigned int size_Mat_SP;

	//size of Plane btw. 2 GPUs//////
	unsigned int sizePlaneSB, sizePlaneRB, startB, endB;
	unsigned int sizePlaneST, sizePlaneRT, startT, endT;
	bool isSetSendB, isSetRecvB, isSetSendT, isSetRecvT;
	int *SendT, *SendB, *RecvT, *RecvB;

	//memsizeSP/////////////////
	unsigned int mem_size_doubflo_SP;
	unsigned int mem_size_int_SP;

	//memsize/////////////////
	unsigned int mem_size_doubflo;
	unsigned int mem_size_int;
	unsigned int mem_size_bool;
	unsigned int mem_size_doubflo_yz;

	//print///////////////////
	unsigned int startz, endz;
	doubflo Lx,Ly,Lz,dx;
	doubflo distX, distY, distZ;

	//BC's////////////////////
	QforBoundaryConditions  QWall,   Qinflow,      Qoutflow,      QPress;
	unsigned int            kQ,      kInflowQ,     kOutflowQ,     kPressQ;
	unsigned int            kQread,  kInflowQread, kOutflowQread, kPressQread;

	QforBoundaryConditions  QGeom;
	QforBoundaryConditions  QPeriodic;

	//deltaPhi
	doubflo deltaPhi;

	///////////////////////////////////////////////////////
	//3D domain decomposition
	std::vector< ProcessNeighbor27 > sendProcessNeighborX;
	std::vector< ProcessNeighbor27 > sendProcessNeighborY;
	std::vector< ProcessNeighbor27 > sendProcessNeighborZ;
	std::vector< ProcessNeighbor27 > recvProcessNeighborX;
	std::vector< ProcessNeighbor27 > recvProcessNeighborY;
	std::vector< ProcessNeighbor27 > recvProcessNeighborZ;
	////////////////////////////////////////////////////////////////////////////
};

class VirtualFluids_GPU_EXPORT Parameter
{
public:
	static std::shared_ptr<Parameter> make();
	ParameterStruct* getParH(int level);
	ParameterStruct* getParD(int level);
	void initParameter();

	//////////////////////////////////////////////////////////////////////////
	//setter
	void setForcing(doubflo forcingX, doubflo forcingY, doubflo forcingZ);
	void setlimitOfNodesForVTK(unsigned int limitOfNodesForVTK);
	void setD3Qxx(int d3qxx);
	void setMaxLevel(int maxlevel);
	void setTEnd(unsigned int tend);
	void setTOut(unsigned int tout);
	void setMaxDev(int maxdev);
	void setMyID(int myid);
	void setNumprocs(int numprocs);
	void setFName(std::string fname);
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
	void setgeomBoundaryBcQs(std::string geomBoundaryBcQs);
	void setgeomBoundaryBcValues(std::string geomBoundaryBcValues);
	void setnoSlipBcPos(std::string noSlipBcPos);
	void setnoSlipBcQs(std::string noSlipBcQs);
	void setnoSlipBcValue(std::string noSlipBcValue);
	void setnoSlipBcValues(std::string noSlipBcValues);
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
	void setnumberNodes(std::string numberNodes);
	void setLBMvsSI(std::string LBMvsSI);
	void setPrintFiles(bool printfiles);
	void setViscosity(doubflo Viscosity);
	void setVelocity(doubflo Velocity);
	void setViscosityRatio(doubflo ViscosityRatio);
	void setVelocityRatio(doubflo VelocityRatio);
	void setDensityRatio(doubflo DensityRatio);
	void setRe(doubflo Re);
	void setFactorPressBC(doubflo factorPressBC);
	void setIsGeo(bool isGeo);
	void setDevices(std::vector<int> devices);
	void setGridX(std::vector<int> GridX);
	void setGridY(std::vector<int> GridY);
	void setGridZ(std::vector<int> GridZ);
	void setDistX(std::vector<int> DistX);
	void setDistY(std::vector<int> DistY);
	void setDistZ(std::vector<int> DistZ);
	void setMinCoordX(std::vector<doubflo> MinCoordX);
	void setMinCoordY(std::vector<doubflo> MinCoordY);
	void setMinCoordZ(std::vector<doubflo> MinCoordZ);
	void setMaxCoordX(std::vector<doubflo> MaxCoordX);
	void setMaxCoordY(std::vector<doubflo> MaxCoordY);
	void setMaxCoordZ(std::vector<doubflo> MaxCoordZ);
	void setObj(std::string str, bool isObj);
	void setMemsizeGPU(double admem, bool reset);
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

	//getter
	double* getForcesDouble();
	doubflo* getForcesHost();
	doubflo* getForcesDev();
	unsigned int getlimitOfNodesForVTK();
	bool getEvenOrOdd(int level);
	bool getPrintFiles();
	int getFine();
	int getCoarse();
	int getFactorNZ();
	int getD3Qxx();
	int getMaxLevel();
	int getMaxDev();
	int getMyID();
	int getNumprocs();
	std::string getFName();
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
	std::string getgeomBoundaryBcQs();
	std::string getgeomBoundaryBcValues();
	std::string getnoSlipBcPos();
	std::string getnoSlipBcQs();
	std::string getnoSlipBcValue();
	std::string getnoSlipBcValues();
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
	std::string getnumberNodes();
	std::string getLBMvsSI();
	unsigned int getMemSizedoubflo(int level);
	unsigned int getMemSizeInt(int level);
	unsigned int getMemSizeBool(int level);
	unsigned int getMemSizedoubfloYZ(int level);
	unsigned int getSizeMat(int level);
	unsigned int getTStart();
	unsigned int getTInit();
	unsigned int getTEnd();
	unsigned int getTOut();
	doubflo getViscosity();
	doubflo getVelocity();
	doubflo getViscosityRatio();
	doubflo getVelocityRatio();
	doubflo getDensityRatio();
	doubflo getPressRatio();
	doubflo getRe();
	doubflo getFactorPressBC();
	std::vector<int> getDevices();
	std::vector<int> getGridX();
	std::vector<int> getGridY();
	std::vector<int> getGridZ();
	std::vector<int> getDistX();
	std::vector<int> getDistY();
	std::vector<int> getDistZ();
	std::vector<doubflo> getMinCoordX();
	std::vector<doubflo> getMinCoordY();
	std::vector<doubflo> getMinCoordZ();
	std::vector<doubflo> getMaxCoordX();
	std::vector<doubflo> getMaxCoordY();
	std::vector<doubflo> getMaxCoordZ();
	bool getIsGeo();
	double getMemsizeGPU();
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


public:
        //Forcing///////////////
        doubflo *forcingH, *forcingD;
        double hostForcing[3];

private:
	int coarse, fine, maxlevel;
	int factor_gridNZ;
	int D3Qxx;
	InitCondition ic;
	double memsizeGPU;
	unsigned int limitOfNodesForVTK;


	std::vector<ParameterStruct*> parH;
	std::vector<ParameterStruct*> parD;

	Parameter();
	Parameter(const Parameter&);

	//3D domain decomposition
	std::vector<std::string> possNeighborFilesSendX, possNeighborFilesSendY, possNeighborFilesSendZ;
	std::vector<std::string> possNeighborFilesRecvX, possNeighborFilesRecvY, possNeighborFilesRecvZ;
	bool isNeigborX, isNeigborY, isNeigborZ;
};

#endif


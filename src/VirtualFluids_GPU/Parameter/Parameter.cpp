#include "Parameter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

Parameter::Parameter()
{

}

std::shared_ptr<Parameter> Parameter::make()
{
    return std::shared_ptr<Parameter>(new Parameter());
}

void Parameter::initParameter()
{
	factor_gridNZ  = 2;
	coarse         = 0;
	fine           = this->maxlevel;
	parH.resize(this->maxlevel+1);
	parD.resize(this->maxlevel+1);

	//host
	for (int i = coarse; i <= fine; i++)
	{
		parH[i]                        = new ParameterStruct;
		parH[i]->numberofthreads       = 128;
		parH[i]->gridNX                = getGridX().at(i);
		parH[i]->gridNY                = getGridY().at(i);
		parH[i]->gridNZ                = getGridZ().at(i);
		parH[i]->vis                   = ic.vis*pow(2.f,i);
		parH[i]->omega                 = 1.0f/(3.0f*parH[i]->vis+0.5f);//omega :-) not s9 = -1.0f/(3.0f*parH[i]->vis+0.5f);//
		parH[i]->nx                    = parH[i]->gridNX + 2 * STARTOFFX;
		parH[i]->ny                    = parH[i]->gridNY + 2 * STARTOFFY;
		parH[i]->nz                    = parH[i]->gridNZ + 2 * STARTOFFZ;
		parH[i]->size_Mat              = parH[i]->nx * parH[i]->ny * parH[i]->nz;
		parH[i]->sizePlaneXY           = parH[i]->nx * parH[i]->ny;
		parH[i]->sizePlaneYZ           = parH[i]->ny * parH[i]->nz;
		parH[i]->sizePlaneXZ           = parH[i]->nx * parH[i]->nz;
		parH[i]->mem_size_doubflo      = sizeof(doubflo     ) * parH[i]->size_Mat;
		parH[i]->mem_size_int          = sizeof(unsigned int) * parH[i]->size_Mat;
		parH[i]->mem_size_bool         = sizeof(bool        ) * parH[i]->size_Mat;
		parH[i]->mem_size_doubflo_yz   = sizeof(doubflo     ) * parH[i]->ny * parH[i]->nz;
		parH[i]->evenOrOdd             = true;
		parH[i]->startz                = parH[i]->gridNZ * ic.myid;
		parH[i]->endz                  = parH[i]->gridNZ * ic.myid + parH[i]->gridNZ;
		parH[i]->Lx                    = (doubflo)((1.f*parH[i]->gridNX - 1.f)/(pow(2.f,i)));
		parH[i]->Ly                    = (doubflo)((1.f*parH[i]->gridNY - 1.f)/(pow(2.f,i)));
		parH[i]->Lz                    = (doubflo)((1.f*parH[i]->gridNZ - 1.f)/(pow(2.f,i)));
		parH[i]->dx                    = (doubflo)(1.f/(pow(2.f,i)));
		parH[i]->distX                 = (doubflo)getDistX().at(i);
		parH[i]->distY                 = (doubflo)getDistY().at(i);
		parH[i]->distZ                 = (doubflo)getDistZ().at(i);
	}

	//device
	for (int i = coarse; i <= fine; i++)
	{
		parD[i]                        = new ParameterStruct;
		parD[i]->numberofthreads       = parH[i]->numberofthreads;
		parD[i]->gridNX                = parH[i]->gridNX;
		parD[i]->gridNY                = parH[i]->gridNY;
		parD[i]->gridNZ                = parH[i]->gridNZ;
		parD[i]->vis                   = parH[i]->vis;
		parD[i]->omega                 = parH[i]->omega;
		parD[i]->nx                    = parH[i]->nx;
		parD[i]->ny                    = parH[i]->ny;
		parD[i]->nz                    = parH[i]->nz;
		parD[i]->size_Mat              = parH[i]->size_Mat;
		parD[i]->sizePlaneXY           = parH[i]->sizePlaneXY;
		parD[i]->sizePlaneYZ           = parH[i]->sizePlaneYZ;
		parD[i]->sizePlaneXZ           = parH[i]->sizePlaneXZ;
		parD[i]->mem_size_doubflo      = sizeof(doubflo     ) * parD[i]->size_Mat;
		parD[i]->mem_size_int          = sizeof(unsigned int) * parD[i]->size_Mat;
		parD[i]->mem_size_bool         = sizeof(bool        ) * parD[i]->size_Mat;
		parD[i]->mem_size_doubflo_yz   = sizeof(doubflo     ) * parD[i]->ny * parD[i]->nz;
		parD[i]->evenOrOdd             = parH[i]->evenOrOdd;
		parD[i]->startz                = parH[i]->startz;
		parD[i]->endz                  = parH[i]->endz;
		parD[i]->Lx                    = parH[i]->Lx;
		parD[i]->Ly                    = parH[i]->Ly;
		parD[i]->Lz                    = parH[i]->Lz;
		parD[i]->dx                    = parH[i]->dx;
		parD[i]->distX                 = parH[i]->distX;
		parD[i]->distY                 = parH[i]->distY;
		parD[i]->distZ                 = parH[i]->distZ;
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//set-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::setForcing(doubflo forcingX, doubflo forcingY, doubflo forcingZ)
{
	this->hostForcing[0] = forcingX;
	this->hostForcing[1] = forcingY;
	this->hostForcing[2] = forcingZ;
}
void Parameter::setlimitOfNodesForVTK(unsigned int limitOfNodesForVTK)
{
	this->limitOfNodesForVTK = limitOfNodesForVTK;
}
void Parameter::setD3Qxx(int d3qxx)
{
	this->D3Qxx = d3qxx;
}
void Parameter::setMaxLevel(int maxlevel)
{
	this->maxlevel = maxlevel-1;
}
void Parameter::setTEnd(unsigned int tend)
{
	ic.tend = tend;
}
void Parameter::setTOut(unsigned int tout)
{
	ic.tout = tout;
}
void Parameter::setFName(std::string fname)
{
	//std::string test = fname;
	ic.fname = fname;
}
void Parameter::setPrintFiles(bool printfiles)
{
	ic.printFiles = printfiles;
}
void Parameter::setViscosity(doubflo Viscosity)
{
	ic.vis = Viscosity;
}
void Parameter::setVelocity(doubflo Velocity)
{
	ic.u0 = Velocity;
}
void Parameter::setViscosityRatio(doubflo ViscosityRatio)
{
	ic.vis_ratio = ViscosityRatio;
}
void Parameter::setVelocityRatio(doubflo VelocityRatio)
{
	ic.u0_ratio = VelocityRatio;
}
void Parameter::setDensityRatio(doubflo DensityRatio)
{
	ic.delta_rho = DensityRatio;
}
void Parameter::setMaxDev(int maxdev)
{
	ic.maxdev = maxdev;
}
void Parameter::setMyID(int myid)
{
	ic.myid = myid;
}
void Parameter::setNumprocs(int numprocs)
{
	ic.numprocs = numprocs;
}
void Parameter::setDevices(std::vector<int> devices)
{
	ic.devices = devices;
}
void Parameter::setRe(doubflo Re)
{
	ic.Re = Re;
}
void Parameter::setFactorPressBC(doubflo factorPressBC)
{
	ic.factorPressBC = factorPressBC;
}
void Parameter::setIsGeo(bool isGeo)
{
	ic.isGeo = isGeo;
}
void Parameter::setGridX(std::vector<int> GridX)
{
	ic.GridX = GridX;
}
void Parameter::setGridY(std::vector<int> GridY)
{
	ic.GridY = GridY;
}
void Parameter::setGridZ(std::vector<int> GridZ)
{
	ic.GridZ = GridZ;
}
void Parameter::setDistX(std::vector<int> DistX)
{
	ic.DistX = DistX;
}
void Parameter::setDistY(std::vector<int> DistY)
{
	ic.DistY = DistY;
}
void Parameter::setDistZ(std::vector<int> DistZ)
{
	ic.DistZ = DistZ;
}
void Parameter::setMinCoordX(std::vector<doubflo> MinCoordX)
{
	ic.minCoordX = MinCoordX;
}
void Parameter::setMinCoordY(std::vector<doubflo> MinCoordY)
{
	ic.minCoordY = MinCoordY;
}
void Parameter::setMinCoordZ(std::vector<doubflo> MinCoordZ)
{
	ic.minCoordZ = MinCoordZ;
}
void Parameter::setMaxCoordX(std::vector<doubflo> MaxCoordX)
{
	ic.maxCoordX = MaxCoordX;
}
void Parameter::setMaxCoordY(std::vector<doubflo> MaxCoordY)
{
	ic.maxCoordY = MaxCoordY;
}
void Parameter::setMaxCoordZ(std::vector<doubflo> MaxCoordZ)
{
	ic.maxCoordZ = MaxCoordZ;
}
void Parameter::setkFull(std::string kFull)
{
	ic.kFull = kFull;
}
void Parameter::setgeoFull(std::string geoFull)
{
	ic.geoFull = geoFull;
}
void Parameter::setgeoVec(std::string geoVec)
{
	ic.geoVec = geoVec;
}
void Parameter::setcoordX(std::string coordX)
{
	ic.coordX = coordX;
}
void Parameter::setcoordY(std::string coordY)
{
	ic.coordY = coordY;
}
void Parameter::setcoordZ(std::string coordZ)
{
	ic.coordZ = coordZ;
}
void Parameter::setneighborX(std::string neighborX)
{
	ic.neighborX = neighborX;
}
void Parameter::setneighborY(std::string neighborY)
{
	ic.neighborY = neighborY;
}
void Parameter::setneighborZ(std::string neighborZ)
{
	ic.neighborZ = neighborZ;
}
void Parameter::setneighborWSB(std::string neighborWSB)
{
	ic.neighborWSB = neighborWSB;
}
void Parameter::setgeomBoundaryBcQs(std::string geomBoundaryBcQs)
{
	ic.geomBoundaryBcQs = geomBoundaryBcQs;
}
void Parameter::setgeomBoundaryBcValues(std::string geomBoundaryBcValues)
{
	ic.geomBoundaryBcValues = geomBoundaryBcValues;
}
void Parameter::setnoSlipBcPos(std::string noSlipBcPos)
{
	ic.noSlipBcPos = noSlipBcPos;
}
void Parameter::setnoSlipBcQs(std::string noSlipBcQs)
{
	ic.noSlipBcQs = noSlipBcQs;
}
void Parameter::setnoSlipBcValue(std::string noSlipBcValue)
{
	ic.noSlipBcValue = noSlipBcValue;
}
void Parameter::setnoSlipBcValues(std::string noSlipBcValues)
{
	ic.noSlipBcValues = noSlipBcValues;
}
void Parameter::setpressBcPos(std::string pressBcPos)
{
	ic.pressBcPos = pressBcPos;
}
void Parameter::setpressBcQs(std::string pressBcQs)
{
	ic.pressBcQs = pressBcQs;
}
void Parameter::setpressBcValue(std::string pressBcValue)
{
	ic.pressBcValue = pressBcValue;
}
void Parameter::setpressBcValues(std::string pressBcValues)
{
	ic.pressBcValues = pressBcValues;
}
void Parameter::setvelBcQs(std::string velBcQs)
{
	ic.velBcQs = velBcQs;
}
void Parameter::setvelBcValues(std::string velBcValues)
{
	ic.velBcValues = velBcValues;
}
void Parameter::setinletBcQs(std::string inletBcQs)
{
	ic.inletBcQs = inletBcQs;
}
void Parameter::setinletBcValues(std::string inletBcValues)
{
	ic.inletBcValues = inletBcValues;
}
void Parameter::setoutletBcQs(std::string outletBcQs)
{
	ic.outletBcQs = outletBcQs;
}
void Parameter::setoutletBcValues(std::string outletBcValues)
{
	ic.outletBcValues = outletBcValues;
}
void Parameter::settopBcQs(std::string topBcQs)
{
	ic.topBcQs = topBcQs;
}
void Parameter::settopBcValues(std::string topBcValues)
{
	ic.topBcValues = topBcValues;
}
void Parameter::setbottomBcQs(std::string bottomBcQs)
{
	ic.bottomBcQs = bottomBcQs;
}
void Parameter::setbottomBcValues(std::string bottomBcValues)
{
	ic.bottomBcValues = bottomBcValues;
}
void Parameter::setfrontBcQs(std::string frontBcQs)
{
	ic.frontBcQs = frontBcQs;
}
void Parameter::setfrontBcValues(std::string frontBcValues)
{
	ic.frontBcValues = frontBcValues;
}
void Parameter::setbackBcQs(std::string backBcQs)
{
	ic.backBcQs = backBcQs;
}
void Parameter::setbackBcValues(std::string backBcValues)
{
	ic.backBcValues = backBcValues;
}
void Parameter::setwallBcQs(std::string wallBcQs)
{
	ic.wallBcQs = wallBcQs;
}
void Parameter::setwallBcValues(std::string wallBcValues)
{
	ic.wallBcValues = wallBcValues;
}
void Parameter::setperiodicBcQs(std::string periodicBcQs)
{
	ic.periodicBcQs = periodicBcQs;
}
void Parameter::setperiodicBcValues(std::string periodicBcValues)
{
	ic.periodicBcValues = periodicBcValues;
}
void Parameter::setnumberNodes(std::string numberNodes)
{
	ic.numberNodes = numberNodes;
}
void Parameter::setLBMvsSI(std::string LBMvsSI)
{
	ic.LBMvsSI = LBMvsSI;
}
void Parameter::setObj(std::string str, bool isObj)
{
	if (str == "geo")
	{
		this->setIsGeo(isObj);
	}
}
void Parameter::setMemsizeGPU(double admem, bool reset)
{
	if (reset == true)
	{
		this->memsizeGPU = 0.;
	} 
	else
	{
		this->memsizeGPU += admem;
	}
}
//3D domain decomposition
void Parameter::setPossNeighborFilesX(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSendX = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecvX = possNeighborFiles;
	}
}
void Parameter::setPossNeighborFilesY(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSendY = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecvY = possNeighborFiles;
	}
}
void Parameter::setPossNeighborFilesZ(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSendZ = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecvZ = possNeighborFiles;
	}
}
void Parameter::setNumberOfProcessNeighborsX(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighborX.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighborX.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
	}
}
void Parameter::setNumberOfProcessNeighborsY(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighborY.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighborY.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
	}
}
void Parameter::setNumberOfProcessNeighborsZ(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighborZ.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighborZ.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double* Parameter::getForcesDouble()
{
	return this->hostForcing;
}
doubflo* Parameter::getForcesHost()
{
	return this->forcingH;
}
doubflo* Parameter::getForcesDev()
{
	return this->forcingD;
}
unsigned int Parameter::getlimitOfNodesForVTK()
{
	return this->limitOfNodesForVTK;
}
ParameterStruct* Parameter::getParD(int level)
{
	return parD[level];
}

ParameterStruct* Parameter::getParH(int level)
{
	return parH[level];
}
unsigned int Parameter::getSizeMat(int level)
{
	return parH[level]->size_Mat;
}
unsigned int Parameter::getMemSizedoubflo(int level)
{
	return parH[level]->mem_size_doubflo;
}
unsigned int Parameter::getMemSizeInt(int level)
{
	return parH[level]->mem_size_int;
}
unsigned int Parameter::getMemSizeBool(int level)
{
	return parH[level]->mem_size_bool;
}
unsigned int Parameter::getMemSizedoubfloYZ(int level)
{
	return parH[level]->mem_size_doubflo_yz;
}
int Parameter::getFine()
{
	return fine;
}
int Parameter::getCoarse()
{
	return coarse;
}
bool Parameter::getEvenOrOdd(int level)
{
	return parH[level]->evenOrOdd;
}
int Parameter::getFactorNZ()
{
	return factor_gridNZ;
}
int Parameter::getD3Qxx()
{
	return this->D3Qxx;
}
int Parameter::getMaxLevel()
{
	return this->maxlevel;
}
unsigned int Parameter::getTStart()
{
	return 1;
}
unsigned int Parameter::getTInit()
{
	return 0;
}
unsigned int Parameter::getTEnd()
{
	return ic.tend;
}
unsigned int Parameter::getTOut()
{
	return ic.tout;
}
std::string Parameter::getFName()
{
	return ic.fname;
}
bool Parameter::getPrintFiles()
{
	return ic.printFiles;
}
doubflo Parameter::getViscosity()
{
	return ic.vis;
}
doubflo Parameter::getVelocity()
{
	return ic.u0;
}
doubflo Parameter::getViscosityRatio()
{
	return ic.vis_ratio;
}
doubflo Parameter::getVelocityRatio()
{
	return ic.u0_ratio;
}
doubflo Parameter::getDensityRatio()
{
	return ic.delta_rho;
}
doubflo Parameter::getPressRatio()
{
	return ic.delta_press;
}
int Parameter::getMaxDev()
{
	return ic.maxdev;
}
int Parameter::getMyID()
{
	return ic.myid;
}
int Parameter::getNumprocs()
{
	return ic.numprocs;
}
std::vector<int> Parameter::getDevices()
{
	return ic.devices;
}
doubflo Parameter::getRe()
{
	return ic.Re;
}
doubflo Parameter::getFactorPressBC()
{
	return ic.factorPressBC;
}
std::vector<int> Parameter::getGridX()
{
	return ic.GridX;
}
std::vector<int> Parameter::getGridY()
{
	return ic.GridY;
}
std::vector<int> Parameter::getGridZ()
{
	return ic.GridZ;
}
std::vector<int> Parameter::getDistX()
{
	return ic.DistX;
}
std::vector<int> Parameter::getDistY()
{
	return ic.DistY;
}
std::vector<int> Parameter::getDistZ()
{
	return ic.DistZ;
}
std::vector<doubflo> Parameter::getMinCoordX()
{
	return ic.minCoordX;
}
std::vector<doubflo> Parameter::getMinCoordY()
{
	return ic.minCoordY;
}
std::vector<doubflo> Parameter::getMinCoordZ()
{
	return ic.minCoordZ;
}
std::vector<doubflo> Parameter::getMaxCoordX()
{
	return ic.maxCoordX;
}
std::vector<doubflo> Parameter::getMaxCoordY()
{
	return ic.maxCoordY;
}
std::vector<doubflo> Parameter::getMaxCoordZ()
{
	return ic.maxCoordZ;
}
std::string Parameter::getkFull()
{
	return ic.kFull;
}
std::string Parameter::getgeoFull()
{
	return ic.geoFull;
}
std::string Parameter::getgeoVec()
{
	return ic.geoVec;
}
std::string Parameter::getcoordX()
{
	return ic.coordX;
}
std::string Parameter::getcoordY()
{
	return ic.coordY;
}
std::string Parameter::getcoordZ()
{
	return ic.coordZ;
}
std::string Parameter::getneighborX()
{
	return ic.neighborX;
}
std::string Parameter::getneighborY()
{
	return ic.neighborY;
}
std::string Parameter::getneighborZ()
{
	return ic.neighborZ;
}
std::string Parameter::getneighborWSB()
{
	return ic.neighborWSB;
}
std::string Parameter::getgeomBoundaryBcQs()
{
	return ic.geomBoundaryBcQs;
}
std::string Parameter::getgeomBoundaryBcValues()
{
	return ic.geomBoundaryBcValues;
}
std::string Parameter::getnoSlipBcPos()
{
	return ic.noSlipBcPos;
}
std::string Parameter::getnoSlipBcQs()
{
	return ic.noSlipBcQs;
}
std::string Parameter::getnoSlipBcValue()
{
	return ic.noSlipBcValue;
}
std::string Parameter::getnoSlipBcValues()
{
	return ic.noSlipBcValues;
}
std::string Parameter::getpressBcPos()
{
	return ic.pressBcPos;
}
std::string Parameter::getpressBcQs()
{
	return ic.pressBcQs;
}
std::string Parameter::getpressBcValue()
{
	return ic.pressBcValue;
}
std::string Parameter::getpressBcValues()
{
	return ic.pressBcValues;
}
std::string Parameter::getvelBcQs()
{
	return ic.velBcQs;
}
std::string Parameter::getvelBcValues()
{
	return ic.velBcValues;
}
std::string Parameter::getinletBcQs()
{
	return ic.inletBcQs;
}
std::string Parameter::getinletBcValues()
{
	return ic.inletBcValues;
}
std::string Parameter::getoutletBcQs()
{
	return ic.outletBcQs;
}
std::string Parameter::getoutletBcValues()
{
	return ic.outletBcValues;
}
std::string Parameter::gettopBcQs()
{
	return ic.topBcQs;
}
std::string Parameter::gettopBcValues()
{
	return ic.topBcValues;
}
std::string Parameter::getbottomBcQs()
{
	return ic.bottomBcQs;
}
std::string Parameter::getbottomBcValues()
{
	return ic.bottomBcValues;
}
std::string Parameter::getfrontBcQs()
{
	return ic.frontBcQs;
}
std::string Parameter::getfrontBcValues()
{
	return ic.frontBcValues;
}
std::string Parameter::getbackBcQs()
{
	return ic.backBcQs;
}
std::string Parameter::getbackBcValues()
{
	return ic.backBcValues;
}
std::string Parameter::getwallBcQs()
{
	return ic.wallBcQs;
}
std::string Parameter::getwallBcValues()
{
	return ic.wallBcValues;
}
std::string Parameter::getperiodicBcQs()
{
	return ic.periodicBcQs;
}
std::string Parameter::getperiodicBcValues()
{
	return ic.periodicBcValues;
}
std::string Parameter::getLBMvsSI()
{
	return ic.LBMvsSI;
}
std::string Parameter::getnumberNodes()
{
	return ic.numberNodes;
}
bool Parameter::getIsGeo()
{
	return ic.isGeo;
}
double Parameter::getMemsizeGPU()
{
	return this->memsizeGPU;
}
//3D domain decomposition
std::vector<std::string> Parameter::getPossNeighborFilesX(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSendX;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecvX;
	}
}
std::vector<std::string> Parameter::getPossNeighborFilesY(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSendY;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecvY;
	}
}
std::vector<std::string> Parameter::getPossNeighborFilesZ(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSendZ;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecvZ;
	}
}
unsigned int Parameter::getNumberOfProcessNeighborsX(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighborX.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighborX.size();
	}
}
unsigned int Parameter::getNumberOfProcessNeighborsY(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighborY.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighborY.size();
	}
}
unsigned int Parameter::getNumberOfProcessNeighborsZ(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighborZ.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighborZ.size();
	}
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

#include "ConfigDataImp.h"


std::shared_ptr<ConfigDataImp> ConfigDataImp::getNewInstance()
{
	return std::shared_ptr<ConfigDataImp>(new ConfigDataImp());
}

ConfigDataImp::ConfigDataImp()
{
	this->isViscosity=false;
	this->isNumberOfDevices=false;
	this->isDevices=false;
	this->isOutputPath=false;
	this->isPrefix=false;
	this->isGridPath=false;
	this->isPrintOutputFiles=false;
	this->isGeometryValues=false;
	this->isCalc2ndOrderMoments=false;
	this->isCalc3rdOrderMoments=false;
	this->isCalcHighOrderMoments=false;
	this->isReadGeo=false;
	this->isCalcMedian=false;
	this->isConcFile=false;
	this->isUseMeasurePoints=false;
	this->isUseWale=false;
	this->isSimulatePorousMedia=false;
	this->isD3Qxx=false;
	this->isTEnd=false;
	this->isTOut=false;
	this->isTStartOut=false;
	this->isTimeCalcMedStart=false;
	this->isTimeCalcMedEnd=false;
	this->isPressInID=false;
	this->isPressOutID=false;
	this->isPressInZ=false;
	this->isPressOutZ=false;
	this->isDiffOn=false;
	this->isDiffMod=false;
	this->isDiffusivity=false;
	this->isTemperatureInit=false;
	this->isTemperatureBC=false;
	this->isVelocity=false;
	this->isViscosityRatio=false;
	this->isVelocityRatio=false;
	this->isDensityRatio=false;
	this->isPressRatio=false;
	this->isRealX=false;
	this->isRealY=false;
	this->isFactorPressBC=false;
	this->isGeometryFileC=false;
	this->isGeometryFileM=false;
	this->isGeometryFileF=false;
	this->isClockCycleForMP=false;
	this->isTimestepForMP=false;
	this->isForcingX=false;
	this->isForcingY=false;
	this->isForcingZ=false;
	this->isCalcParticles=false;
	this->isParticleBasicLevel=false;
	this->isParticleInitLevel=false;
	this->isNumberOfParticles=false;
	this->isNeighborWSB=false;
	this->isStartXHotWall=false;
	this->isEndXHotWall=false;
	this->isPossNeighborFilesX=false;
	this->isPossNeighborFilesY=false;
	this->isPossNeighborFilesZ=false;
	this->isTimeDoCheckPoint=false;
	this->isTimeDoRestart=false;
	this->isDoCheckPoint=false;
	this->isDoRestart=false;
	this->isMaxLevel=false;
	this->isGridX=false;
	this->isGridY=false;
	this->isGridZ=false;
	this->isDistX=false;
	this->isDistY=false;
	this->isDistZ=false;
	this->isNeedInterface=false;
	this->isMainKernel=false;
	this->isMultiKernelOn=false;
	this->isMultiKernelLevel=false;
	this->isMultiKernelName=false;

}

ConfigDataImp::~ConfigDataImp(void)
{
}

real ConfigDataImp::getViscosity()
{
	return this->viscosity;
}

uint ConfigDataImp::getNumberOfDevices()
{
	return this->numberOfDevices;
}

std::vector<uint> ConfigDataImp::getDevices()
{
	return this->devices;
}

std::string ConfigDataImp::getOutputPath()
{
	return this->outputPath;
}

std::string ConfigDataImp::getPrefix()
{
	return this->prefix;
}

std::string ConfigDataImp::getGridPath()
{
	return this->gridPath;
}

bool ConfigDataImp::getPrintOutputFiles()
{
	return this->printOutputFiles;
}

bool ConfigDataImp::getGeometryValues()
{
	return this->geometryValues;
}

bool ConfigDataImp::getCalc2ndOrderMoments()
{
	return this->calc2ndOrderMoments;
}

bool ConfigDataImp::getCalc3rdOrderMoments()
{
	return this->calc3rdOrderMoments;
}

bool ConfigDataImp::getCalcHighOrderMoments()
{
	return this->calcHighOrderMoments;
}

bool ConfigDataImp::getReadGeo()
{
	return this->readGeo;
}

bool ConfigDataImp::getCalcMedian()
{
	return this->calcMedian;
}

bool ConfigDataImp::getConcFile()
{
	return this->concFile;
}

bool ConfigDataImp::getUseMeasurePoints()
{
	return this->useMeasurePoints;
}

bool ConfigDataImp::getUseWale()
{
	return this->useWale;
}

bool ConfigDataImp::getSimulatePorousMedia()
{
	return this->simulatePorousMedia;
}

uint ConfigDataImp::getD3Qxx()
{
	return this->d3Qxx;
}

uint ConfigDataImp::getTEnd()
{
	return this->tEnd;
}

uint ConfigDataImp::getTOut()
{
	return this->tOut;
}

uint ConfigDataImp::getTStartOut()
{
	return this->tStartOut;
}

uint ConfigDataImp::getTimeCalcMedStart()
{
	return this->timeCalcMedStart;
}

uint ConfigDataImp::getTimeCalcMedEnd()
{
	return this->timeCalcMedEnd;
}

uint ConfigDataImp::getPressInID()
{
	return this->pressInID;
}

uint ConfigDataImp::getPressOutID()
{
	return this->pressOutID;
}

uint ConfigDataImp::getPressInZ()
{
	return this->pressInZ;
}

uint ConfigDataImp::getPressOutZ()
{
	return this->pressOutZ;
}

bool ConfigDataImp::getDiffOn()
{
	return this->diffOn;
}

uint ConfigDataImp::getDiffMod()
{
	return this->diffMod;
}

real ConfigDataImp::getDiffusivity()
{
	return this->diffusivity;
}

real ConfigDataImp::getTemperatureInit()
{
	return this->temperatureInit;
}

real ConfigDataImp::getTemperatureBC()
{
	return this->temperatureBC;
}

real ConfigDataImp::getVelocity()
{
	return this->velocity;
}

real ConfigDataImp::getViscosityRatio()
{
	return this->viscosityRatio;
}

real ConfigDataImp::getVelocityRatio()
{
	return this->velocityRatio;
}

real ConfigDataImp::getDensityRatio()
{
	return this->densityRatio;
}

real ConfigDataImp::getPressRatio()
{
	return this->pressRatio;
}

real ConfigDataImp::getRealX()
{
	return this->realX;
}

real ConfigDataImp::getRealY()
{
	return this->realY;
}

real ConfigDataImp::getFactorPressBC()
{
	return this->factorPressBC;
}

std::string ConfigDataImp::getGeometryFileC()
{
	return this->geometryFileC;
}

std::string ConfigDataImp::getGeometryFileM()
{
	return this->geometryFileM;
}

std::string ConfigDataImp::getGeometryFileF()
{
	return this->geometryFileF;
}

uint ConfigDataImp::getClockCycleForMP()
{
	return this->clockCycleForMP;
}

uint ConfigDataImp::getTimestepForMP()
{
	return this->timestepForMP;
}

real ConfigDataImp::getForcingX()
{
	return this->forcingX;
}

real ConfigDataImp::getForcingY()
{
	return this->forcingY;
}

real ConfigDataImp::getForcingZ()
{
	return this->forcingZ;
}

bool ConfigDataImp::getCalcParticles()
{
	return this->calcParticles;
}

int ConfigDataImp::getParticleBasicLevel()
{
	return this->particleBasicLevel;
}

int ConfigDataImp::getParticleInitLevel()
{
	return this->particleInitLevel;
}

int ConfigDataImp::getNumberOfParticles()
{
	return this->numberOfParticles;
}

real ConfigDataImp::getStartXHotWall()
{
	return this->startXHotWall;
}

real ConfigDataImp::getEndXHotWall()
{
	return this->endXHotWall;
}

std::vector<std::string> ConfigDataImp::getPossNeighborFilesX()
{
	return this->possNeighborFilesX;
}

std::vector<std::string> ConfigDataImp::getPossNeighborFilesY()
{
	return this->possNeighborFilesY;
}

std::vector<std::string> ConfigDataImp::getPossNeighborFilesZ()
{
	return this->possNeighborFilesZ;
}

int ConfigDataImp::getTimeDoCheckPoint()
{
	return this->timeDoCheckPoint;
}

int ConfigDataImp::getTimeDoRestart()
{
	return this->timeDoRestart;
}

bool ConfigDataImp::getDoCheckPoint()
{
	return this->doCheckPoint;
}

bool ConfigDataImp::getDoRestart()
{
	return this->doRestart;
}

uint ConfigDataImp::getMaxLevel()
{
	return this->maxLevel;
}

std::vector<int> ConfigDataImp::getGridX()
{
	return this->gridX;
}

std::vector<int> ConfigDataImp::getGridY()
{
	return this->gridY;
}

std::vector<int> ConfigDataImp::getGridZ()
{
	return this->gridZ;
}

std::vector<int> ConfigDataImp::getDistX()
{
	return this->distX;
}

std::vector<int> ConfigDataImp::getDistY()
{
	return this->distY;
}

std::vector<int> ConfigDataImp::getDistZ()
{
	return this->distZ;
}

std::vector<bool> ConfigDataImp::getNeedInterface()
{
	return this->needInterface;
}

std::string ConfigDataImp::getMainKernel()
{
	return this->mainKernel;
}

bool ConfigDataImp::getMultiKernelOn()
{
	return this->multiKernelOn;
}

std::vector<int> ConfigDataImp::getMultiKernelLevel()
{
	return this->multiKernelLevel;
}

std::vector<std::string> ConfigDataImp::getMultiKernelName()
{
	return this->multiKernelName;
}

void ConfigDataImp::setViscosity(real viscosity)
{
	this->viscosity = viscosity;
	this->isViscosity = true;
}

void ConfigDataImp::setNumberOfDevices(uint numberOfDevices)
{
	this->numberOfDevices = numberOfDevices;
	this->isNumberOfDevices = true;
}

void ConfigDataImp::setDevices(std::vector<uint> devices)
{
	this->devices = devices;
	this->isDevices = true;
}

void ConfigDataImp::setOutputPath(std::string outputPath)
{
	this->outputPath = outputPath;
	this->isOutputPath = true;
}

void ConfigDataImp::setPrefix(std::string prefix)
{
	this->prefix = prefix;
	this->isPrefix = true;
}

void ConfigDataImp::setGridPath(std::string gridPath)
{
	this->gridPath = gridPath;
	this->isGridPath = true;
}

void ConfigDataImp::setPrintOutputFiles(bool printOutputFiles)
{
	this->printOutputFiles = printOutputFiles;
	this->isPrintOutputFiles = true;
}

void ConfigDataImp::setGeometryValues(bool geometryValues)
{
	this->geometryValues = geometryValues;
	this->isGeometryValues = true;
}

void ConfigDataImp::setCalc2ndOrderMoments(bool calc2ndOrderMoments)
{
	this->calc2ndOrderMoments = calc2ndOrderMoments;
	this->isCalc2ndOrderMoments = true;
}

void ConfigDataImp::setCalc3rdOrderMoments(bool calc3rdOrderMoments)
{
	this->calc3rdOrderMoments = calc3rdOrderMoments;
	this->isCalc3rdOrderMoments = true;
}

void ConfigDataImp::setCalcHighOrderMoments(bool calcHighOrderMoments)
{
	this->calcHighOrderMoments = calcHighOrderMoments;
	this->isCalcHighOrderMoments = true;
}

void ConfigDataImp::setReadGeo(bool readGeo)
{
	this->readGeo = readGeo;
	this->isReadGeo = true;
}

void ConfigDataImp::setCalcMedian(bool calcMedian)
{
	this->calcMedian = calcMedian;
	this->isCalcMedian = true;
}

void ConfigDataImp::setConcFile(bool concFile)
{
	this->concFile = concFile;
	this->isConcFile = true;
}

void ConfigDataImp::setUseMeasurePoints(bool useMeasurePoints)
{
	this->useMeasurePoints = useMeasurePoints;
	this->isUseMeasurePoints = true;
}

void ConfigDataImp::setUseWale(bool useWale)
{
	this->useWale = useWale;
	this->isUseWale = true;
}

void ConfigDataImp::setSimulatePorousMedia(bool simulatePorousMedia)
{
	this->simulatePorousMedia = simulatePorousMedia;
	this->isSimulatePorousMedia = true;
}

void ConfigDataImp::setD3Qxx(uint d3Qxx)
{
	this->d3Qxx = d3Qxx;
	this->isD3Qxx = true;
}

void ConfigDataImp::setTEnd(uint tEnd)
{
	this->tEnd = tEnd;
	this->isTEnd = true;
}

void ConfigDataImp::setTOut(uint tOut)
{
	this->tOut = tOut;
	this->isTOut = true;
}

void ConfigDataImp::setTStartOut(uint tStartOut)
{
	this->tStartOut = tStartOut;
	this->isTStartOut = true;
}

void ConfigDataImp::setTimeCalcMedStart(uint timeCalcMedStart)
{
	this->timeCalcMedStart = timeCalcMedStart;
	this->isTimeCalcMedStart = true;
}

void ConfigDataImp::setTimeCalcMedEnd(uint timeCalcMedEnd)
{
	this->timeCalcMedEnd = timeCalcMedEnd;
	this->isTimeCalcMedEnd = true;
}

void ConfigDataImp::setPressInID(uint pressInID)
{
	this->pressInID = pressInID;
	this->isPressInID = true;
}

void ConfigDataImp::setPressOutID(uint pressOutID)
{
	this->pressOutID = pressOutID;
	this->isPressOutID = true;
}

void ConfigDataImp::setPressInZ(uint pressInZ)
{
	this->pressInZ = pressInZ;
	this->isPressInZ = true;
}

void ConfigDataImp::setPressOutZ(uint pressOutZ)
{
	this->pressOutZ = pressOutZ;
	this->isPressOutZ = true;
}

void ConfigDataImp::setDiffOn(bool diffOn)
{
	this->diffOn = diffOn;
	this->isDiffOn = true;
}

void ConfigDataImp::setDiffMod(uint diffMod)
{
	this->diffMod = diffMod;
	this->isDiffMod = true;
}

void ConfigDataImp::setDiffusivity(real diffusivity)
{
	this->diffusivity = diffusivity;
	this->isDiffusivity = true;
}

void ConfigDataImp::setTemperatureInit(real temperatureInit)
{
	this->temperatureInit = temperatureInit;
	this->isTemperatureInit = true;
}

void ConfigDataImp::setTemperatureBC(real temperatureBC)
{
	this->temperatureBC = temperatureBC;
	this->isTemperatureBC = true;
}

//void ConfigDataImp::setViscosity(real viscosity)
//{
//	this->viscosity = viscosity;
//	this->isViscosity = true;
//}

void ConfigDataImp::setVelocity(real velocity)
{
	this->velocity = velocity;
	this->isVelocity = true;
}

void ConfigDataImp::setViscosityRatio(real viscosityRatio)
{
	this->viscosityRatio = viscosityRatio;
	this->isViscosityRatio = true;
}

void ConfigDataImp::setVelocityRatio(real velocityRatio)
{
	this->velocityRatio = velocityRatio;
	this->isVelocityRatio = true;
}

void ConfigDataImp::setDensityRatio(real densityRatio)
{
	this->densityRatio = densityRatio;
	this->isDensityRatio = true;
}

void ConfigDataImp::setPressRatio(real pressRatio)
{
	this->pressRatio = pressRatio;
	this->isPressRatio = true;
}

void ConfigDataImp::setRealX(real realX)
{
	this->realX = realX;
	this->isRealX = true;
}

void ConfigDataImp::setRealY(real realY)
{
	this->realY = realY;
	this->isRealY = true;
}

void ConfigDataImp::setFactorPressBC(real factorPressBC)
{
	this->factorPressBC = factorPressBC;
	this->isFactorPressBC = true;
}

void ConfigDataImp::setGeometryFileC(std::string geometryFileC)
{
	this->geometryFileC = geometryFileC;
	this->isGeometryFileC = true;
}

void ConfigDataImp::setGeometryFileM(std::string geometryFileM)
{
	this->geometryFileM = geometryFileM;
	this->isGeometryFileM = true;
}

void ConfigDataImp::setGeometryFileF(std::string geometryFileF)
{
	this->geometryFileF = geometryFileF;
	this->isGeometryFileF = true;
}

void ConfigDataImp::setClockCycleForMP(uint clockCycleForMP)
{
	this->clockCycleForMP = clockCycleForMP;
	this->isClockCycleForMP = true;
}

void ConfigDataImp::setTimestepForMP(uint timestepForMP)
{
	this->timestepForMP = timestepForMP;
	this->isTimestepForMP = true;
}

void ConfigDataImp::setForcingX(real forcingX)
{
	this->forcingX = forcingX;
	this->isForcingX = true;
}

void ConfigDataImp::setForcingY(real forcingY)
{
	this->forcingY = forcingY;
	this->isForcingY = true;
}

void ConfigDataImp::setForcingZ(real forcingZ)
{
	this->forcingZ = forcingZ;
	this->isForcingZ = true;
}

void ConfigDataImp::setCalcParticles(bool calcParticles)
{
	this->calcParticles = calcParticles;
	this->isCalcParticles = true;
}

void ConfigDataImp::setParticleBasicLevel(int particleBasicLevel)
{
	this->particleBasicLevel = particleBasicLevel;
	this->isParticleBasicLevel = true;
}

void ConfigDataImp::setParticleInitLevel(int particleInitLevel)
{
	this->particleInitLevel = particleInitLevel;
	this->isParticleInitLevel = true;
}

void ConfigDataImp::setNumberOfParticles(int numberOfParticles)
{
	this->numberOfParticles = numberOfParticles;
	this->isNumberOfParticles = true;
}

void ConfigDataImp::setStartXHotWall(real startXHotWall)
{
	this->startXHotWall = startXHotWall;
	this->isStartXHotWall = true;
}

void ConfigDataImp::setEndXHotWall(real endXHotWall)
{
	this->endXHotWall = endXHotWall;
	this->isEndXHotWall = true;
}

void ConfigDataImp::setPossNeighborFilesX(std::vector<std::string> possNeighborFilesX)
{
	this->possNeighborFilesX = possNeighborFilesX;
	this->isPossNeighborFilesX = true;
}

void ConfigDataImp::setPossNeighborFilesY(std::vector<std::string> possNeighborFilesY)
{
	this->possNeighborFilesY = possNeighborFilesY;
	this->isPossNeighborFilesY = true;
}

void ConfigDataImp::setPossNeighborFilesZ(std::vector<std::string> possNeighborFilesZ)
{
	this->possNeighborFilesZ = possNeighborFilesZ;
	this->isPossNeighborFilesZ = true;
}

void ConfigDataImp::setTimeDoCheckPoint(int timeDoCheckPoint)
{
	this->timeDoCheckPoint = timeDoCheckPoint;
	this->isTimeDoCheckPoint = true;
}

void ConfigDataImp::setTimeDoRestart(int timeDoRestart)
{
	this->timeDoRestart = timeDoRestart;
	this->isTimeDoRestart = true;
}

void ConfigDataImp::setDoCheckPoint(bool doCheckPoint)
{
	this->doCheckPoint = doCheckPoint;
	this->isDoCheckPoint = true;
}

void ConfigDataImp::setDoRestart(bool doRestart)
{
	this->doRestart = doRestart;
	this->isDoRestart = true;
}

void ConfigDataImp::setMaxLevel(uint maxLevel)
{
	this->maxLevel = maxLevel;
	this->isMaxLevel = true;
}

void ConfigDataImp::setGridX(std::vector<int> gridX)
{
	this->gridX = gridX;
	this->isGridX = true;
}

void ConfigDataImp::setGridY(std::vector<int> gridY)
{
	this->gridY = gridY;
	this->isGridY = true;
}

void ConfigDataImp::setGridZ(std::vector<int> gridZ)
{
	this->gridZ = gridZ;
	this->isGridZ = true;
}

void ConfigDataImp::setDistX(std::vector<int> distX)
{
	this->distX = distX;
	this->isDistX = true;
}

void ConfigDataImp::setDistY(std::vector<int> distY)
{
	this->distY = distY;
	this->isDistY = true;
}

void ConfigDataImp::setDistZ(std::vector<int> distZ)
{
	this->distZ = distZ;
	this->isDistZ = true;
}

void ConfigDataImp::setNeedInterface(std::vector<bool> needInterface)
{
	this->needInterface = needInterface;
	this->isNeedInterface = true;
}

void ConfigDataImp::setMainKernel(std::string mainKernel)
{
	this->mainKernel = mainKernel;
	this->isMainKernel = true;
}

void ConfigDataImp::setMultiKernelOn(bool multiKernelOn)
{
	this->multiKernelOn = multiKernelOn;
	this->isMultiKernelOn = true;
}

void ConfigDataImp::setMultiKernelLevel(std::vector<int> multiKernelLevel)
{
	this->multiKernelLevel = multiKernelLevel;
	this->isMultiKernelLevel = true;
}

void ConfigDataImp::setMultiKernelName(std::vector<std::string> multiKernelName)
{
	this->multiKernelName = multiKernelName;
	this->isMultiKernelName = true;
}

bool ConfigDataImp::isCalc2ndOrderMomentsInConfigFile()
{
	return this->isCalc2ndOrderMoments;
}

bool ConfigDataImp::isCalc3rdOrderMomentsInConfigFile()
{
	return this->isCalc2ndOrderMoments;
}

bool ConfigDataImp::isCalcHighOrderMomentsInConfigFile()
{
	return this->isCalcHighOrderMoments;
}

bool ConfigDataImp::isReadGeoInConfigFile()
{
	return this->isReadGeo;
}

bool ConfigDataImp::isCalcMedianInConfigFile()
{
	return this->isCalcMedian;
}

bool ConfigDataImp::isConcFileInConfigFile()
{
	return this->isConcFile;
}

bool ConfigDataImp::isUseMeasurePointsInConfigFile()
{
	return this->isUseMeasurePoints;
}

bool ConfigDataImp::isUseWaleInConfigFile()
{
	return this->isUseWale;
}

bool ConfigDataImp::isSimulatePorousMediaInConfigFile()
{
	return this->isSimulatePorousMedia;
}

bool ConfigDataImp::isD3QxxInConfigFile()
{
	return this->isD3Qxx;
}

bool ConfigDataImp::isTEndInConfigFile()
{
	return this->isTEnd;
}

bool ConfigDataImp::isTOutInConfigFile()
{
	return this->isTOut;
}

bool ConfigDataImp::isTStartOutInConfigFile()
{
	return this->isTStartOut;
}

bool ConfigDataImp::isTimeCalcMedStartInConfigFile()
{
	return this->isTimeCalcMedStart;
}

bool ConfigDataImp::isTimeCalcMedEndInConfigFile()
{
	return this->isTimeCalcMedEnd;
}

bool ConfigDataImp::isPressInIDInConfigFile()
{
	return this->isPressInID;
}

bool ConfigDataImp::isPressOutIDInConfigFile()
{
	return this->isPressOutID;
}

bool ConfigDataImp::isPressInZInConfigFile()
{
	return this->isPressInZ;
}

bool ConfigDataImp::isPressOutZInConfigFile()
{
	return this->isPressOutZ;
}

bool ConfigDataImp::isDiffOnInConfigFile()
{
	return this->isDiffOn;
}

bool ConfigDataImp::isDiffModInConfigFile()
{
	return this->isDiffMod;
}

bool ConfigDataImp::isDiffusivityInConfigFile()
{
	return this->isDiffusivity;
}

bool ConfigDataImp::isTemperatureInitInConfigFile()
{
	return this->isTemperatureInit;
}

bool ConfigDataImp::isTemperatureBCInConfigFile()
{
	return this->isTemperatureBC;
}

bool ConfigDataImp::isViscosityInConfigFile()
{
	return this->isViscosity;
}

bool ConfigDataImp::isNumberOfDevicesInConfigFile()
{
	return this->isNumberOfDevices;
}

bool ConfigDataImp::isDevicesInConfigFile()
{
	return this->isDevices;
}

bool ConfigDataImp::isOutputPathInConfigFile()
{
	return this->isOutputPath;
}

bool ConfigDataImp::isPrefixInConfigFile()
{
	return this->isPrefix;
}

bool ConfigDataImp::isGridPathInConfigFile()
{
	return this->isGridPath;
}

bool ConfigDataImp::isPrintOutputFilesInConfigFile()
{
	return this->isPrintOutputFiles;
}

bool ConfigDataImp::isGeometryValuesInConfigFile()
{
	return this->isGeometryValues;
}

bool ConfigDataImp::isVelocityInConfigFile()
{
	return this->isVelocity;
}

bool ConfigDataImp::isViscosityRatioInConfigFile()
{
	return this->isViscosityRatio;
}

bool ConfigDataImp::isVelocityRatioInConfigFile()
{
	return this->isVelocityRatio;
}

bool ConfigDataImp::isDensityRatioInConfigFile()
{
	return this->isDensityRatio;
}

bool ConfigDataImp::isPressRatioInConfigFile()
{
	return this->isPressRatio;
}

bool ConfigDataImp::isRealXInConfigFile()
{
	return this->isRealX;
}

bool ConfigDataImp::isRealYInConfigFile()
{
	return this->isRealY;
}

bool ConfigDataImp::isFactorPressBCInConfigFile()
{
	return this->isFactorPressBC;
}

bool ConfigDataImp::isGeometryFileCInConfigFile()
{
	return this->isGeometryFileC;
}

bool ConfigDataImp::isGeometryFileMInConfigFile()
{
	return this->isGeometryFileM;
}

bool ConfigDataImp::isGeometryFileFInConfigFile()
{
	return this->isGeometryFileF;
}

bool ConfigDataImp::isClockCycleForMPInConfigFile()
{
	return this->isClockCycleForMP;
}

bool ConfigDataImp::isTimestepForMPInConfigFile()
{
	return this->isTimestepForMP;
}

bool ConfigDataImp::isForcingXInConfigFile()
{
	return this->isForcingX;
}

bool ConfigDataImp::isForcingYInConfigFile()
{
	return this->isForcingY;
}

bool ConfigDataImp::isForcingZInConfigFile()
{
	return this->isForcingZ;
}

bool ConfigDataImp::isCalcParticlesInConfigFile()
{
	return this->isCalcParticles;
}

bool ConfigDataImp::isParticleBasicLevelInConfigFile()
{
	return this->isParticleBasicLevel;
}

bool ConfigDataImp::isParticleInitLevelInConfigFile()
{
	return this->isParticleInitLevel;
}

bool ConfigDataImp::isNumberOfParticlesInConfigFile()
{
	return this->isNumberOfParticles;
}

bool ConfigDataImp::isNeighborWSBInConfigFile()
{
	return this->isNeighborWSB;
}

bool ConfigDataImp::isStartXHotWallInConfigFile()
{
	return this->isStartXHotWall;
}

bool ConfigDataImp::isEndXHotWallInConfigFile()
{
	return this->isEndXHotWall;
}

bool ConfigDataImp::isPossNeighborFilesXInConfigFile()
{
	return this->isPossNeighborFilesX;
}

bool ConfigDataImp::isPossNeighborFilesYInConfigFile()
{
	return this->isPossNeighborFilesY;
}

bool ConfigDataImp::isPossNeighborFilesZInConfigFile()
{
	return this->isPossNeighborFilesZ;
}

bool ConfigDataImp::isTimeDoCheckPointInConfigFile()
{
	return this->isTimeDoCheckPoint;
}

bool ConfigDataImp::isTimeDoRestartInConfigFile()
{
	return this->isTimeDoCheckPoint;
}

bool ConfigDataImp::isDoCheckPointInConfigFile()
{
	return this->isDoCheckPoint;
}

bool ConfigDataImp::isDoRestartInConfigFile()
{
	return this->isDoRestart;
}

bool ConfigDataImp::isMaxLevelInConfigFile()
{
	return this->isMaxLevel;
}

bool ConfigDataImp::isGridXInConfigFile()
{
	return this->isGridX;
}

bool ConfigDataImp::isGridYInConfigFile()
{
	return this->isGridY;
}

bool ConfigDataImp::isGridZInConfigFile()
{
	return this->isGridZ;
}

bool ConfigDataImp::isDistXInConfigFile()
{
	return this->isDistX;
}

bool ConfigDataImp::isDistYInConfigFile()
{
	return this->isDistY;
}

bool ConfigDataImp::isDistZInConfigFile()
{
	return this->isDistZ;
}

bool ConfigDataImp::isNeedInterfaceInConfigFile()
{
	return this->isNeedInterface;
}

bool ConfigDataImp::isMainKernelInConfigFile()
{
	return this->isMainKernel;
}

bool ConfigDataImp::isMultiKernelOnInConfigFile()
{
	return this->isMultiKernelOn;
}

bool ConfigDataImp::isMultiKernelLevelInConfigFile()
{
	return this->isMultiKernelLevel;
}

bool ConfigDataImp::isMultiKernelNameInConfigFile()
{
	return this->isMultiKernelName;
}



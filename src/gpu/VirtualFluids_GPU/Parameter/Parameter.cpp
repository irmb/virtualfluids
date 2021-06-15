//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include "Parameter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <curand_kernel.h>


#include "Core/Input/ConfigData/ConfigData.h"
#include "Core/StringUtilities/StringUtil.h"
#include "Communication/Communicator.h"
//#ifdef WIN32
//   #include <Winsock2.h>
//#endif
//lib for windows Ws2_32.lib

#include <basics/config/ConfigurationFile.h>


SPtr<Parameter> Parameter::make(SPtr<ConfigData> configData, vf::gpu::Communicator* comm)
{
	return SPtr<Parameter>(new Parameter(configData, comm));
}


Parameter::Parameter(const vf::gpu::Communicator& comm)
{
    ic.numprocs = comm.getNummberOfProcess();
    ic.myid = comm.getPID();
}


Parameter::Parameter(const vf::basics::ConfigurationFile& configData,
                     const vf::gpu::Communicator& comm) :
                     Parameter(comm)

{
    if (configData.contains("NumberOfDevices"))
        this->setMaxDev(configData.getValue<int>("NumberOfDevices"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("Devices"))
        this->setDevices(configData.getVector<uint>("Devices"));
    //////////////////////////////////////////////////////////////////////////
	if (configData.contains("Path"))
		this->setOutputPath(configData.getValue<std::string>("Path"));
	else
		this->setOutputPath("C:/Output/"); //TODO: Shouldnt we throw an exception here?
    //////////////////////////////////////////////////////////////////////////
	if (configData.contains("Prefix"))
		this->setOutputPrefix(configData.getValue<std::string>("Prefix"));
    //////////////////////////////////////////////////////////////////////////
    if (configData.contains("WriteGrid"))
		this->setPrintFiles(configData.getValue<bool>("WriteGrid"));
    //////////////////////////////////////////////////////////////////////////
	if (configData.contains("GeometryValues"))
		this->setGeometryValues(configData.getValue<bool>("GeometryValues"));
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
	if (configData.contains("calcMedian"))
		this->setCalcMedian(configData.getValue<bool>("calcMedian"));
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
	if (configData.contains("UseStreetVelocityFile"))
		this->setStreetVelocityFile(configData.getValue<bool>("UseStreetVelocityFile"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("UseMeasurePoints"))
		this->setUseMeasurePoints(configData.getValue<bool>("UseMeasurePoints"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("UseWale"))
		this->setUseWale(configData.getValue<bool>("UseWale"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("UseInitNeq"))
		this->setUseInitNeq(configData.getValue<bool>("UseInitNeq"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("SimulatePorousMedia"))
		this->setSimulatePorousMedia(configData.getValue<bool>("SimulatePorousMedia"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("D3Qxx"))
		this->setD3Qxx(configData.getValue<int>("D3Qxx"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("TimeEnd"))
		this->setTEnd(configData.getValue<int>("TimeEnd"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("TimeOut"))
		this->setTOut(configData.getValue<int>("TimeOut"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("TimeStartOut"))
		this->setTStartOut(configData.getValue<int>("TimeStartOut"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("TimeStartCalcMedian"))
		this->setTimeCalcMedStart(configData.getValue<int>("TimeStartCalcMedian"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("TimeEndCalcMedian"))
		this->setTimeCalcMedEnd(configData.getValue<int>("TimeEndCalcMedian"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("PressInID"))
		this->setTOut(configData.getValue<int>("PressInID"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("PressOutID"))
		this->setTStartOut(configData.getValue<int>("PressOutID"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("PressInZ"))
		this->setTimeCalcMedStart(configData.getValue<int>("PressInZ"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("PressOutZ"))
		this->setTimeCalcMedEnd(configData.getValue<int>("PressOutZ"));

	// //////////////////////////////////////////////////////////////////////////
	// //second component

	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("DiffOn"))
		this->setDiffOn(configData.getValue<bool>("DiffOn"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("DiffMod"))
		this->setDiffMod(configData.getValue<int>("DiffMod"));
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
		this->setViscosity(configData.getValue<real>("Viscosity_LB"));
	//////////////////////////////////////////////////////////////////////////
	if (configData.contains("Velocity_LB"))
		this->setVelocity(configData.getValue<real>("Velocity_LB"));
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
	// //////////////////////////////////////////////////////////////////////////
	// //read Geometry (STL)
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

	std::string gridPath {""};
	if (configData.contains("GridPath"))
		gridPath = configData.getValue<std::string>("GridPath");
	else
		throw std::runtime_error("GridPath has to be defined in config file!");

	if (this->getNumprocs() == 1)
		gridPath += "/";
	else
		gridPath += "/" + StringUtil::toString(this->getMyID()) + "/";
	
	// //////////////////////////////////////////////////////////////////////////
	this->setFName(this->getOutputPath() + "/" + this->getOutputPrefix());
	//////////////////////////////////////////////////////////////////////////
	this->setgeoVec(				gridPath + "geoVec.dat");
	this->setcoordX(				gridPath + "coordX.dat");
	this->setcoordY(				gridPath + "coordY.dat");
	this->setcoordZ(				gridPath + "coordZ.dat");
	this->setneighborX(				gridPath + "neighborX.dat");
	this->setneighborY(				gridPath + "neighborY.dat");
	this->setneighborZ(				gridPath + "neighborZ.dat");
	this->setneighborWSB(			gridPath + "neighborWSB.dat");
	this->setscaleCFC(				gridPath + "scaleCFC.dat");
	this->setscaleCFF(				gridPath + "scaleCFF.dat");
	this->setscaleFCC(				gridPath + "scaleFCC.dat");
	this->setscaleFCF(				gridPath + "scaleFCF.dat");
	this->setscaleOffsetCF(			gridPath + "offsetVecCF.dat");
	this->setscaleOffsetFC(			gridPath + "offsetVecFC.dat");
	this->setgeomBoundaryBcQs(		gridPath + "geomBoundaryQs.dat");
	this->setgeomBoundaryBcValues(	gridPath + "geomBoundaryValues.dat");
	this->setinletBcQs(				gridPath + "inletBoundaryQs.dat");
	this->setinletBcValues(			gridPath + "inletBoundaryValues.dat");
	this->setoutletBcQs(			gridPath + "outletBoundaryQs.dat");
	this->setoutletBcValues(		gridPath + "outletBoundaryValues.dat");
	this->settopBcQs(				gridPath + "topBoundaryQs.dat");
	this->settopBcValues(			gridPath + "topBoundaryValues.dat");
	this->setbottomBcQs(			gridPath + "bottomBoundaryQs.dat");
	this->setbottomBcValues(		gridPath + "bottomBoundaryValues.dat");
	this->setfrontBcQs(				gridPath + "frontBoundaryQs.dat");
	this->setfrontBcValues(			gridPath + "frontBoundaryValues.dat");
	this->setbackBcQs(				gridPath + "backBoundaryQs.dat");
	this->setbackBcValues(			gridPath + "backBoundaryValues.dat");
	this->setnumberNodes(			gridPath + "numberNodes.dat");
	this->setLBMvsSI(				gridPath + "LBMvsSI.dat");
	this->setmeasurePoints(			gridPath + "measurePoints.dat");
	this->setpropellerValues(		gridPath + "propellerValues.dat");
	this->setcpTop(					gridPath + "cpTop.dat");
	this->setcpBottom(				gridPath + "cpBottom.dat");
	this->setcpBottom2(				gridPath + "cpBottom2.dat");
	this->setConcentration(			gridPath + "conc.dat");
	this->setStreetVelocity(		gridPath + "streetVector.dat");
	//////////////////////////////////////////////////////////////////////////
	//Normals - Geometry
	this->setgeomBoundaryNormalX(gridPath + "geomBoundaryNormalX.dat");
	this->setgeomBoundaryNormalY(gridPath + "geomBoundaryNormalY.dat");
	this->setgeomBoundaryNormalZ(gridPath + "geomBoundaryNormalZ.dat");
	//Normals - Inlet
	this->setInflowBoundaryNormalX(gridPath + "inletBoundaryNormalX.dat");
	this->setInflowBoundaryNormalY(gridPath + "inletBoundaryNormalY.dat");
	this->setInflowBoundaryNormalZ(gridPath + "inletBoundaryNormalZ.dat");
	//Normals - Outlet
	this->setOutflowBoundaryNormalX(gridPath + "outletBoundaryNormalX.dat");
	this->setOutflowBoundaryNormalY(gridPath + "outletBoundaryNormalY.dat");
	this->setOutflowBoundaryNormalZ(gridPath + "outletBoundaryNormalZ.dat");
	//////////////////////////////////////////////////////////////////////////
	// //Forcing
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
	//quadricLimiters
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
	//////////////////////////////////////////////////////////////////////////
	//Particles
	if (configData.contains("calcParticles"))
		this->setCalcParticles(configData.getValue<bool>("calcParticles"));

	if (configData.contains("baseLevel"))
		this->setParticleBasicLevel(configData.getValue<int>("baseLevel"));

	if (configData.contains("initLevel"))
		this->setParticleInitLevel(configData.getValue<int>("initLevel"));

	if (configData.contains("numberOfParticles"))
		this->setNumberOfParticles(configData.getValue<int>("numberOfParticles"));

	if (configData.contains("startXHotWall"))
		this->setEndXHotWall(configData.getValue<real>("startXHotWall"));

	if (configData.contains("endXHotWall"))
		this->setCalcParticles(configData.getValue<real>("endXHotWall"));
	//////////////////////////////////////////////////////////////////////////
	//for Multi GPU
	if (this->getNumprocs() > 1)
	{
		//////////////////////////////////////////////////////////////////////////
		//3D domain decomposition
		std::vector<std::string> sendProcNeighborsX, sendProcNeighborsY, sendProcNeighborsZ;
		std::vector<std::string> recvProcNeighborsX, recvProcNeighborsY, recvProcNeighborsZ;
		for (int i = 0; i < this->getNumprocs(); i++)
		{
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
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Restart
	if (configData.contains("TimeDoCheckPoint"))
		this->setTimeDoCheckPoint(configData.getValue<uint>("TimeDoCheckPoint"));

	if (configData.contains("TimeDoRestart"))
		this->setTimeDoRestart(configData.getValue<uint>("TimeDoRestart"));

	if (configData.contains("DoCheckPoint"))
		this->setDoCheckPoint(configData.getValue<bool>("DoCheckPoint"));

	if (configData.contains("DoRestart"))
		this->setDoRestart(configData.getValue<bool>("DoRestart"));
	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (configData.contains("NOGL"))
		this->setMaxLevel(configData.getValue<int>("NOGL"));
	
	this->setGridX(std::vector<int>(this->getMaxLevel()+1, 32));
	this->setGridY(std::vector<int>(this->getMaxLevel()+1, 32));
	this->setGridZ(std::vector<int>(this->getMaxLevel()+1, 32));

	this->setDistX(std::vector<int>(this->getMaxLevel()+1, 32));
	this->setDistY(std::vector<int>(this->getMaxLevel()+1, 32));
	this->setDistZ(std::vector<int>(this->getMaxLevel()+1, 32));

	this->setNeedInterface(std::vector<bool>(6, true));

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

	
	if (configData.contains("NeedInterface"))
		this->setNeedInterface(configData.getVector<bool>("NeedInterface"));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Kernel
	if (configData.contains("MainKernelName"))
		this->setMainKernel(configData.getValue<std::string>("MainKernelName"));

	if (configData.contains("MultiKernelOn"))
		this->setMultiKernelOn(configData.getValue<bool>("MultiKernelOn"));

	if (configData.contains("MultiKernelLevel"))
		this->setMultiKernelLevel(configData.getVector<int>("MultiKernelLevel"));
	else if (this->getMultiKernelOn())
	{
		std::vector<int> tmp;
		for (int i = 0; i < this->getMaxLevel()+1; i++)
		{
			tmp.push_back(i);
		}
		this->setMultiKernelLevel(tmp);
	} 

	if (configData.contains("MultiKernelName"))
		this->setMultiKernel(StringUtil::toStringVector(configData.getValue<std::string>("MultiKernelName")));
	else if (this->getMultiKernelOn())
	{
        std::vector<std::string> tmp;
		for (int i = 0; i < this->getMaxLevel()+1; i++)
		{
			tmp.push_back("CumulantK17Comp");
		}
		this->setMultiKernel(tmp);
	}
}

Parameter::Parameter(SPtr<ConfigData> configData, vf::gpu::Communicator* comm)
{
	//////////////////////////////////////////////////////////////////////////
	this->setNumprocs(comm->getNummberOfProcess());
	this->setMyID(comm->getPID());
	//////////////////////////////////////////////////////////////////////////
	if (configData->isNumberOfDevicesInConfigFile())
		this->setMaxDev(configData->getNumberOfDevices());
	else
		this->setMaxDev((int)1);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isDevicesInConfigFile())
		this->setDevices(configData->getDevices());
	else
		this->setDevices(std::vector<uint>{(uint)0});
	//////////////////////////////////////////////////////////////////////////
	if (configData->isOutputPathInConfigFile())
		this->setOutputPath(configData->getOutputPath());
	else
		this->setOutputPath("C:/Output/");
	//////////////////////////////////////////////////////////////////////////
	if (configData->isPrefixInConfigFile())
		this->setOutputPrefix(configData->getPrefix());
	else
		this->setOutputPrefix("MyFile");
	//////////////////////////////////////////////////////////////////////////
	if (configData->isPrintOutputFilesInConfigFile())
		this->setPrintFiles(configData->getPrintOutputFiles());
	else
		this->setPrintFiles(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isGeometryValuesInConfigFile())
		this->setGeometryValues(configData->getGeometryValues());
	else
		this->setGeometryValues(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isCalc2ndOrderMomentsInConfigFile())
		this->setCalc2ndOrderMoments(configData->getCalc2ndOrderMoments());
	else
		this->setCalc2ndOrderMoments(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isCalc3rdOrderMomentsInConfigFile())
		this->setCalc3rdOrderMoments(configData->getCalc3rdOrderMoments());
	else
		this->setCalc3rdOrderMoments(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isCalcHighOrderMomentsInConfigFile())
		this->setCalcHighOrderMoments(configData->getCalcHighOrderMoments());
	else
		this->setCalcHighOrderMoments(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isCalcMedianInConfigFile())
		this->setCalcMedian(configData->getCalcMedian());
	else
		this->setCalcMedian(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isCalcDragLiftInConfigFile())
		this->setCalcDragLift(configData->getCalcDragLift());
	else
		this->setCalcDragLift(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isCalcCpInConfigFile())
		this->setCalcCp(configData->getCalcCp());
	else
		this->setCalcCp(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isWriteVeloASCIIfilesInConfigFile())
		this->setWriteVeloASCIIfiles(configData->getWriteVeloASCIIfiles());
	else
		this->setWriteVeloASCIIfiles(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isCalcPlaneConcInConfigFile())
		this->setCalcPlaneConc(configData->getCalcPlaneConc());
	else
		this->setCalcPlaneConc(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isConcFileInConfigFile())
		this->setConcFile(configData->getConcFile());
	else
		this->setConcFile(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isStreetVelocityFileInConfigFile())
		this->setStreetVelocityFile(configData->getStreetVelocityFile());
	else
		this->setStreetVelocityFile(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isUseMeasurePointsInConfigFile())
		this->setUseMeasurePoints(configData->getUseMeasurePoints());
	else
		this->setUseMeasurePoints(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isUseWaleInConfigFile())
		this->setUseWale(configData->getUseWale());
	else
		this->setUseWale(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isUseInitNeqInConfigFile())
		this->setUseInitNeq(configData->getUseInitNeq());
	else
		this->setUseInitNeq(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isSimulatePorousMediaInConfigFile())
		this->setSimulatePorousMedia(configData->getSimulatePorousMedia());
	else
		this->setSimulatePorousMedia(false);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isD3QxxInConfigFile())
		this->setD3Qxx(configData->getD3Qxx());
	else
		this->setD3Qxx((int)27);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isTEndInConfigFile())
		this->setTEnd(configData->getTEnd());
	else
		this->setTEnd((uint)10);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isTOutInConfigFile())
		this->setTOut(configData->getTOut());
	else
		this->setTOut((uint)1);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isTStartOutInConfigFile())
		this->setTStartOut(configData->getTStartOut());
	else
		this->setTStartOut((uint)0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isTimeCalcMedStartInConfigFile())
		this->setTimeCalcMedStart(configData->getTimeCalcMedStart());
	else
		this->setTimeCalcMedStart((int)0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isTimeCalcMedEndInConfigFile())
		this->setTimeCalcMedEnd(configData->getTimeCalcMedEnd());
	else
		this->setTimeCalcMedEnd((int)10);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isPressInIDInConfigFile())
		this->setPressInID(configData->getPressInID());
	else
		this->setPressInID((uint)0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isPressOutIDInConfigFile())
		this->setPressOutID(configData->getPressOutID());
	else
		this->setPressOutID((uint)0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isPressInZInConfigFile())
		this->setPressInZ(configData->getPressInZ());
	else
		this->setPressInZ((uint)1);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isPressOutZInConfigFile())
		this->setPressOutZ(configData->getPressOutZ());
	else
		this->setPressOutZ((uint)2);
	//////////////////////////////////////////////////////////////////////////
	//second component
	if (configData->isDiffOnInConfigFile())
		this->setDiffOn(configData->getDiffOn());
	else
		this->setDiffOn(false);

	if (configData->isDiffModInConfigFile())
		this->setDiffMod(configData->getDiffMod());
	else
		this->setDiffMod((int)27);

	if (configData->isDiffusivityInConfigFile())
		this->setDiffusivity(configData->getDiffusivity());
	else
		this->setDiffusivity((real)0.001);

	if (configData->isTemperatureInitInConfigFile())
		this->setTemperatureInit(configData->getTemperatureInit());
	else
		this->setTemperatureInit((real)0.0);

	if (configData->isTemperatureBCInConfigFile())
		this->setTemperatureBC(configData->getTemperatureBC());
	else
		this->setTemperatureBC((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isViscosityInConfigFile())
		this->setViscosity(configData->getViscosity());
	else
		this->setViscosity((real)0.001);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isVelocityInConfigFile())
		this->setVelocity(configData->getVelocity());
	else
		this->setVelocity((real)0.01);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isViscosityRatioInConfigFile())
		this->setViscosityRatio(configData->getViscosityRatio());
	else
		this->setViscosityRatio((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isVelocityRatioInConfigFile())
		this->setVelocityRatio(configData->getVelocityRatio());
	else
		this->setVelocityRatio((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isDensityRatioInConfigFile())
		this->setDensityRatio(configData->getDensityRatio());
	else
		this->setDensityRatio((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isPressRatioInConfigFile())
		this->setPressRatio(configData->getPressRatio());
	else
		this->setPressRatio((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isRealXInConfigFile())
		this->setRealX(configData->getRealX());
	else
		this->setRealX((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isRealYInConfigFile())
		this->setRealY(configData->getRealY());
	else
		this->setRealY((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	if (configData->isFactorPressBCInConfigFile())
		this->setFactorPressBC(configData->getFactorPressBC());
	else
		this->setFactorPressBC((real)1.0);
	//////////////////////////////////////////////////////////////////////////
	//read Geometry (STL)
	if (configData->isReadGeoInConfigFile())
		this->setReadGeo(configData->getReadGeo());
	else
		this->setReadGeo(false);

	if (configData->isGeometryFileCInConfigFile())
		this->setGeometryFileC(configData->getGeometryFileC());
	else if (this->getReadGeo())
	{
		std::cout << "GeometryFileC has to be defined!" << std::endl;
		exit(1);
	}
	else
		this->setGeometryFileC("");

	if (configData->isGeometryFileMInConfigFile())
		this->setGeometryFileM(configData->getGeometryFileM());
	else if (this->getReadGeo())
	{
		std::cout << "GeometryFileM has to be defined!" << std::endl;
		exit(1);
	}
	else
		this->setGeometryFileM("");

	if (configData->isGeometryFileFInConfigFile())
		this->setGeometryFileF(configData->getGeometryFileF());
	else if (this->getReadGeo())
	{
		std::cout << "GeometryFileF has to be defined!" << std::endl;
		exit(1);
	}
	else
		this->setGeometryFileF("");
	//////////////////////////////////////////////////////////////////////////
	if (configData->isClockCycleForMPInConfigFile())
		this->setclockCycleForMP(configData->getClockCycleForMP());
	else
		this->setclockCycleForMP((real)1.0);

	if (configData->isTimestepForMPInConfigFile())
		this->settimestepForMP(configData->getTimestepForMP());
	else
		this->settimestepForMP((uint)10);
	//////////////////////////////////////////////////////////////////////////
	std::string gridPath = "";
	if (configData->isGridPathInConfigFile())
		gridPath = configData->getGridPath();
	else
	{
		std::cout << "GridPath has to be defined!" << std::endl;
		exit(1);
	}

	if (this->getNumprocs() == 1)
		gridPath += "/";
	else
		gridPath += "/" + StringUtil::toString(this->getMyID()) + "/";
	//////////////////////////////////////////////////////////////////////////
	this->setFName(this->getOutputPath() + "/" + this->getOutputPrefix());
	//////////////////////////////////////////////////////////////////////////
	this->setgeoVec(				gridPath + "geoVec.dat");
	this->setcoordX(				gridPath + "coordX.dat");
	this->setcoordY(				gridPath + "coordY.dat");
	this->setcoordZ(				gridPath + "coordZ.dat");
	this->setneighborX(				gridPath + "neighborX.dat");
	this->setneighborY(				gridPath + "neighborY.dat");
	this->setneighborZ(				gridPath + "neighborZ.dat");
	this->setneighborWSB(			gridPath + "neighborWSB.dat");
	this->setscaleCFC(				gridPath + "scaleCFC.dat");
	this->setscaleCFF(				gridPath + "scaleCFF.dat");
	this->setscaleFCC(				gridPath + "scaleFCC.dat");
	this->setscaleFCF(				gridPath + "scaleFCF.dat");
	this->setscaleOffsetCF(			gridPath + "offsetVecCF.dat");
	this->setscaleOffsetFC(			gridPath + "offsetVecFC.dat");
	this->setgeomBoundaryBcQs(		gridPath + "geomBoundaryQs.dat");
	this->setgeomBoundaryBcValues(	gridPath + "geomBoundaryValues.dat");
	this->setinletBcQs(				gridPath + "inletBoundaryQs.dat");
	this->setinletBcValues(			gridPath + "inletBoundaryValues.dat");
	this->setoutletBcQs(			gridPath + "outletBoundaryQs.dat");
	this->setoutletBcValues(		gridPath + "outletBoundaryValues.dat");
	this->settopBcQs(				gridPath + "topBoundaryQs.dat");
	this->settopBcValues(			gridPath + "topBoundaryValues.dat");
	this->setbottomBcQs(			gridPath + "bottomBoundaryQs.dat");
	this->setbottomBcValues(		gridPath + "bottomBoundaryValues.dat");
	this->setfrontBcQs(				gridPath + "frontBoundaryQs.dat");
	this->setfrontBcValues(			gridPath + "frontBoundaryValues.dat");
	this->setbackBcQs(				gridPath + "backBoundaryQs.dat");
	this->setbackBcValues(			gridPath + "backBoundaryValues.dat");
	this->setnumberNodes(			gridPath + "numberNodes.dat");
	this->setLBMvsSI(				gridPath + "LBMvsSI.dat");
	this->setmeasurePoints(			gridPath + "measurePoints.dat");
	this->setpropellerValues(		gridPath + "propellerValues.dat");
	this->setcpTop(					gridPath + "cpTop.dat");
	this->setcpBottom(				gridPath + "cpBottom.dat");
	this->setcpBottom2(				gridPath + "cpBottom2.dat");
	this->setConcentration(			gridPath + "conc.dat");
	this->setStreetVelocity(		gridPath + "streetVector.dat");
	//////////////////////////////////////////////////////////////////////////
	//Normals - Geometry
	this->setgeomBoundaryNormalX(gridPath + "geomBoundaryNormalX.dat");
	this->setgeomBoundaryNormalY(gridPath + "geomBoundaryNormalY.dat");
	this->setgeomBoundaryNormalZ(gridPath + "geomBoundaryNormalZ.dat");
	//Normals - Inlet
	this->setInflowBoundaryNormalX(gridPath + "inletBoundaryNormalX.dat");
	this->setInflowBoundaryNormalY(gridPath + "inletBoundaryNormalY.dat");
	this->setInflowBoundaryNormalZ(gridPath + "inletBoundaryNormalZ.dat");
	//Normals - Outlet
	this->setOutflowBoundaryNormalX(gridPath + "outletBoundaryNormalX.dat");
	this->setOutflowBoundaryNormalY(gridPath + "outletBoundaryNormalY.dat");
	this->setOutflowBoundaryNormalZ(gridPath + "outletBoundaryNormalZ.dat");
	//////////////////////////////////////////////////////////////////////////
	//Forcing
	real forcingX = 0.0;
	real forcingY = 0.0;
	real forcingZ = 0.0;

	if (configData->isForcingXInConfigFile())
		forcingX = configData->getForcingX();
	if (configData->isForcingYInConfigFile())
		forcingY = configData->getForcingY();
	if (configData->isForcingZInConfigFile())
		forcingZ = configData->getForcingZ();

	this->setForcing(forcingX, forcingY, forcingZ);
	//////////////////////////////////////////////////////////////////////////
	//quadricLimiters
	real quadricLimiterP = (real)0.01;
	real quadricLimiterM = (real)0.01;
	real quadricLimiterD = (real)0.01;

	if (configData->isQuadricLimiterPInConfigFile())
		quadricLimiterP = configData->getQuadricLimiterP();
	if (configData->isQuadricLimiterMInConfigFile())
		quadricLimiterM = configData->getQuadricLimiterM();
	if (configData->isQuadricLimiterDInConfigFile())
		quadricLimiterD = configData->getQuadricLimiterD();

	this->setQuadricLimiters(quadricLimiterP, quadricLimiterM, quadricLimiterD);
	//////////////////////////////////////////////////////////////////////////
	//Particles
	if (configData->isCalcParticlesInConfigFile())
		this->setCalcParticles(configData->getCalcParticles());
	else
		this->setCalcParticles(false);

	if (configData->isParticleBasicLevelInConfigFile())
		this->setParticleBasicLevel(configData->getParticleBasicLevel());
	else
		this->setParticleBasicLevel((int)0);

	if (configData->isParticleInitLevelInConfigFile())
		this->setParticleInitLevel(configData->getParticleInitLevel());
	else
		this->setParticleInitLevel((int)0);

	if (configData->isNumberOfParticlesInConfigFile())
		this->setNumberOfParticles(configData->getNumberOfParticles());
	else
		this->setNumberOfParticles((int)0);

	if (configData->isStartXHotWallInConfigFile())
		this->setStartXHotWall(configData->getStartXHotWall());
	else
		this->setStartXHotWall((real)0);

	if (configData->isEndXHotWallInConfigFile())
		this->setEndXHotWall(configData->getEndXHotWall());
	else
		this->setEndXHotWall((real)0);
	//////////////////////////////////////////////////////////////////////////
	//for Multi GPU
	if (this->getNumprocs() > 1)
	{
		//////////////////////////////////////////////////////////////////////////
		//3D domain decomposition
		std::vector<std::string> sendProcNeighborsX, sendProcNeighborsY, sendProcNeighborsZ;
		std::vector<std::string> recvProcNeighborsX, recvProcNeighborsY, recvProcNeighborsZ;
		for (int i = 0; i < this->getNumprocs(); i++)
		{
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
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Restart
	if (configData->isTimeDoCheckPointInConfigFile())
		this->setTimeDoCheckPoint(configData->getTimeDoCheckPoint());
	else
		this->setTimeDoCheckPoint((uint)0);

	if (configData->isTimeDoRestartInConfigFile())
		this->setTimeDoRestart(configData->getTimeDoRestart());
	else
		this->setTimeDoRestart((uint)0);

	if (configData->isDoCheckPointInConfigFile())
		this->setDoCheckPoint(configData->getDoCheckPoint());
	else
		this->setDoCheckPoint(false);

	if (configData->isDoRestartInConfigFile())
		this->setDoRestart(configData->getDoRestart());
	else
		this->setDoRestart(false);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (configData->isMaxLevelInConfigFile())
		this->setMaxLevel(configData->getMaxLevel());
	else
		this->setMaxLevel((int)1);

	if (configData->isGridXInConfigFile())
		this->setGridX(configData->getGridX());
	else
		this->setGridX(std::vector<int>(this->getMaxLevel()+1, 32));

	if (configData->isGridYInConfigFile())
		this->setGridY(configData->getGridY());
	else
		this->setGridY(std::vector<int>(this->getMaxLevel()+1, 32));

	if (configData->isGridZInConfigFile())
		this->setGridZ(configData->getGridZ());
	else
		this->setGridZ(std::vector<int>(this->getMaxLevel()+1, 32));

	if (configData->isDistXInConfigFile())
		this->setDistX(configData->getDistX());
	else
		this->setDistX(std::vector<int>(this->getMaxLevel()+1, 32));

	if (configData->isDistYInConfigFile())
		this->setDistY(configData->getDistY());
	else
		this->setDistY(std::vector<int>(this->getMaxLevel()+1, 32));

	if (configData->isDistZInConfigFile())
		this->setDistZ(configData->getDistZ());
	else
		this->setDistZ(std::vector<int>(this->getMaxLevel()+1, 32));

	if (configData->isNeedInterfaceInConfigFile())
		this->setNeedInterface(configData->getNeedInterface());
	else
		this->setNeedInterface(std::vector<bool>(6, true));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Kernel
	if (configData->isMainKernelInConfigFile())
		this->setMainKernel(configData->getMainKernel());
	else
		this->setMainKernel("CumulantK15Comp");

	if (configData->isMultiKernelOnInConfigFile())
		this->setMultiKernelOn(configData->getMultiKernelOn());
	else
		this->setMultiKernelOn(false);

	if (configData->isMultiKernelLevelInConfigFile())
		this->setMultiKernelLevel(configData->getMultiKernelLevel());
	else if (this->getMultiKernelOn())
	{
		std::vector<int> tmp;
		for (int i = 0; i < this->getMaxLevel()+1; i++)
		{
			tmp.push_back(i);
		}
		this->setMultiKernelLevel(tmp);
	} 
	else
		this->setMultiKernelLevel(std::vector<int>(0));

	if (configData->isMultiKernelNameInConfigFile()) {
        std::vector<std::string> kernels;
		for (std::size_t i = 0; i < configData->getMultiKernelName().size(); i++) {
			kernels.push_back(configData->getMultiKernelName().at(i));
		}
		this->setMultiKernel(kernels);
	}
	else if (this->getMultiKernelOn())
	{
        std::vector<std::string> tmp;
		for (int i = 0; i < this->getMaxLevel()+1; i++)
		{
			tmp.push_back("CumulantK15Comp");
		}
		this->setMultiKernel(tmp);
	}
	else {
        std::vector<std::string> tmp;
		this->setMultiKernel(tmp);
	}		
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//init-method
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		parH[i]->numberofthreads       = 64;// 128;
		parH[i]->gridNX                = getGridX().at(i);
		parH[i]->gridNY                = getGridY().at(i);
		parH[i]->gridNZ                = getGridZ().at(i);
		parH[i]->vis                   = ic.vis*pow(2.f,i);
		parH[i]->diffusivity           = ic.Diffusivity*pow(2.f,i);
		parH[i]->omega                 = 1.0f/(3.0f*parH[i]->vis+0.5f);//omega :-) not s9 = -1.0f/(3.0f*parH[i]->vis+0.5f);//
		parH[i]->nx                    = parH[i]->gridNX + 2 * STARTOFFX;
		parH[i]->ny                    = parH[i]->gridNY + 2 * STARTOFFY;
		parH[i]->nz                    = parH[i]->gridNZ + 2 * STARTOFFZ;
		parH[i]->size_Mat              = parH[i]->nx * parH[i]->ny * parH[i]->nz;
		parH[i]->sizePlaneXY           = parH[i]->nx * parH[i]->ny;
		parH[i]->sizePlaneYZ           = parH[i]->ny * parH[i]->nz;
		parH[i]->sizePlaneXZ           = parH[i]->nx * parH[i]->nz;
		parH[i]->mem_size_real         = sizeof(real     ) * parH[i]->size_Mat;
		parH[i]->mem_size_int          = sizeof(unsigned int) * parH[i]->size_Mat;
		parH[i]->mem_size_bool         = sizeof(bool        ) * parH[i]->size_Mat;
		parH[i]->mem_size_real_yz      = sizeof(real     ) * parH[i]->ny * parH[i]->nz;
		parH[i]->evenOrOdd             = true;
		parH[i]->startz                = parH[i]->gridNZ * ic.myid;
		parH[i]->endz                  = parH[i]->gridNZ * ic.myid + parH[i]->gridNZ;
		parH[i]->Lx                    = (real)((1.f*parH[i]->gridNX - 1.f)/(pow(2.f,i)));
		parH[i]->Ly                    = (real)((1.f*parH[i]->gridNY - 1.f)/(pow(2.f,i)));
		parH[i]->Lz                    = (real)((1.f*parH[i]->gridNZ - 1.f)/(pow(2.f,i)));
		parH[i]->dx                    = (real)(1.f/(pow(2.f,i)));
		parH[i]->XdistKn               = getDistX().at(i);
		parH[i]->YdistKn               = getDistY().at(i);
		parH[i]->ZdistKn               = getDistZ().at(i);
		if (i==coarse)
		{
			parH[i]->distX                 = (real)getDistX().at(i);
			parH[i]->distY                 = (real)getDistY().at(i);
			parH[i]->distZ                 = (real)getDistZ().at(i);
			parH[i]->mTtoWx                = (real)1.0f;
			parH[i]->mTtoWy                = (real)1.0f;
			parH[i]->mTtoWz                = (real)1.0f;
			parH[i]->cTtoWx                = (real)0.0f;
			parH[i]->cTtoWy                = (real)0.0f;
			parH[i]->cTtoWz                = (real)0.0f;
			////MGs Trafo///////////////////////////////////////////////////////////////
			//parH[i]->cStartx               = (real)parH[i]->XdistKn;
			//parH[i]->cStarty               = (real)parH[i]->XdistKn;
			//parH[i]->cStartz               = (real)parH[i]->XdistKn;
			////////////////////////////////////////////////////////////////////////////
		} 
		else
		{
			//Geller
			parH[i]->distX                 = ((real)getDistX().at(i) + 0.25f) * parH[i-1]->dx;
			parH[i]->distY                 = ((real)getDistY().at(i) + 0.25f) * parH[i-1]->dx;
			parH[i]->distZ                 = ((real)getDistZ().at(i) + 0.25f) * parH[i-1]->dx;
			//parH[i]->distX                 = ((real)getDistX().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distX;
			//parH[i]->distY                 = ((real)getDistY().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distY;
			//parH[i]->distZ                 = ((real)getDistZ().at(i) + 0.25f) * parH[i-1]->dx + parH[i-1]->distZ;
			parH[i]->mTtoWx                = (real)pow(0.5f,i);
			parH[i]->mTtoWy                = (real)pow(0.5f,i);
			parH[i]->mTtoWz                = (real)pow(0.5f,i);
			parH[i]->cTtoWx                = (real)(STARTOFFX/2.f + (parH[i]->gridNX+1.f)/4.f); //funzt nur fuer zwei level
			parH[i]->cTtoWy                = (real)(STARTOFFY/2.f + (parH[i]->gridNY+1.f)/4.f); //funzt nur fuer zwei level
			parH[i]->cTtoWz                = (real)(STARTOFFZ/2.f + (parH[i]->gridNZ+1.f)/4.f); //funzt nur fuer zwei level
			////MGs Trafo///////////////////////////////////////////////////////////////
			//parH[i]->cStartx               = (real)parH[i]->XdistKn;
			//parH[i]->cStarty               = (real)parH[i]->XdistKn;
			//parH[i]->cStartz               = (real)parH[i]->XdistKn;
			////////////////////////////////////////////////////////////////////////////
		}
		parH[i]->need_interface[INTERFACE_E]=getNeedInterface().at(INTERFACE_E);
		parH[i]->need_interface[INTERFACE_W]=getNeedInterface().at(INTERFACE_W);
		parH[i]->need_interface[INTERFACE_N]=getNeedInterface().at(INTERFACE_N);
		parH[i]->need_interface[INTERFACE_S]=getNeedInterface().at(INTERFACE_S);
		parH[i]->need_interface[INTERFACE_T]=getNeedInterface().at(INTERFACE_T);
		parH[i]->need_interface[INTERFACE_B]=getNeedInterface().at(INTERFACE_B);
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
		parD[i]->diffusivity           = parH[i]->diffusivity;
		parD[i]->omega                 = parH[i]->omega;
		parD[i]->nx                    = parH[i]->nx;
		parD[i]->ny                    = parH[i]->ny;
		parD[i]->nz                    = parH[i]->nz;
		parD[i]->size_Mat              = parH[i]->size_Mat;
		parD[i]->sizePlaneXY           = parH[i]->sizePlaneXY;
		parD[i]->sizePlaneYZ           = parH[i]->sizePlaneYZ;
		parD[i]->sizePlaneXZ           = parH[i]->sizePlaneXZ;
		parD[i]->mem_size_real         = sizeof(real     ) * parD[i]->size_Mat;
		parD[i]->mem_size_int          = sizeof(unsigned int) * parD[i]->size_Mat;
		parD[i]->mem_size_bool         = sizeof(bool        ) * parD[i]->size_Mat;
		parD[i]->mem_size_real_yz      = sizeof(real     ) * parD[i]->ny * parD[i]->nz;
		parD[i]->evenOrOdd             = parH[i]->evenOrOdd;
		parD[i]->startz                = parH[i]->startz;
		parD[i]->endz                  = parH[i]->endz;
		parD[i]->Lx                    = parH[i]->Lx;
		parD[i]->Ly                    = parH[i]->Ly;
		parD[i]->Lz                    = parH[i]->Lz;
		parD[i]->dx                    = parH[i]->dx;
		parD[i]->XdistKn               = parH[i]->XdistKn;
		parD[i]->YdistKn               = parH[i]->YdistKn;
		parD[i]->ZdistKn               = parH[i]->ZdistKn;
		parD[i]->distX                 = parH[i]->distX;
		parD[i]->distY                 = parH[i]->distY;
		parD[i]->distZ                 = parH[i]->distZ;
	}

	//Interface
	//comment out for geller
	//for (int i = coarse; i < fine; i++)
	//{
	//   initInterfaceParameter(i);
	//}
}
void Parameter::setSizeMatSparse(int level)
{
	parH[level]->size_Mat_SP = 1;
	parD[level]->size_Mat_SP = 1;
	parH[level]->sizePlaneSB = 0;
	parH[level]->sizePlaneST = 0;
	parH[level]->sizePlaneRB = 0;
	parH[level]->sizePlaneRT = 0;
	parH[level]->isSetSendB  = false;
	parH[level]->isSetSendT  = false;
	parH[level]->isSetRecvB  = false;
	parH[level]->isSetRecvT  = false;
	unsigned int mm[8];

	for (unsigned int k=1; k<parH[level]->gridNZ + 2 * STARTOFFZ - 1; k++)
	{
		for (unsigned int j=1; j<parH[level]->gridNY + 2 * STARTOFFY - 1; j++)
		{
			for (unsigned int i=1; i<parH[level]->gridNX + 2 * STARTOFFX - 1; i++)
			{
				mm[0]= parH[level]->nx*(parH[level]->ny*k + j) + i;
				mm[1]= mm[0]                                                  -1; //W
				mm[2]= mm[0]                                  -parH[level]->nx-1; //SW
				mm[3]= mm[0]                                  -parH[level]->nx;   //S
				mm[4]= mm[0]-(parH[level]->nx*parH[level]->ny);                   //B
				mm[5]= mm[0]-(parH[level]->nx*parH[level]->ny)                -1; //BW
				mm[6]= mm[0]-(parH[level]->nx*parH[level]->ny)-parH[level]->nx;   //BS
				mm[7]= mm[0]-(parH[level]->nx*parH[level]->ny)-parH[level]->nx-1; //BSW

				if ( parH[level]->geo[mm[0]] != GEO_VOID ||
					parH[level]->geo[mm[1]] != GEO_VOID ||
					parH[level]->geo[mm[2]] != GEO_VOID ||
					parH[level]->geo[mm[3]] != GEO_VOID ||
					parH[level]->geo[mm[4]] != GEO_VOID ||
					parH[level]->geo[mm[5]] != GEO_VOID ||
					parH[level]->geo[mm[6]] != GEO_VOID ||
					parH[level]->geo[mm[7]] != GEO_VOID )
				{
					//////////////////////////////////////////////////////////////////////////
					//add some stuff for the data exchange between the GPUs //////////////////
					if (k == STARTOFFZ)
					{
						parH[level]->sizePlaneSB  += 1;
						if (parH[level]->isSetSendB == false)
						{
							parH[level]->startB = mm[0];
							parH[level]->isSetSendB = true;
						}
					} 
					else if (k == parH[level]->gridNZ + STARTOFFZ - 1)
					{
						parH[level]->sizePlaneST  += 1;
						if (parH[level]->isSetSendT == false)
						{
							parH[level]->startT = mm[0];
							parH[level]->isSetSendT = true;
						}
					}
					else if (k == parH[level]->gridNZ + STARTOFFZ)
					{
						parH[level]->sizePlaneRB  += 1;
						if (parH[level]->isSetRecvB == false)
						{
							parH[level]->endB = mm[0];
							parH[level]->isSetRecvB = true;
						}
					}
					else if (k == STARTOFFZ-1)
					{
						parH[level]->sizePlaneRT  += 1;
						if (parH[level]->isSetRecvT == false)
						{
							parH[level]->endT = mm[0];
							parH[level]->isSetRecvT = true;
						}
					}
					//////////////////////////////////////////////////////////////////////////
					parH[level]->k[mm[0]]    = parH[level]->size_Mat_SP;
					parH[level]->size_Mat_SP = parH[level]->size_Mat_SP + 1;               
					parD[level]->size_Mat_SP = parD[level]->size_Mat_SP + 1;  
				}
				else parH[level]->k[mm[0]] = 0;
			}
		}
	}
	parH[level]->mem_size_real_SP    = sizeof(real     ) * parH[level]->size_Mat_SP;
	parH[level]->mem_size_int_SP        = sizeof(unsigned int) * parH[level]->size_Mat_SP;
	parD[level]->mem_size_real_SP    = sizeof(real     ) * parD[level]->size_Mat_SP;
	parD[level]->mem_size_int_SP        = sizeof(unsigned int) * parD[level]->size_Mat_SP;
}
void Parameter::fillSparse(int level)
{
    //nsigned int li = ((parH[level]->gridNX+STARTOFFX-2)-(STARTOFFX+1)-1);
    //unsigned int lj = ((parH[level]->gridNY+STARTOFFY-2)-(STARTOFFY+1)-1);
	// real globalX, globalY, globalZ;

	real PI = 3.141592653589793238462643383279f;

	for (unsigned int k=1; k<parH[level]->gridNZ + 2 * STARTOFFZ - 1; k++)
	{
		for (unsigned int j=1; j<parH[level]->gridNY + 2 * STARTOFFY - 1; j++)
		{
			for (unsigned int i=1; i<parH[level]->gridNX + 2 * STARTOFFX - 1; i++)
			{
				int m = parH[level]->nx*(parH[level]->ny*k + j) + i;
				if ((k < parH[level]->gridNZ + 2 * STARTOFFZ - 2) && (j < parH[level]->gridNY + 2 * STARTOFFY - 2) && (i < parH[level]->gridNX + 2 * STARTOFFX - 2))
				{
					if ((X1PERIODIC == true) && (level==coarse) && (i==parH[level]->gridNX + STARTOFFX - 1)) 
					{
						int mm = parH[level]->nx*(parH[level]->ny*k + j) + STARTOFFX;
						parH[level]->neighborX_SP[parH[level]->k[m]] = parH[level]->k[mm];
					}
					else
					{
						parH[level]->neighborX_SP[parH[level]->k[m]] = parH[level]->k[m+1];
					}
					if ((X2PERIODIC == true) && (level==coarse) && (j==parH[level]->gridNY + STARTOFFY - 1)) 
					{
						int mm = parH[level]->nx*(parH[level]->ny*k + STARTOFFY) + i;
						parH[level]->neighborY_SP[parH[level]->k[m]] = parH[level]->k[mm];
					}
					else
					{
						parH[level]->neighborY_SP[parH[level]->k[m]] = parH[level]->k[m+parH[level]->nx];
					}
					if ((X3PERIODIC == true) && (level==coarse) && (k==parH[level]->gridNZ + STARTOFFZ - 1)) 
					{
						int mm = parH[level]->nx*(parH[level]->ny*STARTOFFZ + j) + i;
						parH[level]->neighborZ_SP[parH[level]->k[m]] = parH[level]->k[mm];
					}
					else
					{
						parH[level]->neighborZ_SP[parH[level]->k[m]] = parH[level]->k[m+(parH[level]->nx*parH[level]->ny)];
					}
				}
				parH[level]->geoSP[parH[level]->k[m]]        = parH[level]->geo[m];
				////////////////////////////////////////////////////////////////////////////
				////Coordinates
				//parH[level]->coordX_SP[parH[level]->k[m]]    = i;
				//parH[level]->coordY_SP[parH[level]->k[m]]    = j;
				//parH[level]->coordZ_SP[parH[level]->k[m]]    = k;
				////////////////////////////////////////////////////////////////////////////
				if (diffOn==true)
				{
					parH[level]->Conc[parH[level]->k[m]]         = parH[level]->Conc_Full[m];
				}
				////////////////////////////////////////////////////////////////////////////
				////set pressure in the middle of the fine grid
				//if (level == getFine())
				//{
				//   if(   i == parH[level]->gridNX/2 + STARTOFFX
				//      && j == parH[level]->gridNY/2 + STARTOFFY 
				//      && k == parH[level]->gridNZ/2 + STARTOFFZ) 
				//      parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.1f;             
				//   else 
				//      parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.0f;
				//} 
				//else
				//{
				//   parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.0f;
				//}
				// globalX = TrafoXtoWorld(i,level);
				// globalY = TrafoYtoWorld(j,level);
				// globalZ = TrafoZtoWorld(k,level);
				//without setting a pressure
				parH[level]->rho_SP[parH[level]->k[m]]       = (real)0.0f;       //parH[level]->Conc_Full[m];//bitte schnell wieder entfernen!!!
				//////////////////////////////////////////////////////////////////////////
				parH[level]->vx_SP[parH[level]->k[m]]        = (real)0.0f;
				//parH[level]->vx_SP[parH[level]->k[m]]        = u0/3.0;
				parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0f;
				//parH[level]->vy_SP[parH[level]->k[m]]        = u0/3.0;
				parH[level]->vz_SP[parH[level]->k[m]]        = (real)0.0f;
				//parH[level]->vz_SP[parH[level]->k[m]]        = u0/3.0;
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(u0*2.f)*((-4.f*globalX*globalX + parH[level]->gridNX*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + globalX*(-4.f + 4.f*parH[level]->gridNX + 8.f*STARTOFFX))*(-4.f*globalY*globalY + parH[level]->gridNY*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + globalY*(-4.f + 4.f*parH[level]->gridNY + 8.f*STARTOFFY)))/((2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNY)*(2.f - parH[level]->gridNY));
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(u0*2.f)*((-4.f*i*i + parH[level]->gridNX*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + i*(-4.f + 4.f*parH[level]->gridNX + 8.f*STARTOFFX))*(-4.f*j*j + parH[level]->gridNY*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + j*(-4.f + 4.f*parH[level]->gridNY + 8.f*STARTOFFY)))/((2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNX)*(2.f - parH[level]->gridNY)*(2.f - parH[level]->gridNY));
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(16.f*(u0*2.f)*(i-(STARTOFFX+1)-0.5f)*(li-1.5f-(i-(STARTOFFX+1)))*(j-(STARTOFFY+1)-0.5f)*(lj-1.5f-(j-(STARTOFFY+1))))/(li*lj*li*lj);//(16.f*(u0*2.f)*i*j*(parH[level]->nx-i)*(parH[level]->ny-j))/(parH[level]->nx*parH[level]->nx*parH[level]->ny*parH[level]->ny); //u0;
				//////////////////////////////////////////////////////////////////////////
				////gerade
				//parH[level]->vx_SP[parH[level]->k[m]] = (real)((32. * 32. * 3.) / (1000.*(real)parH[level]->gridNX));//(real)parH[level]->gridNX / (real)1000 * 3.0;
				//parH[level]->vy_SP[parH[level]->k[m]] = (real)((getVelocity() * sin(2.0 * i / parH[level]->gridNX * PI) * cos(2.0 * k / parH[level]->gridNZ * PI)) * (32. / (real)parH[level]->gridNX));
				//parH[level]->vz_SP[parH[level]->k[m]] = (real)0.0f;
				//schraeg x
				// 			parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNX) + (getVelocity() * cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));
				// 			parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				// 			parH[level]->vz_SP[parH[level]->k[m]]        = (real)(getVelocity() * cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				//schraeg z
				//parH[level]->vx_SP[parH[level]->k[m]]        = (real)(getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				//parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNZ) + (getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));

				  			//Taylor Green Vortex uniform
				  			parH[level]->rho_SP[parH[level]->k[m]]       = (real)((getVelocity()*getVelocity())*3.0/4.0*(cos((i)*4.0*PI/(real)parH[level]->gridNX)+cos((k)*4.0*PI/(real)parH[level]->gridNZ)))*(real)(parH[level]->gridNZ)/(real)(parH[level]->gridNX);
				  			//inkl. ueberlagerter Geschwindigkeit
				  // 			parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNX) + getVelocity()*sin(((i)*2.0*PI/(real)parH[level]->gridNX))*cos((k)*2.0*PI/(real)parH[level]->gridNZ));
				  			parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000. * 32.) * getVelocity() / 0.001 + getVelocity()*sin(((i)*2.0*PI/(real)parH[level]->gridNX))*cos((k)*2.0*PI/(real)parH[level]->gridNZ));
				  			//ohne ueberlagerter Geschwindigkeit
				  //			parH[level]->vx_SP[parH[level]->k[m]]        = (real)(getVelocity()*sin(((i)*2.0*PI/(real)parH[level]->gridNX))*cos((k)*2.0*PI/(real)parH[level]->gridNZ));
				  			parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				  			parH[level]->vz_SP[parH[level]->k[m]]        = (real)(-getVelocity()*cos(((i)*2.0*PI/(real)parH[level]->gridNX))*sin((k)*2.0*PI/(real)parH[level]->gridNZ))*(real)(parH[level]->gridNZ)/(real)(parH[level]->gridNX);            

				//Kernel Fix Test
				//parH[level]->vx_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNX) + (getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));
				//parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)(getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				////parH[level]->vx_SP[parH[level]->k[m]]        = (real)(getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI)));
				////parH[level]->vy_SP[parH[level]->k[m]]        = (real)0.0;
				////parH[level]->vz_SP[parH[level]->k[m]]        = (real)((32. * 32. * 3.)/(1000.*(real)parH[level]->gridNZ) + (getVelocity() * std::cos((2.0 * k / parH[level]->gridNZ * PI) + (2.0 * i / parH[level]->gridNX * PI))));
				//////////////////////////////////////////////////////////////////////////
				//Taylor Green Vortex
				//InitglobalX = TrafoXtoMGsWorld(i,level);
				//InitglobalY = TrafoYtoMGsWorld(j,level);
				//InitglobalZ = TrafoZtoMGsWorld(k,level);
				//parH[level]->rho_SP[parH[level]->k[m]]       = (real)((u0*u0)*3.f/4.f*(cos((InitglobalX)*4.f*PI/parH[level]->gridNX)+cos((InitglobalY)*4.f*PI/parH[level]->gridNY)));
				//parH[level]->vx_SP[parH[level]->k[m]]        = (real)( u0*sin(((InitglobalX)*2.f*PI/parH[level]->gridNX))*cos((InitglobalY)*2.f*PI/parH[level]->gridNY));
				//parH[level]->vy_SP[parH[level]->k[m]]        = (real)(-u0*cos(((InitglobalX)*2.f*PI/parH[level]->gridNX))*sin((InitglobalY)*2.f*PI/parH[level]->gridNY));
				//parH[level]->vz_SP[parH[level]->k[m]]        = (real)0.0f;            
				//////////////////////////////////////////////////////////////////////////
			}
		}
	}
	parH[level]->neighborX_SP[parH[level]->k[0]] = 0;
	parH[level]->neighborY_SP[parH[level]->k[0]] = 0;
	parH[level]->neighborZ_SP[parH[level]->k[0]] = 0;
	parH[level]->geoSP[       parH[level]->k[0]] = GEO_VOID;
	parH[level]->rho_SP[      parH[level]->k[0]] = (real)0.f;
	parH[level]->vx_SP[       parH[level]->k[0]] = (real)0.f;
	parH[level]->vy_SP[       parH[level]->k[0]] = (real)0.f;
	parH[level]->vz_SP[       parH[level]->k[0]] = (real)0.f;
	////////////////////////////////////////////////////////////////////////////
	////Coordinates
	//parH[level]->coordX_SP[parH[level]->k[0]]    = 0;
	//parH[level]->coordY_SP[parH[level]->k[0]]    = 0;
	//parH[level]->coordZ_SP[parH[level]->k[0]]    = 0;
	////////////////////////////////////////////////////////////////////////////
}
void Parameter::copyMeasurePointsArrayToVector(int lev)
{
	int valuesPerClockCycle = (int)(getclockCycleForMP()/getTimestepForMP());
	for(int i = 0; i < (int)parH[lev]->MP.size(); i++)
	{
		for(int j = 0; j < valuesPerClockCycle; j++)
		{
			int index = i*valuesPerClockCycle+j;
			parH[lev]->MP[i].Vx.push_back(parH[lev]->VxMP[index]);
			parH[lev]->MP[i].Vy.push_back(parH[lev]->VyMP[index]);
			parH[lev]->MP[i].Vz.push_back(parH[lev]->VzMP[index]);
			parH[lev]->MP[i].Rho.push_back(parH[lev]->RhoMP[index]);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//set-methods
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
void Parameter::setlimitOfNodesForVTK(unsigned int limitOfNodesForVTK)
{
	this->limitOfNodesForVTK = limitOfNodesForVTK;
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
void Parameter::setDiffMod(int DiffMod)
{
	diffMod = DiffMod;
}
void Parameter::setD3Qxx(int d3qxx)
{
	this->D3Qxx = d3qxx;
}
void Parameter::setMaxLevel(int maxlevel)
{
	this->maxlevel = maxlevel-1;
}
void Parameter::setParticleBasicLevel(int pbl)
{
	this->particleBasicLevel = pbl;
}
void Parameter::setParticleInitLevel(int pil)
{
	this->particleInitLevel = pil;
}
void Parameter::setNumberOfParticles(int nop)
{
	this->numberOfParticles = nop;
}
void Parameter::setCalcParticles(bool calcParticles)
{
	this->calcParticles = calcParticles;
}
void Parameter::setStartXHotWall(real startXHotWall)
{
	this->startXHotWall = startXHotWall;
}
void Parameter::setEndXHotWall(real endXHotWall)
{
	this->endXHotWall = endXHotWall;
}
void Parameter::setTEnd(unsigned int tend)
{
	ic.tend = tend;
}
void Parameter::setTOut(unsigned int tout)
{
	ic.tout = tout;
}
void Parameter::setTStartOut(unsigned int tStartOut)
{
	ic.tStartOut = tStartOut;
}
void Parameter::setTimestepOfCoarseLevel(unsigned int timestep)
{
	this->timestep = timestep;
}
void Parameter::setCalcMedian(bool calcMedian)
{
	ic.calcMedian = calcMedian;
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
	ic.tCalcMedStart = CalcMedStart;
}
void Parameter::setTimeCalcMedEnd(int CalcMedEnd)
{
	ic.tCalcMedEnd = CalcMedEnd;
}
void Parameter::setOutputPath(std::string oPath)
{
	ic.oPath = oPath;
}
void Parameter::setOutputPrefix(std::string oPrefix)
{
	//std::string test = fname;
	ic.oPrefix = oPrefix;
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
void Parameter::setReadGeo(bool readGeo)
{
	ic.readGeo = readGeo;
}
void Parameter::setDiffusivity(real Diffusivity)
{
	ic.Diffusivity = Diffusivity;
}
void Parameter::setTemperatureInit(real Temp)
{
	ic.Temp = Temp;
}
void Parameter::setTemperatureBC(real TempBC)
{
	ic.TempBC = TempBC;
}
void Parameter::setViscosity(real Viscosity)
{
	ic.vis = Viscosity;
}
void Parameter::setVelocity(real Velocity)
{
	ic.u0 = Velocity;
}
void Parameter::setViscosityRatio(real ViscosityRatio)
{
	ic.vis_ratio = ViscosityRatio;
}
void Parameter::setVelocityRatio(real VelocityRatio)
{
	ic.u0_ratio = VelocityRatio;
}
void Parameter::setDensityRatio(real DensityRatio)
{
	ic.delta_rho = DensityRatio;
}
void Parameter::setPressRatio(real PressRatio)
{
	ic.delta_press = PressRatio;
}
void Parameter::setRealX(real RealX)
{
	ic.RealX = RealX;
}
void Parameter::setRealY(real RealY)
{
	ic.RealY = RealY;
}
void Parameter::setPressInID(unsigned int PressInID)
{
	ic.PressInID = PressInID;
}
void Parameter::setPressOutID(unsigned int PressOutID)
{
	ic.PressOutID = PressOutID;
}
void Parameter::setPressInZ(unsigned int PressInZ)
{
	ic.PressInZ = PressInZ;
}
void Parameter::setPressOutZ(unsigned int PressOutZ)
{
	ic.PressOutZ = PressOutZ;
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
void Parameter::setDevices(std::vector<uint> devices)
{
	ic.devices = devices;
}
void Parameter::setGeometryFileC(std::string GeometryFileC)
{
	ic.geometryFileC = GeometryFileC;
}
void Parameter::setGeometryFileM(std::string GeometryFileM)
{
	ic.geometryFileM = GeometryFileM;
}
void Parameter::setGeometryFileF(std::string GeometryFileF)
{
	ic.geometryFileF = GeometryFileF;
}
void Parameter::setNeedInterface(std::vector<bool> NeedInterface)
{
	ic.NeedInterface = NeedInterface;
}
void Parameter::setRe(real Re)
{
	ic.Re = Re;
}
void Parameter::setFactorPressBC(real factorPressBC)
{
	ic.factorPressBC = factorPressBC;
}
void Parameter::setIsGeo(bool isGeo)
{
	ic.isGeo = isGeo;
}
void Parameter::setIsGeoNormal(bool isGeoNormal)
{
	ic.isGeoNormal = isGeoNormal;
}
void Parameter::setIsInflowNormal(bool isInflowNormal)
{
	ic.isInflowNormal = isInflowNormal;
}
void Parameter::setIsOutflowNormal(bool isOutflowNormal)
{
	ic.isOutflowNormal = isOutflowNormal;
}
void Parameter::setIsProp(bool isProp)
{
	ic.isProp = isProp;
}
void Parameter::setIsCp(bool isCp)
{
	ic.isCp = isCp;
}
void Parameter::setConcFile(bool concFile)
{
	ic.isConc = concFile;
}
void Parameter::setStreetVelocityFile(bool streetVelocityFile)
{
	ic.streetVelocityFile = streetVelocityFile;
}
void Parameter::setUseMeasurePoints(bool useMeasurePoints)
{
	ic.isMeasurePoints = useMeasurePoints;
}
void Parameter::setUseWale(bool useWale)
{
	ic.isWale = useWale;
}
void Parameter::setUseInitNeq(bool useInitNeq)
{
	ic.isInitNeq = useInitNeq;
}
void Parameter::setSimulatePorousMedia(bool simulatePorousMedia)
{
	ic.simulatePorousMedia = simulatePorousMedia;
}

void Parameter::setIsF3(bool isF3)
{
	this->isF3 = isF3; 
}

void Parameter::setIsBodyForce(bool isBodyForce) 
{
	this->isBodyForce = isBodyForce;
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
void Parameter::setScaleLBMtoSI(std::vector<real> scaleLBMtoSI)
{
	ic.scaleLBMtoSI = scaleLBMtoSI;
}
void Parameter::setTranslateLBMtoSI(std::vector<real> translateLBMtoSI)
{
	ic.translateLBMtoSI = translateLBMtoSI;
}
void Parameter::setMinCoordX(std::vector<real> MinCoordX)
{
	ic.minCoordX = MinCoordX;
}
void Parameter::setMinCoordY(std::vector<real> MinCoordY)
{
	ic.minCoordY = MinCoordY;
}
void Parameter::setMinCoordZ(std::vector<real> MinCoordZ)
{
	ic.minCoordZ = MinCoordZ;
}
void Parameter::setMaxCoordX(std::vector<real> MaxCoordX)
{
	ic.maxCoordX = MaxCoordX;
}
void Parameter::setMaxCoordY(std::vector<real> MaxCoordY)
{
	ic.maxCoordY = MaxCoordY;
}
void Parameter::setMaxCoordZ(std::vector<real> MaxCoordZ)
{
	ic.maxCoordZ = MaxCoordZ;
}
void Parameter::setTempH(TempforBoundaryConditions* TempH)
{
	this->TempH = TempH;
}
void Parameter::setTempD(TempforBoundaryConditions* TempD)
{
	this->TempD = TempD;
}
void Parameter::setTempVelH(TempVelforBoundaryConditions* TempVelH)
{
	this->TempVelH = TempVelH;
}
void Parameter::setTempVelD(TempVelforBoundaryConditions* TempVelD)
{
	this->TempVelD = TempVelD;
}
void Parameter::setTempPressH(TempPressforBoundaryConditions* TempPressH)
{
	this->TempPressH = TempPressH;
}
void Parameter::setTempPressD(TempPressforBoundaryConditions* TempPressD)
{
	this->TempPressD = TempPressD;
}
//void Parameter::setkInflowQ(unsigned int kInflowQ)
//{
//   this->kInflowQ = kInflowQ;
//}
//void Parameter::setkOutflowQ(unsigned int kOutflowQ)
//{
//   this->kOutflowQ = kOutflowQ;
//}
//void Parameter::setQinflowH(QforBoundaryConditions* QinflowH)
//{
//   this->QinflowH = QinflowH;
//}
//void Parameter::setQinflowD(QforBoundaryConditions* QinflowD)
//{
//   this->QinflowD = QinflowD;
//}
//void Parameter::setQoutflowH(QforBoundaryConditions* QoutflowH)
//{
//   this->QoutflowH = QoutflowH;
//}
//void Parameter::setQoutflowD(QforBoundaryConditions* QoutflowD)
//{
//   this->QoutflowD = QoutflowD;
//}
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
void Parameter::setscaleCFC(std::string scaleCFC)
{
	ic.scaleCFC = scaleCFC;
}
void Parameter::setscaleCFF(std::string scaleCFF)
{
	ic.scaleCFF = scaleCFF;
}
void Parameter::setscaleFCC(std::string scaleFCC)
{
	ic.scaleFCC = scaleFCC;
}
void Parameter::setscaleFCF(std::string scaleFCF)
{
	ic.scaleFCF = scaleFCF;
}
void Parameter::setscaleOffsetCF(std::string scaleOffsetCF)
{
	ic.scaleOffsetCF = scaleOffsetCF;
}
void Parameter::setscaleOffsetFC(std::string scaleOffsetFC)
{
	ic.scaleOffsetFC = scaleOffsetFC;
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
void Parameter::setslipBcPos(std::string slipBcPos)
{
	ic.slipBcPos = slipBcPos;
}
void Parameter::setslipBcQs(std::string slipBcQs)
{
	ic.slipBcQs = slipBcQs;
}
void Parameter::setslipBcValue(std::string slipBcValue)
{
	ic.slipBcValue = slipBcValue;
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
void Parameter::setpropellerQs(std::string propellerQs)
{
	ic.propellerQs = propellerQs;
}
void Parameter::setpropellerValues(std::string propellerValues)
{
	ic.propellerValues = propellerValues;
}
void Parameter::setpropellerCylinder(std::string propellerCylinder)
{
	ic.propellerCylinder = propellerCylinder;
}
void Parameter::setmeasurePoints(std::string measurePoints)
{
	ic.measurePoints = measurePoints;
}
void Parameter::setnumberNodes(std::string numberNodes)
{
	ic.numberNodes = numberNodes;
}
void Parameter::setLBMvsSI(std::string LBMvsSI)
{
	ic.LBMvsSI = LBMvsSI;
}
void Parameter::setcpTop(std::string cpTop)
{
	ic.cpTop = cpTop;
}
void Parameter::setcpBottom(std::string cpBottom)
{
	ic.cpBottom = cpBottom;
}
void Parameter::setcpBottom2(std::string cpBottom2)
{
	ic.cpBottom2 = cpBottom2;
}
void Parameter::setConcentration(std::string concFile)
{
	ic.concentration = concFile;
}
void Parameter::setStreetVelocity(std::string streetVelocity)
{
	ic.streetVelocity = streetVelocity;
}
void Parameter::setclockCycleForMP(real clockCycleForMP)
{
	ic.clockCycleForMP = clockCycleForMP;
}
void Parameter::setTimeDoCheckPoint(unsigned int tDoCheckPoint)
{
	ic.tDoCheckPoint = tDoCheckPoint;
}
void Parameter::setTimeDoRestart(unsigned int tDoRestart)
{
	ic.tDoRestart = tDoRestart;
}
void Parameter::setDoCheckPoint(bool doCheckPoint)
{
	ic.doCheckPoint = doCheckPoint;
}
void Parameter::setDoRestart(bool doRestart)
{
	ic.doRestart = doRestart;
}
void Parameter::settimestepForMP(unsigned int timestepForMP)
{
	ic.timeStepForMP = timestepForMP;
}
void Parameter::setObj(std::string str, bool isObj)
{
	if (str == "geo")
	{
		this->setIsGeo(isObj);
	}
	else if (str == "prop")
	{
		this->setIsProp(isObj);
	}
	else if (str == "cp")
	{
		this->setIsCp(isObj);
	}
	else if (str == "geoNormal")
	{
		this->setIsGeoNormal(isObj);
	}
	else if (str == "inflowNormal")
	{
		this->setIsInflowNormal(isObj);
	}
	else if (str == "outflowNormal")
	{
		this->setIsOutflowNormal(isObj);
	}
}
void Parameter::setGeometryValues(bool GeometryValues)
{
	ic.GeometryValues = GeometryValues;
}
void Parameter::setCalc2ndOrderMoments(bool is2ndOrderMoments)
{
	ic.is2ndOrderMoments = is2ndOrderMoments;
}
void Parameter::setCalc3rdOrderMoments(bool is3rdOrderMoments)
{
	ic.is3rdOrderMoments = is3rdOrderMoments;
}
void Parameter::setCalcHighOrderMoments(bool isHighOrderMoments)
{
	ic.isHighOrderMoments = isHighOrderMoments;
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
//1D domain decomposition
void Parameter::setPossNeighborFiles(std::vector<std::string> possNeighborFiles, std::string sor)
{
	if (sor=="send")
	{
		this->possNeighborFilesSend = possNeighborFiles;
	} 
	else if (sor == "recv")
	{
		this->possNeighborFilesRecv = possNeighborFiles;
	}
}
void Parameter::setNumberOfProcessNeighbors(unsigned int numberOfProcessNeighbors, int level, std::string sor)
{
	if (sor=="send")
	{
		parH[level]->sendProcessNeighbor.resize(numberOfProcessNeighbors);
		parD[level]->sendProcessNeighbor.resize(numberOfProcessNeighbors);
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighbor.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighbor.resize(numberOfProcessNeighbors);
	}
}
void Parameter::setIsNeighbor(bool isNeigbor)
{
	this->isNeigbor = isNeigbor;
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
		if (getDiffOn()==true){
			parH[level]->sendProcessNeighborADX.resize(numberOfProcessNeighbors);
			parD[level]->sendProcessNeighborADX.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborX.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->recvProcessNeighborADX.resize(numberOfProcessNeighbors);
			parD[level]->recvProcessNeighborADX.resize(numberOfProcessNeighbors);
		}
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
		if (getDiffOn()==true){
			parH[level]->sendProcessNeighborADY.resize(numberOfProcessNeighbors);
			parD[level]->sendProcessNeighborADY.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborY.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
			parH[level]->recvProcessNeighborADY.resize(numberOfProcessNeighbors);
			parD[level]->recvProcessNeighborADY.resize(numberOfProcessNeighbors);
		}
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
		if (getDiffOn()==true){
			parH[level]->sendProcessNeighborADZ.resize(numberOfProcessNeighbors);
			parD[level]->sendProcessNeighborADZ.resize(numberOfProcessNeighbors);
		}
		//////////////////////////////////////////////////////////////////////////
	} 
	else if (sor == "recv")
	{
		parH[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
		parD[level]->recvProcessNeighborZ.resize(numberOfProcessNeighbors);
		//////////////////////////////////////////////////////////////////////////
		if (getDiffOn()==true){
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
void Parameter::setgeomBoundaryNormalX(std::string geomNormalX)
{
	ic.geomNormalX = geomNormalX;
}
void Parameter::setgeomBoundaryNormalY(std::string geomNormalY)
{
	ic.geomNormalY = geomNormalY;
}
void Parameter::setgeomBoundaryNormalZ(std::string geomNormalZ)
{
	ic.geomNormalZ = geomNormalZ;
}
void Parameter::setInflowBoundaryNormalX(std::string inflowNormalX)
{
	ic.inflowNormalX = inflowNormalX;
}
void Parameter::setInflowBoundaryNormalY(std::string inflowNormalY)
{
	ic.inflowNormalY = inflowNormalY;
}
void Parameter::setInflowBoundaryNormalZ(std::string inflowNormalZ)
{
	ic.inflowNormalZ = inflowNormalZ;
}
void Parameter::setOutflowBoundaryNormalX(std::string outflowNormalX)
{
	ic.outflowNormalX = outflowNormalX;
}
void Parameter::setOutflowBoundaryNormalY(std::string outflowNormalY)
{
	ic.outflowNormalY = outflowNormalY;
}
void Parameter::setOutflowBoundaryNormalZ(std::string outflowNormalZ)
{
	ic.outflowNormalZ = outflowNormalZ;
}
void Parameter::setMainKernel(std::string kernel)
{
	this->mainKernel = kernel;
}
void Parameter::setMultiKernelOn(bool isOn)
{
	this->multiKernelOn = isOn;
}
void Parameter::setMultiKernelLevel(std::vector< int> kernelLevel)
{
	this->multiKernelLevel = kernelLevel;
}
void Parameter::setMultiKernel(std::vector< std::string> kernel)
{
	this->multiKernel = kernel;
}
void Parameter::setADKernel(std::string adKernel)
{
	this->adKernel = adKernel;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double* Parameter::getForcesDouble()
{
	return this->hostForcing;
}
real* Parameter::getForcesHost()
{
	return this->forcingH;
}
real* Parameter::getForcesDev()
{
	return this->forcingD;
}
double * Parameter::getQuadricLimitersDouble()
{
    return this->hostQuadricLimiters;
}
real * Parameter::getQuadricLimitersHost()
{
    return this->quadricLimitersH;
}
real * Parameter::getQuadricLimitersDev()
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
real Parameter::getStartXHotWall()
{
	return this->startXHotWall;
}
real Parameter::getEndXHotWall()
{
	return this->endXHotWall;
}
unsigned int Parameter::getStepEnsight()
{
	return this->stepEnsight;
}
unsigned int Parameter::getOutputCount()
{
	return this->outputCount;
}
unsigned int Parameter::getlimitOfNodesForVTK()
{
	return this->limitOfNodesForVTK;
}
unsigned int Parameter::getStartTurn()
{
	return startTurn;
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
unsigned int Parameter::getMemSizereal(int level)
{
	return parH[level]->mem_size_real;
}
unsigned int Parameter::getMemSizeInt(int level)
{
	return parH[level]->mem_size_int;
}
unsigned int Parameter::getMemSizeBool(int level)
{
	return parH[level]->mem_size_bool;
}
unsigned int Parameter::getMemSizerealYZ(int level)
{
	return parH[level]->mem_size_real_yz;
}
int Parameter::getFine()
{
	return fine;
}
int Parameter::getCoarse()
{
	return coarse;
}
int Parameter::getParticleBasicLevel()
{
	return this->particleBasicLevel;
}
int Parameter::getParticleInitLevel()
{
	return this->particleInitLevel;
}
int Parameter::getNumberOfParticles()
{
	return this->numberOfParticles;
}
bool Parameter::getEvenOrOdd(int level)
{
	return parH[level]->evenOrOdd;
}
bool Parameter::getDiffOn()
{
	return diffOn;
}
bool Parameter::getCompOn()
{
	return compOn;
}
int Parameter::getDiffMod()
{
	return diffMod;
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
	if (getDoRestart())
	{
		return getTimeDoRestart() + 1;
	} 
	else
	{
		return 1;
	}
}
unsigned int Parameter::getTInit()
{
	if (getDoRestart())
	{
		return getTimeDoRestart();
	} 
	else
	{
		return 0;
	}
}
unsigned int Parameter::getTEnd()
{
	return ic.tend;
}
unsigned int Parameter::getTOut()
{
	return ic.tout;
}
unsigned int Parameter::getTStartOut()
{
	return ic.tStartOut;
}
bool Parameter::getCalcMedian()
{
	return ic.calcMedian;
}
bool Parameter::getCalcDragLift()
{
	return this->calcDragLift;
}
bool Parameter::getCalcCp()
{
	return this->calcCp;
}
bool Parameter::getCalcParticle()
{
	return this->calcParticles;
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
	return ic.tCalcMedStart;
}
int Parameter::getTimeCalcMedEnd()
{
	return ic.tCalcMedEnd;
}
std::string Parameter::getOutputPath()
{
	return ic.oPath;
}
std::string Parameter::getOutputPrefix()
{
	return ic.oPrefix;
}
std::string Parameter::getFName()
{
	return ic.fname;
}
bool Parameter::getPrintFiles()
{
	return ic.printFiles;
}
bool Parameter::getReadGeo()
{
	return ic.readGeo;
}
real Parameter::getDiffusivity()
{
	return ic.Diffusivity;
}
real Parameter::getTemperatureInit()
{
	return ic.Temp;
}
real Parameter::getTemperatureBC()
{
	return ic.TempBC;
}
real Parameter::getViscosity()
{
	return ic.vis;
}
real Parameter::getVelocity()
{
	return ic.u0;
}
real Parameter::getViscosityRatio()
{
	return ic.vis_ratio;
}
real Parameter::getVelocityRatio()
{
	return ic.u0_ratio;
}
real Parameter::getDensityRatio()
{
	return ic.delta_rho;
}
real Parameter::getPressRatio()
{
	return ic.delta_press;
}
real Parameter::getRealX()
{
	return ic.RealX;
}
real Parameter::getRealY()
{
	return ic.RealY;
}
unsigned int Parameter::getPressInID()
{
	return ic.PressInID;
}
unsigned int Parameter::getPressOutID()
{
	return ic.PressOutID;
}
unsigned int Parameter::getPressInZ()
{
	return ic.PressInZ;
}
unsigned int Parameter::getPressOutZ()
{
	return ic.PressOutZ;
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
std::vector<uint> Parameter::getDevices()
{
	return ic.devices;
}
std::string Parameter::getGeometryFileC()
{
	return ic.geometryFileC;
}
std::string Parameter::getGeometryFileM()
{
	return ic.geometryFileM;
}
std::string Parameter::getGeometryFileF()
{
	return ic.geometryFileF;
}
std::vector<bool> Parameter::getNeedInterface()
{
	return ic.NeedInterface;
}
real Parameter::getRe()
{
	return ic.Re;
}
real Parameter::getFactorPressBC()
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
std::vector<real> Parameter::getScaleLBMtoSI()
{
	return ic.scaleLBMtoSI;
}
std::vector<real> Parameter::getTranslateLBMtoSI()
{
	return ic.translateLBMtoSI;
}
std::vector<real> Parameter::getMinCoordX()
{
	return ic.minCoordX;
}
std::vector<real> Parameter::getMinCoordY()
{
	return ic.minCoordY;
}
std::vector<real> Parameter::getMinCoordZ()
{
	return ic.minCoordZ;
}
std::vector<real> Parameter::getMaxCoordX()
{
	return ic.maxCoordX;
}
std::vector<real> Parameter::getMaxCoordY()
{
	return ic.maxCoordY;
}
std::vector<real> Parameter::getMaxCoordZ()
{
	return ic.maxCoordZ;
}
TempforBoundaryConditions* Parameter::getTempH()
{
	return this->TempH;
}
TempforBoundaryConditions* Parameter::getTempD()
{
	return this->TempD;
}
TempVelforBoundaryConditions* Parameter::getTempVelH()
{
	return this->TempVelH;
}
TempVelforBoundaryConditions* Parameter::getTempVelD()
{
	return this->TempVelD;
}
TempPressforBoundaryConditions* Parameter::getTempPressH()
{
	return this->TempPressH;
}
TempPressforBoundaryConditions* Parameter::getTempPressD()
{
	return this->TempPressD;
}
//unsigned int Parameter::getkInflowQ()
//{
//   return this->kInflowQ;
//}
//unsigned int Parameter::getkOutflowQ()
//{
//   return this->kOutflowQ;
//}
//QforBoundaryConditions* Parameter::getQinflowH()
//{
//   return this->QinflowH;
//}
//QforBoundaryConditions* Parameter::getQinflowD()
//{
//   return this->QinflowD;
//}
//QforBoundaryConditions* Parameter::getQoutflowH()
//{
//   return this->QoutflowH;
//}
//QforBoundaryConditions* Parameter::getQoutflowD()
//{
//   return this->QoutflowD;
//}
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
std::string Parameter::getscaleCFC()
{
	return ic.scaleCFC;
}
std::string Parameter::getscaleCFF()
{
	return ic.scaleCFF;
}
std::string Parameter::getscaleFCC()
{
	return ic.scaleFCC;
}
std::string Parameter::getscaleFCF()
{
	return ic.scaleFCF;
}
std::string Parameter::getscaleOffsetCF()
{
	return ic.scaleOffsetCF;
}
std::string Parameter::getscaleOffsetFC()
{
	return ic.scaleOffsetFC;
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
std::string Parameter::getslipBcPos()
{
	return ic.slipBcPos;
}
std::string Parameter::getslipBcQs()
{
	return ic.slipBcQs;
}
std::string Parameter::getslipBcValue()
{
	return ic.slipBcValue;
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
std::string Parameter::getpropellerQs()
{
	return ic.propellerQs;
}
std::string Parameter::getpropellerValues()
{
	return ic.propellerValues;
}
std::string Parameter::getpropellerCylinder()
{
	return ic.propellerCylinder;
}
std::string Parameter::getmeasurePoints()
{
	return ic.measurePoints;
}
std::string Parameter::getLBMvsSI()
{
	return ic.LBMvsSI;
}
std::string Parameter::getnumberNodes()
{
	return ic.numberNodes;
}
std::string Parameter::getcpTop()
{
	return ic.cpTop;
}
std::string Parameter::getcpBottom()
{
	return ic.cpBottom;
}
std::string Parameter::getcpBottom2()
{
	return ic.cpBottom2;
}
std::string Parameter::getConcentration()
{
	return ic.concentration;
}
std::string Parameter::getStreetVelocityFilePath()
{
	return ic.streetVelocity;
}
real Parameter::getclockCycleForMP()
{
	return ic.clockCycleForMP;
}
unsigned int Parameter::getTimeDoCheckPoint()
{
	return ic.tDoCheckPoint;
}
unsigned int Parameter::getTimeDoRestart()
{
	return ic.tDoRestart;
}
bool Parameter::getDoCheckPoint()
{
	return ic.doCheckPoint;
}
bool Parameter::getDoRestart()
{
	return ic.doRestart;
}
bool Parameter::getIsGeo()
{
	return ic.isGeo;
}
bool Parameter::getIsGeoNormal()
{
	return ic.isGeoNormal;
}
bool Parameter::getIsInflowNormal()
{
	return ic.isInflowNormal;
}
bool Parameter::getIsOutflowNormal()
{
	return ic.isOutflowNormal;
}
bool Parameter::getIsCp()
{
	return ic.isCp;
}
bool Parameter::getConcFile()
{
	return ic.isConc;
}
bool Parameter::isStreetVelocityFile()
{
	return ic.streetVelocityFile;
}
bool Parameter::getUseMeasurePoints()
{
	return ic.isMeasurePoints;
}
bool Parameter::getUseWale()
{
	return ic.isWale;
}
bool Parameter::getUseInitNeq()
{
	return ic.isInitNeq;
}
bool Parameter::getSimulatePorousMedia()
{
	return ic.simulatePorousMedia;
}

bool Parameter::getIsF3()
{
	return this->isF3; 
}

bool Parameter::getIsBodyForce() 
{ 
	return this->isBodyForce; 
}

bool Parameter::getIsGeometryValues()
{
	return ic.GeometryValues;
}
bool Parameter::getCalc2ndOrderMoments()
{
	return ic.is2ndOrderMoments;
}
bool Parameter::getCalc3rdOrderMoments()
{
	return ic.is3rdOrderMoments;
}
bool Parameter::getCalcHighOrderMoments()
{
	return ic.isHighOrderMoments;
}
bool Parameter::getIsProp()
{
	return ic.isProp;
}
bool Parameter::overWritingRestart(unsigned int t)
{
	if (t == getTimeDoRestart())
	{
		return true;
	} 
	else
	{
		return false;
	}
}
unsigned int Parameter::getTimestepForMP()
{
	return ic.timeStepForMP;
}
unsigned int Parameter::getTimestepOfCoarseLevel()
{
	return this->timestep;
}
double Parameter::getMemsizeGPU()
{
	return this->memsizeGPU;
}
//1D domain decomposition
std::vector<std::string> Parameter::getPossNeighborFiles(std::string sor)
{
	if (sor=="send")
	{
		return this->possNeighborFilesSend;
	} 
	else if (sor == "recv")
	{
		return this->possNeighborFilesRecv;
	}
    throw std::runtime_error("Parameter string invalid.");
}
unsigned int Parameter::getNumberOfProcessNeighbors(int level, std::string sor)
{
	if (sor=="send")
	{
		return (unsigned int)parH[level]->sendProcessNeighbor.size();
	} 
	else if (sor == "recv")
	{
		return (unsigned int)parH[level]->recvProcessNeighbor.size();
	}
    throw std::runtime_error("Parameter string invalid.");
}
bool Parameter::getIsNeighbor()
{
	return this->isNeigbor;
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
    throw std::runtime_error("Parameter string invalid.");
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
    throw std::runtime_error("Parameter string invalid.");
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
    throw std::runtime_error("Parameter string invalid.");
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
    throw std::runtime_error("Parameter string invalid.");
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
    throw std::runtime_error("Parameter string invalid.");
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
    throw std::runtime_error("Parameter string invalid.");
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
	return ic.geomNormalX;
}
std::string Parameter::getgeomBoundaryNormalY()
{
	return ic.geomNormalY;
}
std::string Parameter::getgeomBoundaryNormalZ()
{
	return ic.geomNormalZ;
}
std::string Parameter::getInflowBoundaryNormalX()
{
	return ic.inflowNormalX;
}
std::string Parameter::getInflowBoundaryNormalY()
{
	return ic.inflowNormalY;
}
std::string Parameter::getInflowBoundaryNormalZ()
{
	return ic.inflowNormalZ;
}
std::string Parameter::getOutflowBoundaryNormalX()
{
	return ic.outflowNormalX;
}
std::string Parameter::getOutflowBoundaryNormalY()
{
	return ic.outflowNormalY;
}
std::string Parameter::getOutflowBoundaryNormalZ()
{
	return ic.outflowNormalZ;
}
curandState* Parameter::getRandomState()
{
	return this->devState;
}

std::string Parameter::getMainKernel()
{
	return mainKernel;
}
bool Parameter::getMultiKernelOn()
{
	return multiKernelOn;
}
std::vector< int> Parameter::getMultiKernelLevel()
{
	return multiKernelLevel;
}
std::vector< std::string> Parameter::getMultiKernel()
{
	return multiKernel;
}
std::string Parameter::getADKernel()
{
	return adKernel;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Parameter::setInitialCondition(std::function<void(real,real,real,real&,real&,real&,real&)> initialCondition)
{
    this->initialCondition = initialCondition;
}

std::function<void(real,real,real,real&,real&,real&,real&)>& Parameter::getInitialCondition()
{
    return this->initialCondition;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//private methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::initInterfaceParameter(int level)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//host
	parH[level]->K_CF               = 0;
	parH[level]->K_FC               = 0;
	if (parH[level]->need_interface[INTERFACE_E]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNZ   /2);
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNZ-6)/2);
	} 
	if (parH[level]->need_interface[INTERFACE_W]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNZ   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNZ-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_N]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNX   /2)*( parH[level+1]->gridNZ   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNX-6)/2)*((parH[level+1]->gridNZ-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_S]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNX   /2)*( parH[level+1]->gridNZ   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNX-6)/2)*((parH[level+1]->gridNZ-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_T]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNX   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNX-6)/2);
	}
	if (parH[level]->need_interface[INTERFACE_B]==true)
	{
		parH[level]->K_CF           += ( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNX   /2); 
		parH[level]->K_FC           += ((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNX-6)/2);
	}
	//parH[level]->K_CF               = (( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNZ   /2)*2)+
	//                                  (( parH[level+1]->gridNX   /2)*( parH[level+1]->gridNZ   /2)*2)+
	//                                  (( parH[level+1]->gridNY   /2)*( parH[level+1]->gridNX   /2)*2);
	//parH[level]->K_FC               = (((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNZ-6)/2)*2)+
	//                                  (((parH[level+1]->gridNX-6)/2)*((parH[level+1]->gridNZ-6)/2)*2)+
	//                                  (((parH[level+1]->gridNY-6)/2)*((parH[level+1]->gridNX-6)/2)*2);
	parH[level]->mem_size_kCF       = sizeof(unsigned int)*parH[level]->K_CF;
	parH[level]->mem_size_kFC       = sizeof(unsigned int)*parH[level]->K_FC;
	parH[level]->mem_size_kCF_off   = sizeof(real)*parH[level]->K_CF;
	parH[level]->mem_size_kFC_off   = sizeof(real)*parH[level]->K_FC;
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//device
	parD[level]->K_CF               = parH[level]->K_CF;
	parD[level]->K_FC               = parH[level]->K_FC;
	parD[level]->mem_size_kCF       = parH[level]->mem_size_kCF;
	parD[level]->mem_size_kFC       = parH[level]->mem_size_kFC;
	parD[level]->mem_size_kCF_off   = parH[level]->mem_size_kCF_off;
	parD[level]->mem_size_kFC_off   = parH[level]->mem_size_kFC_off;
	///////////////////////////////////////////////////////////////////////////////////////////////////
}
real Parameter::TrafoXtoWorld(int CoordX, int level)
{
	return (parH[level]->mTtoWx*CoordX+parH[level]->cTtoWx);
}
real Parameter::TrafoYtoWorld(int CoordY, int level)
{
	return (parH[level]->mTtoWy*CoordY+parH[level]->cTtoWy);
}
real Parameter::TrafoZtoWorld(int CoordZ, int level)
{
	return (parH[level]->mTtoWz*CoordZ+parH[level]->cTtoWz);
}
real Parameter::TrafoXtoMGsWorld(int CoordX, int level)
{
	real temp = 0;
	for (int i = 0; i <= level; i++)
	{
		temp += (parH[i]->XdistKn + 0.25f) * 2.f * parH[i]->dx;
	}
	temp += (real)((CoordX ) * parH[level]->dx);
	return temp;
}
real Parameter::TrafoYtoMGsWorld(int CoordY, int level)
{
	real temp = 0;
	for (int i = 0; i <= level; i++)
	{
		temp += (parH[i]->YdistKn + 0.25f) * 2.f * parH[i]->dx;
	}
	temp += (real)((CoordY ) * parH[level]->dx);
	return temp;
}
real Parameter::TrafoZtoMGsWorld(int CoordZ, int level)
{
	real temp = 0;
	for (int i = 0; i <= level; i++)
	{
		temp += (parH[i]->ZdistKn + 0.25f) * 2.f * parH[i]->dx;
	}
	temp += (real)((CoordZ) * parH[level]->dx);
	return temp;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

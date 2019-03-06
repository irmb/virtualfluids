#include "ConfigFileReader.h"
#include "../ConfigData/ConfigDataImp.h"
#include "../Input.h"
#include "../../StringUtilities/StringUtil.h"

#include <iostream>
#include <fstream>


VF_PUBLIC std::shared_ptr<ConfigFileReader> ConfigFileReader::getNewInstance()
{
	return std::shared_ptr<ConfigFileReader>(new ConfigFileReader());
}

ConfigFileReader::ConfigFileReader()
{

}

VF_PUBLIC ConfigFileReader::~ConfigFileReader()
{

}

VF_PUBLIC std::shared_ptr<ConfigData> ConfigFileReader::readConfigFile(const std::string &filePath) const
{
	std::shared_ptr<ConfigDataImp> data = ConfigDataImp::getNewInstance();
	std::ifstream stream;
	stream.open(filePath.c_str(), std::ios::in);
	if (stream.fail())
		throw std::runtime_error("can not open config file!"); 
	std::unique_ptr<input::Input> input = input::Input::makeInput(stream, "config");

	if (input->getValue("NumberOfDevices") != "")
		data->setNumberOfDevices(StringUtil::toInt(input->getValue("NumberOfDevices")));

	if (input->getValue("Devices") != "")
		data->setDevices(StringUtil::toUintVector(input->getValue("Devices")));

	if (input->getValue("Path") != "")
		data->setOutputPath(input->getValue("Path"));

	if (input->getValue("Prefix") != "")
		data->setPrefix(input->getValue("Prefix"));

	if (input->getValue("GridPath") != "")
		data->setGridPath(input->getValue("GridPath"));
	else
	{
		std::cout << "GridPath has to be defined!" << std::endl;
		exit(1);
	}


	if (input->getValue("WriteGrid") != "")
		data->setPrintOutputFiles(StringUtil::toBool(input->getValue("WriteGrid")));

	if (input->getValue("GeometryValues") != "")
		data->setGeometryValues(StringUtil::toBool(input->getValue("GeometryValues")));

	if (input->getValue("calc2ndOrderMoments") != "")
		data->setCalc2ndOrderMoments(StringUtil::toBool(input->getValue("calc2ndOrderMoments")));

	if (input->getValue("calc3rdOrderMoments") != "")
		data->setCalc3rdOrderMoments(StringUtil::toBool(input->getValue("calc3rdOrderMoments")));

	if (input->getValue("calcHigherOrderMoments") != "")
		data->setCalcHighOrderMoments(StringUtil::toBool(input->getValue("calcHigherOrderMoments")));

	if (input->getValue("ReadGeometry") != "")
		data->setReadGeo(StringUtil::toBool(input->getValue("ReadGeometry")));

	if (input->getValue("calcMedian") != "")
		data->setCalcMedian(StringUtil::toBool(input->getValue("calcMedian")));

	if (input->getValue("UseConcFile") != "")
		data->setConcFile(StringUtil::toBool(input->getValue("UseConcFile")));

	if (input->getValue("UseMeasurePoints") != "")
		data->setUseMeasurePoints(StringUtil::toBool(input->getValue("UseMeasurePoints")));

	if (input->getValue("UseWale") != "")
		data->setUseWale(StringUtil::toBool(input->getValue("UseWale")));

	if (input->getValue("SimulatePorousMedia") != "")
		data->setSimulatePorousMedia(StringUtil::toBool(input->getValue("SimulatePorousMedia")));

	if (input->getValue("D3Qxx") != "")
		data->setD3Qxx(StringUtil::toInt(input->getValue("D3Qxx")));

	if (input->getValue("TimeEnd") != "")
		data->setTEnd(StringUtil::toInt(input->getValue("TimeEnd")));

	if (input->getValue("TimeOut") != "")
		data->setTOut(StringUtil::toInt(input->getValue("TimeOut")));
	
	if (input->getValue("TimeStartOut") != "")
		data->setTStartOut(StringUtil::toInt(input->getValue("TimeStartOut")));

	if (input->getValue("TimeStartCalcMedian") != "")
		data->setTimeCalcMedStart(StringUtil::toInt(input->getValue("TimeStartCalcMedian")));

	if (input->getValue("TimeEndCalcMedian") != "")
		data->setTimeCalcMedEnd(StringUtil::toInt(input->getValue("TimeEndCalcMedian")));

	if (input->getValue("PressInID") != "")
		data->setPressInID(StringUtil::toInt(input->getValue("PressInID")));

	if (input->getValue("PressOutID") != "")
		data->setPressOutID(StringUtil::toInt(input->getValue("PressOutID")));

	if (input->getValue("PressInZ") != "")
		data->setPressInZ(StringUtil::toInt(input->getValue("PressInZ")));

	if (input->getValue("PressOutZ") != "")
		data->setPressOutZ(StringUtil::toInt(input->getValue("PressOutZ")));
	//////////////////////////////////////////////////////////////////////////
	if (input->getValue("DiffOn") != "")
		data->setDiffOn(StringUtil::toBool(input->getValue("DiffOn")));

	if (input->getValue("DiffMod") != "")
		data->setDiffMod(StringUtil::toInt(input->getValue("DiffMod")));

	if (input->getValue("Diffusivity") != "")
		data->setDiffusivity(StringUtil::toFloat(input->getValue("Diffusivity")));

	if (input->getValue("Temp") != "")
		data->setTemperatureInit(StringUtil::toFloat(input->getValue("Temp")));

	if (input->getValue("TempBC") != "")
		data->setTemperatureBC(StringUtil::toFloat(input->getValue("TempBC")));
	//////////////////////////////////////////////////////////////////////////
	if (input->getValue("Viscosity_LB") != "")
		data->setViscosity(StringUtil::toFloat(input->getValue("Viscosity_LB")));

	if (input->getValue("Velocity_LB") != "")
		data->setVelocity(StringUtil::toFloat(input->getValue("Velocity_LB")));

	if (input->getValue("Viscosity_Ratio_World_to_LB") != "")
		data->setViscosityRatio(StringUtil::toFloat(input->getValue("Viscosity_Ratio_World_to_LB")));

	if (input->getValue("Velocity_Ratio_World_to_LB") != "")
		data->setVelocityRatio(StringUtil::toFloat(input->getValue("Velocity_Ratio_World_to_LB")));

	if (input->getValue("Density_Ratio_World_to_LB") != "")
		data->setDensityRatio(StringUtil::toFloat(input->getValue("Density_Ratio_World_to_LB")));

	if (input->getValue("Delta_Press") != "")
		data->setPressRatio(StringUtil::toFloat(input->getValue("Delta_Press")));

	if (input->getValue("SliceRealX") != "")
		data->setRealX(StringUtil::toFloat(input->getValue("SliceRealX")));

	if (input->getValue("SliceRealY") != "")
		data->setRealY(StringUtil::toFloat(input->getValue("SliceRealY")));

	if (input->getValue("FactorPressBC") != "")
		data->setFactorPressBC(StringUtil::toFloat(input->getValue("FactorPressBC")));

	if (input->getValue("GeometryC") != "")
		data->setGeometryFileC(input->getValue("GeometryC"));

	if (input->getValue("GeometryM") != "")
		data->setGeometryFileM(input->getValue("GeometryM"));

	if (input->getValue("GeometryF") != "")
		data->setGeometryFileF(input->getValue("GeometryF"));
	//////////////////////////////////////////////////////////////////////////
	if (input->getValue("measureClockCycle") != "")
		data->setClockCycleForMP(StringUtil::toFloat(input->getValue("measureClockCycle")));

	if (input->getValue("measureTimestep") != "")
		data->setTimestepForMP(StringUtil::toInt(input->getValue("measureTimestep")));
	//////////////////////////////////////////////////////////////////////////
	//Forcing
	if (input->getValue("ForcingX") != "")
		data->setForcingX(StringUtil::toFloat(input->getValue("ForcingX")));
	if (input->getValue("ForcingY") != "")
		data->setForcingY(StringUtil::toFloat(input->getValue("ForcingY")));
	if (input->getValue("ForcingZ") != "")
		data->setForcingZ(StringUtil::toFloat(input->getValue("ForcingZ")));
	//////////////////////////////////////////////////////////////////////////
	//Particles
	if (input->getValue("calcParticles") != "")
		data->setCalcParticles(StringUtil::toBool(input->getValue("calcParticles")));

	if (input->getValue("baseLevel") != "")
		data->setParticleBasicLevel(StringUtil::toInt(input->getValue("baseLevel")));

	if (input->getValue("initLevel") != "")
		data->setParticleInitLevel(StringUtil::toInt(input->getValue("initLevel")));

	if (input->getValue("numberOfParticles") != "")
		data->setNumberOfParticles(StringUtil::toInt(input->getValue("numberOfParticles")));

	if (input->getValue("startXHotWall") != "")
		data->setStartXHotWall(StringUtil::toDouble(input->getValue("startXHotWall")));

	if (input->getValue("endXHotWall") != "")
		data->setEndXHotWall(StringUtil::toDouble(input->getValue("endXHotWall")));
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Restart
	if (input->getValue("TimeDoCheckPoint") != "")
		data->setTimeDoCheckPoint(StringUtil::toInt(input->getValue("TimeDoCheckPoint")));

	if (input->getValue("TimeDoRestart") != "")
		data->setTimeDoRestart(StringUtil::toInt(input->getValue("TimeDoRestart")));

	if (input->getValue("DoCheckPoint") != "")
		data->setDoCheckPoint(StringUtil::toBool(input->getValue("DoCheckPoint")));

	if (input->getValue("DoRestart") != "")
		data->setDoRestart(StringUtil::toBool(input->getValue("DoRestart")));
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (input->getValue("NOGL") != "")
		data->setMaxLevel(StringUtil::toInt(input->getValue("NOGL")));

	if (input->getValue("GridX") != "")
		data->setGridX(StringUtil::toIntVector(input->getValue("GridX")));

	if (input->getValue("GridY") != "")
		data->setGridY(StringUtil::toIntVector(input->getValue("GridY")));

	if (input->getValue("GridZ") != "")
		data->setGridZ(StringUtil::toIntVector(input->getValue("GridZ")));

	if (input->getValue("DistX") != "")
		data->setDistX(StringUtil::toIntVector(input->getValue("DistX")));

	if (input->getValue("DistY") != "")
		data->setDistY(StringUtil::toIntVector(input->getValue("DistY")));

	if (input->getValue("DistZ") != "")
		data->setDistZ(StringUtil::toIntVector(input->getValue("DistZ")));

	if (input->getValue("NeedInterface") != "")
		data->setNeedInterface(StringUtil::toBoolVector(input->getValue("NeedInterface")));
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Kernel
	if (input->getValue("MainKernelName") != "")
		data->setMainKernel(input->getValue("MainKernelName"));

	if (input->getValue("MultiKernelOn") != "")
		data->setMultiKernelOn(StringUtil::toBool(input->getValue("MultiKernelOn")));

	if (input->getValue("MultiKernelLevel") != "")
		data->setMultiKernelLevel(StringUtil::toIntVector(input->getValue("MultiKernelLevel")));

	if (input->getValue("MultiKernelName") != "")
		data->setMultiKernelName(StringUtil::toStringVector(input->getValue("MultiKernelName")));

	if (StringUtil::toStringVector(input->getValue("MultiKernelName")).size() != StringUtil::toIntVector(input->getValue("MultiKernelLevel")).size())
	{
		std::cout << "MultiKernelName and MultiKernelLevel has to be of same size!" << std::endl;
		exit(1);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return data;

}





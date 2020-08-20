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
//! \file Parameter.cpp
//! \ingroup Parameter
//! \author Martin Schoenherr
//=======================================================================================
#include "Parameter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

SPtr<Parameter> Parameter::make()
{
    return SPtr<Parameter>(new Parameter());
}

Parameter::Parameter()
{
	this->setOutputPath("C:/Output/");

	this->setOutputPrefix("MyFile");

	this->setPrintFiles(false);

	this->setD3Qxx((int)27);

	this->setTimestepEnd((uint)10);

	this->setTimestepOut((uint)1);

	this->setTimestepStartOut((uint)0);

	this->setViscosityLB((real)0.001);

	this->setVelocityLB((real)0.01);

	this->setViscosityRatio((real)1.0);

	this->setVelocityRatio((real)1.0);

	this->setDensityRatio((real)1.0);

	this->setPressureRatio((real)1.0);

	this->setPathAndFilename(this->getOutputPath() + "/" + this->getOutputPrefix());

	this->setlimitOfNodesForVTK((uint)30000000);
}
Parameter::~Parameter()
{
}
Parameter* Parameter::instance = 0;
Parameter* Parameter::getInstanz()
{
	if( instance == 0 )
		instance = new Parameter();
	return instance;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//init-method
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::initParameter()
{
	//host
	parametersOnHost                          = new ParameterStruct;
	parametersOnHost->numberofthreads         = 64;
	parametersOnHost->omega                   = (real)1.0/((real)3.0*this->viscosityLB+(real)0.5);
	parametersOnHost->isEvenTimestep          = true;

	//device
	parametersOnDevice                        = new ParameterStruct;
	parametersOnDevice->numberofthreads       = parametersOnHost->numberofthreads;
	parametersOnDevice->omega                 = parametersOnHost->omega;
	parametersOnDevice->isEvenTimestep        = parametersOnHost->isEvenTimestep;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//set-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameter::setlimitOfNodesForVTK(uint limitOfNodesForVTK)
{
	this->limitOfNodesForVTK = limitOfNodesForVTK;
}
void Parameter::setD3Qxx(int d3qxx)
{
	this->D3Qxx = d3qxx;
}
void Parameter::setTimestepEnd(uint timestepEnd)
{
	this->timestepEnd = timestepEnd;
}
void Parameter::setTimestepOut(uint timestepOut)
{
	this->timestepOut = timestepOut;
}
void Parameter::setTimestepStartOut(uint timestepStartOut)
{
	this->timestepStartOut = timestepStartOut;
}
void Parameter::setOutputPath(std::string outputPath)
{
	this->outputPath = outputPath;
}
void Parameter::setOutputPrefix(std::string outputPrefix)
{
	this->outputPrefix = outputPrefix;
}
void Parameter::setPathAndFilename(std::string pathAndFilename)
{
	this->pathAndFilename = pathAndFilename;
}
void Parameter::setPrintFiles(bool printfiles)
{
	this->printFiles = printfiles;
}
void Parameter::setViscosityLB(real viscosity)
{
	this->viscosityLB = viscosity;
}
void Parameter::setVelocityLB(real velocity)
{
	this->velocityLB = velocity;
}
void Parameter::setViscosityRatio(real viscosityRatio)
{
	this->viscosityRatio = viscosityRatio;
}
void Parameter::setVelocityRatio(real velocityRatio)
{
	this->velocityRatio = velocityRatio;
}
void Parameter::setDensityRatio(real densityRatio)
{
	this->densityRatio = densityRatio;
}
void Parameter::setPressureRatio(real pressureRatio)
{
	this->pressRatio = pressureRatio;
}
void Parameter::setRe(real Re)
{
	this->Re = Re;
}
void Parameter::setMemsizeGPU(double addMemory, bool reset)
{
	if (reset == true)
	{
		this->memsizeGPU = 0.;
	} 
	else
	{
		this->memsizeGPU += addMemory;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get-methods
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
uint Parameter::getlimitOfNodesForVTK()
{
	return this->limitOfNodesForVTK;
}
ParameterStruct* Parameter::getParD()
{
	return this->parametersOnDevice;
}
ParameterStruct* Parameter::getParH()
{
	return this->parametersOnHost;
}
bool Parameter::getIsTimestepEven()
{
	return this->parametersOnHost->isEvenTimestep;
}
int Parameter::getD3Qxx()
{
	return this->D3Qxx;
}
uint Parameter::getTimestepStart()
{
	return 1;
}
uint Parameter::getTimestepEnd()
{
	return this->timestepEnd;
}
uint Parameter::getTimestepOut()
{
	return this->timestepOut;
}
uint Parameter::getTimestepStartOut()
{
	return this->timestepStartOut;
}
std::string Parameter::getOutputPath()
{
	return this->outputPath;
}
std::string Parameter::getOutputPrefix()
{
	return this->outputPrefix;
}
std::string Parameter::getPathAndFilename()
{
	return this->pathAndFilename;
}
bool Parameter::getPrintFiles()
{
	return this->printFiles;
}
real Parameter::getViscosityLB()
{
	return this->viscosityLB;
}
real Parameter::getVelocityLB()
{
	return this->velocityLB;
}
real Parameter::getViscosityRatio()
{
	return this->viscosityRatio;
}
real Parameter::getVelocityRatio()
{
	return this->velocityRatio;
}
real Parameter::getDensityRatio()
{
	return this->densityRatio;
}
real Parameter::getPressureRatio()
{
	return this->pressRatio;
}
real Parameter::getRe()
{
	return this->Re;
}
double Parameter::getMemsizeGPU()
{
	return this->memsizeGPU;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
#ifndef PARAMETER_H
#define PARAMETER_H

#include "LBM/LB.h"
#include "Core/PointerDefinitions.h"
#include "VirtualFluidsDefinitions.h"

//! \struct ParameterStruct
//! \brief struct holds and manages the LB-parameter of the simulation
//! \brief For this purpose it holds structures and pointer for host and device data, respectively.
struct ParameterStruct{
	//////////////////////////////////////////////////////////////////////////
	//! \brief decides if the simulation timestep is even or odd
	//! \brief this information is important for the esoteric twist
	bool isEvenTimestep;
	//////////////////////////////////////////////////////////////////////////
	//! \brief stores the number of threads per GPU block
	uint numberofthreads;
	//////////////////////////////////////////////////////////////////////////
	//! \brief store all distribution functions for the D3Q27
	Distributions27 distributions;
	//////////////////////////////////////////////////////////////////////////
	//! \brief stores the type for every lattice node (f.e. fluid node)
	uint *typeOfGridNode;
	//////////////////////////////////////////////////////////////////////////
	//! \brief store the neighbors in +X, +Y, +Z, and in diagonal negative direction
	//! \brief this information is important because we use an indirect addressing scheme
	uint *neighborX, *neighborY, *neighborZ, *neighborInverse;
	//////////////////////////////////////////////////////////////////////////
	//! \brief store the coordinates for every lattice node 
	real *coordinateX, *coordinateY, *coordinateZ;
	//////////////////////////////////////////////////////////////////////////
	//! \brief store the macroscopic values (velocity, density, pressure)
	//! \brief for every lattice node
	real *velocityX, *velocityY, *velocityZ, *rho, *pressure;
	//! \brief stores the value for omega
	real omega;
	//////////////////////////////////////////////////////////////////////////
	//! \brief stores the number of nodes (based on indirect addressing scheme)
	uint numberOfNodes;
	//! \brief stores the size of the memory consumption for real/int values of the above arrays
	uint mem_size_real, mem_size_int;
	//////////////////////////////////////////////////////////////////////////
	//! \brief stores the velocity boundary condition data
	QforBoundaryConditions inflowBC;
	//! \brief number of lattice nodes for the velocity boundary condition
	uint numberOfInflowBCnodes;
	//////////////////////////////////////////////////////////////////////////
	//! \brief sets the forcing uniform on every fluid node in all three space dimensions 
	real *forcing;
};

//! \class Parameter implements the singleton design pattern.
//! \brief Class for LBM-parameter management
class VIRTUALFLUIDS_GPU_EXPORT Parameter
{
public:
	////////////////////////////////////////////////////////////////////////////
	//! \brief generate a new instance of parameter object
	static SPtr<Parameter> make();
	//! \brief returns the instance generate a new instance of parameter object
	static Parameter* getInstanz();
	//! \brief Pointer to instance of ParameterStruct - stored on Host System
	ParameterStruct* getParH();
	//! \brief Pointer to instance of ParameterStruct - stored on Device (GPU)
	ParameterStruct* getParD();

	//////////////////////////////////////////////////////////////////////////
	//! \brief initialization of necessary parameters at startup
	void initParameter();

	//////////////////////////////////////////////////////////////////////////
	//set methods
	//////////////////////////////////////////////////////////////////////////
	//! \brief sets the limit of nodes, that can be written to a binary unstructured grid VTK file  
	//! \param limitOfNodesForVTK holds the maximum number of nodes
	void setlimitOfNodesForVTK(uint limitOfNodesForVTK);
	//! \brief sets the LBM stencil
	//! \param d3qxx holds the number of neighbors (f.e. 27)
	void setD3Qxx(int d3qxx);
	//! \brief sets timestep to stop the simulation
	//! \param timestepEnd holds the last timestep of the simulation
	void setTimestepEnd(uint timestepEnd);
	//! \brief sets time interval to write output files
	//! \param timestepOut holds the value for the output interval
	void setTimestepOut(uint timestepOut);
	//! \brief sets first timestep to write output files
	//! \param timestepStartOut holds the value for the first output timestep
	void setTimestepStartOut(uint timestepStartOut);
	//! \brief sets the path, where the vtk-files are stored 
	//! \param string "oPath" represents the output path
	void setOutputPath(std::string outputPath);
	//! \brief sets the prefix of the vtk-files name 
	//! \param string "oPrefix" represents the file-name-prefix
	void setOutputPrefix(std::string outputPrefix);
	//! \brief sets the complete string for the vtk-files with results 
	//! \param string "fname" represents the combination of path and prefix
	void setPathAndFilename(std::string pathAndFilename);
	//! \brief sets the status, if the vtk files should be printed
	//! \param if printfiles is true, the vtk files will be printed 
	void setPrintFiles(bool printfiles);
	//! \brief sets the viscosity in LB units
	//! \param viscosity in LB units 
	void setViscosityLB(real viscosity);
	//! \brief sets the velocity in LB units
	//! \param velocity in LB units 
	void setVelocityLB(real velocity);
	//! \brief sets the viscosity ratio between SI and LB units
	//! \param viscosityRatio SI/LB units 
	void setViscosityRatio(real viscosityRatio);
	//! \brief sets the velocity ratio between SI and LB units
	//! \param velocityRatio SI/LB units 
	void setVelocityRatio(real velocityRatio);
	//! \brief sets the density ratio between SI and LB units
	//! \param densityRatio SI/LB units 
	void setDensityRatio(real densityRatio);
	//! \brief sets the pressure ratio between SI and LB units
	//! \param pressureRatio SI/LB units 
	void setPressureRatio(real pressureRatio);
	//! \brief sets the Reynolds number (Re) for the simulation
	//! \param Reynolds number (Re) 
	void setRe(real Re);
	//! \brief sets the necessary memory on the device(s)/GPU(s)
	//! \param addMemory is the amount of additional memory 
	//! \param reset decides if the value for overall GPU memory should be set to zero 
	void setMemsizeGPU(double addMemory, bool reset);

	//////////////////////////////////////////////////////////////////////////
	//get methods
	//////////////////////////////////////////////////////////////////////////
	//! \brief return the limit of nodes, that can be written to a binary unstructured grid VTK file  
	uint getlimitOfNodesForVTK();
	//! \brief return if the timestep is even or odd
	bool getIsTimestepEven();
	//! \brief return if the simulation should write VTK files
	bool getPrintFiles();
	//! \brief return the number of neighbors of a lattice node (stencil)
	int getD3Qxx();
	//! \brief return the output path
	std::string getOutputPath();
	//! \brief return the prefix of the output files
	std::string getOutputPrefix();
	//! \brief return the combination of output path and prefix
	std::string getPathAndFilename();
	//! \brief return the timestep to start the simulation (in this code version = 1)
	uint getTimestepStart();
	//! \brief return the timestep to end the simulation
	uint getTimestepEnd();
	//! \brief return the time interval to write output files
	uint getTimestepOut();
	//! \brief return the timestep to start writing output files
	uint getTimestepStartOut();
	//! \brief return the viscosity in LB units
	real getViscosityLB();
	//! \brief return the velocity in LB units
	real getVelocityLB();
	//! \brief return the viscosity ratio in SI/LB units
	real getViscosityRatio();
	//! \brief return the velocity ratio in SI/LB units
	real getVelocityRatio();
	//! \brief return the density ratio in SI/LB units
	real getDensityRatio();
	//! \brief return the pressure ratio in SI/LB units
	real getPressureRatio();
	//! \brief return the Reynolds number
	real getRe();
	//! \brief return the used device memory
	double getMemsizeGPU();

	//////////////////////////////////////////////////////////////////////////
	//! Class destructor
	~Parameter();
protected:
private:
	//! \brief instance of parameter object 
	static Parameter* instance;
	//! \brief stencil for the LB simulation, number of node neighbors
	int D3Qxx;
	//! \brief limit of nodes, that can be written to a binary unstructured grid VTK file  
	uint limitOfNodesForVTK;
	//! \brief last timestep of the simulation
	uint timestepEnd;
	//! \brief time interval to write output files
	uint timestepOut;
	//! \brief timestep - start writing output files
	uint timestepStartOut;
	//! \brief Reynolds number
	real Re;
	//! \brief viscosity and velocity in LB units
	real viscosityLB, velocityLB;
	//! \brief ratio SI units / LB units for viscosity, velocity, density and pressure
	real viscosityRatio, velocityRatio;
	real densityRatio, pressRatio;
	//! \brief used device memory
	double memsizeGPU;
	//! \brief write output files on/off
	bool printFiles;
	//! \brief strings to store output path, prefix and combination of both
	std::string pathAndFilename, outputPath, outputPrefix;

	//! \brief pointer to LB-parameter struct on host system
	ParameterStruct* parametersOnHost;
	//! \brief pointer to LB-parameter struct on device/GPU
	ParameterStruct* parametersOnDevice;

	//! Class default constructor
	Parameter();
	//Parameter(const Parameter&);
};

#endif


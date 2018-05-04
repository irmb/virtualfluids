#include "FileWriter.h"

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>

#include <utilities/StringUtil/StringUtil.h>

#include "Parameter/Parameter.h"

#include "LBM/LB.h"
#include "LBM/D3Q27.h"

#include <VirtualFluidsBasics/basics/writer/WbWriterVtkXmlBinary.h>

#include "UnstructuredGridWriter.hpp"

void FileWriter::writeInit(std::shared_ptr<Parameter> para)
{
    unsigned int timestep = para->getTInit();
    for (int level = para->getCoarse(); level <= para->getFine(); level++)
        writeTimestep(para, timestep, level);
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    for (int level = para->getCoarse(); level <= para->getFine(); level++)
        writeTimestep(para, timestep, level);
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level)
{
    const unsigned int numberOfParts = para->getParH(level)->size_Mat_SP / para->getlimitOfNodesForVTK() + 1;
    std::vector<std::string> fname;
    for (unsigned int i = 1; i <= numberOfParts; i++)
        fname.push_back(para->getFName() + "_bin_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(para->getMyID()) + "_Part_" + StringUtil::toString<int>(i) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk");
    
    UnstructuredGridWriter::writeUnstrucuredGridLT(para.get(), level, fname);
}

void FileWriter::writeParticle(Parameter* para, unsigned int t)
{
	for (int lev = para->getParticleBasicLevel(); lev <= para->getFine(); lev++)
	{
		//////////////////////////////////////////////////////////////////////////
		//set filename
		std::string fname = para->getFName()+StringUtil::toString<int>(lev)+StringUtil::toString<int>(para->getMyID())+"_t_"+StringUtil::toString<int>(t)+"_Particles.vtk";
		//////////////////////////////////////////////////////////////////////////
		//write particles
		UnstructuredGridWriter::writeUnstrucuredParticles(para, lev, fname);
		//////////////////////////////////////////////////////////////////////////
	}
}
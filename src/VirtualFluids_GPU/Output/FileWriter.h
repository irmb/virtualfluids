#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <memory>

#include <VirtualFluidsDefinitions.h>

#include "DataWriter.h"

class Parameter;
class CudaMemoryManager;

class FileWriter : public DataWriter
{
public:
	VF_PUBLIC FileWriter() {}

	void VF_PUBLIC writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager) override;
	void VF_PUBLIC writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep) override;

private:
	void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level);
	void writeUnstrucuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string>& fname);
	
	FileWriter(const FileWriter& fileWriter) {};
};
#endif
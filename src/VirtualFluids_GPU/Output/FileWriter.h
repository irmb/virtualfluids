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
	void VF_PUBLIC writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level) override;
	//void VF_PUBLIC writeParticle(Parameter* para, unsigned int t);
	void VF_PUBLIC writeUnstrucuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	void VF_PUBLIC writeUnstrucuredGridLTConc(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	void VF_PUBLIC writeUnstrucuredGridMedianLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	void VF_PUBLIC writeUnstrucuredGridMedianLTConc(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	bool VF_PUBLIC isPeriodicCell(std::shared_ptr<Parameter> para, int level, unsigned int number2, unsigned int number1, unsigned int number3, unsigned int number5);

	FileWriter(const FileWriter& fileWriter) {};
};
#endif
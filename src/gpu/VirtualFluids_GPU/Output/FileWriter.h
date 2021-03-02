#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <memory>
#include <vector>
#include <string>



#include "DataWriter.h"

class Parameter;
class CudaMemoryManager;

class FileWriter : public DataWriter
{
public:
	VIRTUALFLUIDS_GPU_EXPORT FileWriter() {}

	void VIRTUALFLUIDS_GPU_EXPORT writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager) override;
	void VIRTUALFLUIDS_GPU_EXPORT writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep) override;

private:
	void VIRTUALFLUIDS_GPU_EXPORT writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level) override;
	//void VIRTUALFLUIDS_GPU_EXPORT writeParticle(Parameter* para, unsigned int t);
	void VIRTUALFLUIDS_GPU_EXPORT writeUnstrucuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	void VIRTUALFLUIDS_GPU_EXPORT writeUnstrucuredGridLTConc(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	void VIRTUALFLUIDS_GPU_EXPORT writeUnstrucuredGridMedianLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	void VIRTUALFLUIDS_GPU_EXPORT writeUnstrucuredGridMedianLTConc(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);
	bool VIRTUALFLUIDS_GPU_EXPORT isPeriodicCell(std::shared_ptr<Parameter> para, int level, unsigned int number2, unsigned int number1, unsigned int number3, unsigned int number5);

	FileWriter(const FileWriter& fileWriter) = default;

    void VIRTUALFLUIDS_GPU_EXPORT writeCollectionFile( std::shared_ptr<Parameter> para, unsigned int timestep );

    void VIRTUALFLUIDS_GPU_EXPORT writeCollectionFileMedian( std::shared_ptr<Parameter> para, unsigned int timestep );

    std::vector< std::string > fileNamesForCollectionFile;
    std::vector< std::string > fileNamesForCollectionFileMedian;
};
#endif
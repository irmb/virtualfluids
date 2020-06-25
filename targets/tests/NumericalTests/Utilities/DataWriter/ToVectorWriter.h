#ifndef WRITE_TO_VECTOR_H
#define WRITE_TO_VECTOR_H

#include "VirtualFluids_GPU/Output/DataWriter.h"

class Parameter;
class FileWriter;
struct VectorWriterInformationStruct;

class ToVectorWriter : public DataWriter
{
public:
	void writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager);
	void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t);
	
protected:
	ToVectorWriter();
	ToVectorWriter(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength);
	virtual void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t, int level) = 0;

	std::shared_ptr<FileWriter> vtkFileWriter;
	bool writeVTKFiles;
	unsigned int timeStepLength;
	unsigned int startTimeVectorWriter, startTimeVTKWriter;
};
#endif 
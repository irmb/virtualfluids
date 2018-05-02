#ifndef WRITE_TO_VECTOR_H
#define WRITE_TO_VECTOR_H

#include "VirtualFluids_GPU/Output/DataWriter.h"

class Parameter;
class Results;
class FileWriter;

class ToVectorWriter : public DataWriter
{
public:
	ToVectorWriter(unsigned int ySliceForCalculation, unsigned int startTimeY2dSliceToVector, unsigned int endTime, unsigned int timeStepLength, bool writeFiles, std::shared_ptr<FileWriter> fileWriter, unsigned int startTimeDataWriter);

	void writeInit(std::shared_ptr<Parameter> para);
	void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t);

protected:
	virtual void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t, int level) = 0;
	void writeUnstrucuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname) {};

	std::shared_ptr<FileWriter> fileWriter;
	bool writeFiles;

	unsigned int ySliceForCalculation;
	unsigned int counterTimeSteps;
	unsigned int startTimeY2dSliceToVector, startTimeDataWriter;
	unsigned int endTime;
	unsigned int maxX, maxY, maxZ;
private:

};
#endif 
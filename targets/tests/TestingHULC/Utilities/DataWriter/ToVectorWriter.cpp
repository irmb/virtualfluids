#include "ToVectorWriter.h"

#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"

ToVectorWriter::ToVectorWriter(unsigned int ySliceForCalculation, unsigned int startTimeY2dSliceToVector, unsigned int endTime, unsigned int timeStepLength, bool writeFiles, std::shared_ptr<FileWriter> fileWriter, unsigned int startTimeDataWriter)
{
	this->writeFiles = writeFiles;
	this->fileWriter = fileWriter;
	this->ySliceForCalculation = ySliceForCalculation;
	this->startTimeY2dSliceToVector = startTimeY2dSliceToVector;
	this->startTimeDataWriter = startTimeDataWriter;
	this->endTime = endTime;
}

void ToVectorWriter::writeInit(std::shared_ptr<Parameter> para)
{
	if (startTimeY2dSliceToVector == 0)
		writeTimestep(para, 0);
	if (writeFiles && startTimeDataWriter == 0)
		fileWriter->writeTimestep(para, 0);
}

void ToVectorWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int t)
{
	if (startTimeY2dSliceToVector <= t && endTime >= t)
	{
		for (int level = para->getCoarse(); level <= para->getFine(); level++)
		{
			writeTimestep(para, t, level);
		}
	}
	if (writeFiles && startTimeDataWriter < t)
		fileWriter->writeTimestep(para, t);
}

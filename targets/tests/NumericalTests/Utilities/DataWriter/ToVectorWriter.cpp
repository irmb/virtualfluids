#include "ToVectorWriter.h"

#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"

#include "Utilities/Structs/VectorWriterInformationStruct.h"

ToVectorWriter::ToVectorWriter(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength)
{
	this->startTimeVectorWriter = vectorWriterInfo->startTimeVectorWriter;
	this->writeVTKFiles = vectorWriterInfo->writeVTKFiles;
	this->startTimeVTKWriter = vectorWriterInfo->startTimeVTKDataWriter;

	this->vtkFileWriter = std::shared_ptr<FileWriter> (new FileWriter());

	this->timeStepLength = timeStepLength;
}

void ToVectorWriter::writeInit(std::shared_ptr<Parameter> para)
{
	if (startTimeVectorWriter == 0)
		writeTimestep(para, 0);
	if (writeVTKFiles && startTimeVTKWriter == 0)
		vtkFileWriter->writeTimestep(para, 0);
}

void ToVectorWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int t)
{
	if (startTimeVectorWriter <= t)
	{
		for (int level = para->getCoarse(); level <= para->getFine(); level++)
		{
			writeTimestep(para, t, level);
		}
	}
	if (writeVTKFiles && startTimeVTKWriter < t)
		vtkFileWriter->writeTimestep(para, t);
}

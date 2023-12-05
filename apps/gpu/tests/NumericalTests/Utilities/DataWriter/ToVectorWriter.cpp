#include "ToVectorWriter.h"

#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/Output/FileWriter.h"

#include "Utilities/Structs/VectorWriterInformationStruct.h"

#include "Cuda/CudaMemoryManager.h"

ToVectorWriter::ToVectorWriter(std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int timeStepLength)
{
    this->startTimeVectorWriter = vectorWriterInfo->startTimeVectorWriter;
    this->writeVTKFiles = vectorWriterInfo->writeVTKFiles;
    this->startTimeVTKWriter = vectorWriterInfo->startTimeVTKDataWriter;

    this->vtkFileWriter = std::shared_ptr<FileWriter> (new FileWriter());

    this->timeStepLength = timeStepLength;
}

void ToVectorWriter::writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager)
{
    if (startTimeVectorWriter == 0) {
        for (int level = para->getCoarse(); level <= para->getFine(); level++)
            cudaManager->cudaCopyPrint(level);
        writeTimestep(para, 0);
    }
        
    if (writeVTKFiles && startTimeVTKWriter == 0)
        vtkFileWriter->writeTimestep(para, 0);
}

void ToVectorWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int t)
{
    if (startTimeVectorWriter <= t)
    {
        for (int level = para->getCoarse(); level <= para->getFine(); level++)                
            writeTimestep(para, t, level);
    }
    if (writeVTKFiles && startTimeVTKWriter < t)
        vtkFileWriter->writeTimestep(para, t);
}

ToVectorWriter::ToVectorWriter()
{
}

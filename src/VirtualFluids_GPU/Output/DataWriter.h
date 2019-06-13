#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#define NOMINMAX

#include <memory>
#include <vector>


class Parameter;
class CudaMemoryManager;

class DataWriter
{
public:
	DataWriter() {}
    virtual ~DataWriter() {}

    virtual void writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager) = 0;
    virtual void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t) = 0;

    DataWriter(const DataWriter& dataWriter) {}
};
#endif
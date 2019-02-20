#ifndef DATA_WRITER_H
#define DATA_WRITER_H


#include <memory>
#include <vector>


class Parameter;

class DataWriter
{
public:
	DataWriter() {}
    virtual ~DataWriter() {}

    virtual void writeInit(std::shared_ptr<Parameter> para) = 0;
    virtual void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep) = 0;
	virtual void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level) = 0;

    DataWriter(const DataWriter& dataWriter) {}
};
#endif

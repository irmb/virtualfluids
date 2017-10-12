#ifndef WRITE_DATA_H
#define WRITE_DATA_H

#define NOMINMAX

#include <memory>
#include <vector>

class Parameter;

class DataWriter
{
public:
    virtual ~DataWriter() {};

    static void writeInit(std::shared_ptr<Parameter> para);
    static void writeTimestep(std::shared_ptr<Parameter> para, unsigned int t);


private:
    static void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level);
    static void writeUnstrucuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);

    DataWriter() {};
    DataWriter(const DataWriter& dataWriter) {};
};

#endif

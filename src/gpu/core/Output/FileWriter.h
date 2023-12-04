#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <memory>
#include <vector>
#include <string>
#include <map>


#include "DataWriter.h"

class Parameter;
class CudaMemoryManager;

class FileWriter : public DataWriter
{
public:
    void writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager) override;
    void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep) override;

private:
    void writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level) override;
    std::vector<std::string> writeUnstructuredGridLT(std::shared_ptr<Parameter> para, int level,
                                                         std::vector<std::string> &fname);
    std::vector<std::string> writeUnstructuredGridMeanLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname);

    std::string writeCollectionFile( std::shared_ptr<Parameter> para, unsigned int timestep );

    std::string writeCollectionFileMean( std::shared_ptr<Parameter> para, unsigned int timestep );

    std::string writePvdCollectionFileForTimeSeries(const Parameter &para);

    std::vector<std::string> getNodeDataNames(std::shared_ptr<Parameter> para);
    std::vector<std::string> getMeanNodeDataNames(std::shared_ptr<Parameter> para);

    std::vector< std::string > fileNamesForCollectionFile;
    std::vector< std::string > fileNamesForCollectionFileMean;

    std::map<uint, std::vector<std::string>> fileNamesForCollectionFileTimeSeries; // key: timeStep, value: fileNames for timeStep
};
#endif

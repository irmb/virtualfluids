#ifndef NUMERICAL_TESTS_GRID_READER_H
#define NUMERICAL_TESTS_GRID_READER_H

#include "gpu/core/DataStructureInitializer/GridReaderFiles/GridReader.h"

#include <memory>

class InitialCondition;

class NumericalTestGridReader : public GridReader
{
public:
    static std::shared_ptr<NumericalTestGridReader> getNewInstance(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager);

protected:
    virtual void setInitialNodeValues(uint numberOfNodes, int level) const override;
    
private:
    NumericalTestGridReader(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition, std::shared_ptr<CudaMemoryManager> cudaManager);

    std::shared_ptr<InitialCondition> initialCondition;
    
};
#endif

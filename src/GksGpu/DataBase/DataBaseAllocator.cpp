#include "DataBaseAllocator.h"

//#include "../../DataBase/DataBaseAllocator/DataBaseAllocatorCPU/DataBaseAllocatorCPU.h"
//#include "../../DataBase/DataBaseAllocator/DataBaseAllocatorGPU/DataBaseAllocatorGPU.h"

#include "DataBaseAllocatorCPU.h"
#include "DataBaseAllocatorGPU.h"

#include <string>

namespace GksGpu {

std::shared_ptr<DataBaseAllocator> DataBaseAllocator::create(std::string type)
{
    if ( type == "GPU" )
        return std::shared_ptr<DataBaseAllocator>( new DataBaseAllocatorGPU() );
    else
        return std::shared_ptr<DataBaseAllocator>( new DataBaseAllocatorCPU() );
}

DataBaseAllocator::~DataBaseAllocator()
{
}

DataBaseAllocator::DataBaseAllocator()
{
}

DataBaseAllocator::DataBaseAllocator(const DataBaseAllocator & orig)
{
}

} // namespace GksGpu

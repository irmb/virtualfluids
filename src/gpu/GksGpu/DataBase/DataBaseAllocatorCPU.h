#ifndef DataBaseAllocatorCPU_H
#define DatabaseAllocatorCPU_H

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "DataBaseAllocator.h"

#include "VirtualFluidsDefinitions.h"

namespace GksGpu {

class VIRTUALFLUIDS_GPU_EXPORT DataBaseAllocatorCPU : public DataBaseAllocator {

public:

    virtual void freeMemory( DataBase& dataBase ) override;

    virtual void allocateMemory( SPtr<DataBase> dataBase ) override;

    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) override;

    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) override;
    
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* dataHost ) override;

    virtual int  getCrashCellIndex( SPtr<DataBase> dataBase ) override;

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( BoundaryCondition& boundaryCondition ) override;

    virtual void allocateMemory( SPtr<BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells ) override;

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( Communicator& communicator ) override;

    virtual void allocateMemory( Communicator& communicator, std::vector<uint>& sendIndices, std::vector<uint>& recvIndices ) override;

    virtual void copyDataDeviceToDevice( SPtr<Communicator> dst, SPtr<Communicator> src ) override;

    virtual void copyBuffersDeviceToHost( SPtr<Communicator> communicator ) override;
    virtual void copyBuffersHostToDevice( SPtr<Communicator> communicator ) override;

    //////////////////////////////////////////////////////////////////////////

    virtual std::string getDeviceType() override;
};

} // namespace GksGpu


#endif
#ifndef DataBaseAllocatorGPU_H
#define DatabaseAllocatorGPU_H

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "DataBaseAllocator.h"

#include "VirtualFluidsDefinitions.h"

class VF_PUBLIC DataBaseAllocatorGPU : public DataBaseAllocator {

public:

    virtual void freeMemory( DataBase& dataBase ) override;

    virtual void allocateMemory( SPtr<DataBase> dataBase ) override;

    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) override;

    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) override;
    
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* dataHost ) override;

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( BoundaryCondition& boundaryCondition ) override;

    virtual void allocateMemory( SPtr<BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells ) override;

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( Communicator& communicator ) override;

    virtual void allocateMemory( Communicator& communicator, std::vector<uint>& sendIndices, std::vector<uint>& recvIndices ) override;

    virtual void copyDataDeviceToDevice( SPtr<Communicator> dst, SPtr<Communicator> src ) override;

    //////////////////////////////////////////////////////////////////////////

    virtual std::string getDeviceType() override;
};


#endif
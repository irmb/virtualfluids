#ifndef DataBaseAllocatorGPU_H
#define DatabaseAllocatorGPU_H

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "DataBaseAllocator.h"

#include <VirtualFluidsDefinitions.h>

class VF_PUBLIC DataBaseAllocatorGPU : public DataBaseAllocator {

public:

    virtual void freeMemory( SPtr<DataBase> dataBase ) override;

    virtual void allocateMemory( SPtr<DataBase> dataBase ) override;

    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) override;

    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) override;
    
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* dataHost ) override;
};


#endif
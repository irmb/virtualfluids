#ifndef DataBaseAllocator_H
#define DataBaseAllocator_H

#include <string>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include <VirtualFluidsDefinitions.h>

class  GksMeshAdapter;
struct DataBase;

class VF_PUBLIC DataBaseAllocator {

public:
    virtual void freeMemory( SPtr<DataBase> dataBase ) = 0;

    virtual void allocateMemory( SPtr<DataBase> dataBase) = 0;

    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) = 0;

    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) = 0;
    
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* hostData ) = 0;

    static std::shared_ptr<DataBaseAllocator> create( std::string type );

    ~DataBaseAllocator();

protected:

    DataBaseAllocator();
    DataBaseAllocator( const DataBaseAllocator& orig );

};


#endif
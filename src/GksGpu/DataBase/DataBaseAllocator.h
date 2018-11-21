#ifndef DataBaseAllocator_H
#define DataBaseAllocator_H

#include <string>
#include <vector>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "VirtualFluidsDefinitions.h"

class  GksMeshAdapter;
struct DataBase;
struct BoundaryCondition;

class VF_PUBLIC DataBaseAllocator {

public:

    static std::shared_ptr<DataBaseAllocator> create( std::string type );

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( DataBase& dataBase ) = 0;

    virtual void allocateMemory( SPtr<DataBase> dataBase) = 0;

    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) = 0;

    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) = 0;
    
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* hostData ) = 0;

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( BoundaryCondition& boundaryCondition ) = 0;

    virtual void allocateMemory( SPtr<BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells ) = 0;

    //////////////////////////////////////////////////////////////////////////

    ~DataBaseAllocator();

    virtual std::string getDeviceType() = 0;

protected:

    DataBaseAllocator();
    DataBaseAllocator( const DataBaseAllocator& orig );

};


#endif
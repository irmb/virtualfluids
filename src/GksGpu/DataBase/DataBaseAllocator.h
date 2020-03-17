#ifndef DataBaseAllocator_H
#define DataBaseAllocator_H

#include <string>
#include <vector>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "VirtualFluidsDefinitions.h"

class  GksMeshAdapter;
struct DataBase;
namespace GksGpu { struct BoundaryCondition; };
struct Communicator;

class VF_PUBLIC DataBaseAllocator {

public:

    static std::shared_ptr<DataBaseAllocator> create( std::string type );

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( DataBase& dataBase ) = 0;

    virtual void allocateMemory( SPtr<DataBase> dataBase) = 0;

    virtual void copyMesh( SPtr<DataBase> dataBase, GksMeshAdapter& adapter ) = 0;

    virtual void copyDataHostToDevice( SPtr<DataBase> dataBase ) = 0;
    
    virtual void copyDataDeviceToHost( SPtr<DataBase> dataBase, real* hostData ) = 0;

    virtual int  getCrashCellIndex( SPtr<DataBase> dataBase ) = 0;

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( GksGpu::BoundaryCondition& boundaryCondition ) = 0;

    virtual void allocateMemory( SPtr<GksGpu::BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells ) = 0;

    //////////////////////////////////////////////////////////////////////////

    virtual void freeMemory( Communicator& communicator ) = 0;

    virtual void allocateMemory( Communicator& communicator, std::vector<uint>& sendIndices, std::vector<uint>& recvIndices ) = 0;

    virtual void copyDataDeviceToDevice( SPtr<Communicator> dst, SPtr<Communicator> src ) = 0;

    virtual void copyBuffersDeviceToHost( SPtr<Communicator> communicator ) = 0;
    virtual void copyBuffersHostToDevice( SPtr<Communicator> communicator ) = 0;

    //////////////////////////////////////////////////////////////////////////

    ~DataBaseAllocator();

    virtual std::string getDeviceType() = 0;

protected:

    DataBaseAllocator();
    DataBaseAllocator( const DataBaseAllocator& orig );

};


#endif
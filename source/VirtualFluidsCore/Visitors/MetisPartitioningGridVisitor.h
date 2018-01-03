#ifndef MetisPartitioningGridVisitor_h 
#define MetisPartitioningGridVisitor_h

#if defined VF_METIS && defined VF_MPI

#include <vector>
#include <memory>

#include "Grid3DVisitor.h"
#include "MetisPartitioner.h"


class Communicator;

class MetisPartitioningGridVisitor;
typedef std::shared_ptr<MetisPartitioningGridVisitor> PartitionMetisGridVisitorPtr;

////////////////////////////////////////////////////////////////////////
//! \brief The class implements domain decomposition with METIS library
//! \author Kostyantyn Kucher
//////////////////////////////////////////////////////////////////////////

class Grid3D;
class Block3D;

class MetisPartitioningGridVisitor : public Grid3DVisitor
{                                             
public:
   //! This describe different types of decomposition   
   enum GraphType{LevelIntersected, LevelBased};

public:
   //! Constructor
   //! \param comm - communicator
   //! \param graphType - type of decomposition
   //! \param numOfDirs - maximum number of neighbors for each process
   //! \param threads - on/off decomposition for threads
   //! \param numberOfThreads - number of threads
   MetisPartitioningGridVisitor(std::shared_ptr<Communicator> comm, GraphType graphType, int numOfDirs, MetisPartitioner::PartType partType = MetisPartitioner::KWAY, bool threads = false, int numberOfThreads = 0);
   virtual ~MetisPartitioningGridVisitor();
   void visit(std::shared_ptr<Grid3D> grid) override;
   void setNumberOfProcesses(int np);

protected:
   enum PartLevel {BUNDLE, PROCESS, THREAD};
   void collectData(std::shared_ptr<Grid3D> grid, int nofSegments, PartLevel level);
   void buildMetisGraphLevelIntersected(std::shared_ptr<Grid3D> grid, int nofSegments, PartLevel level);
   void buildMetisGraphLevelBased(std::shared_ptr<Grid3D> grid, int nofSegments, PartLevel level);
   bool getPartitionCondition(std::shared_ptr<Block3D> block, PartLevel level);
   void distributePartitionData(std::shared_ptr<Grid3D> grid, PartLevel level);
   void clear();
   int  nofSegments;
   int numOfDirs;
   std::vector<int> blockID;
   std::vector<idx_t> parts;
   std::shared_ptr<Communicator> comm;
   int bundleRoot;
   int processRoot;
   int bundleID;
   int processID;
   int numberOfBundles;
   int numberOfThreads;
   bool threads;
   GraphType graphType;
   MetisPartitioner::PartType partType;
   int numberOfProcesses;
};

#endif  //VF_MPI
#endif 

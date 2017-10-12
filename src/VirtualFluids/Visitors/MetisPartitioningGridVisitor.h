#ifndef MetisPartitioningGridVisitor_h 
#define MetisPartitioningGridVisitor_h


#if defined VF_METIS && defined VF_MPI

#include <vector>

#include "MetisPartitioner.h"
#include "Grid3DVisitor.h"
#include "LBMSystem.h"
#include "Block3D.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class MetisPartitioningGridVisitor;
typedef boost::shared_ptr<MetisPartitioningGridVisitor> PartitionMetisGridVisitorPtr;

////////////////////////////////////////////////////////////////////////
//! \brief The class implements domain decomposition with METIS library
//! \author Kostyantyn Kucher
//////////////////////////////////////////////////////////////////////////
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
   MetisPartitioningGridVisitor(CommunicatorPtr comm, GraphType graphType, int numOfDirs, MetisPartitioner::PartType partType = MetisPartitioner::KWAY, bool threads = false, int numberOfThreads = 0);
   virtual ~MetisPartitioningGridVisitor();
   void visit(Grid3DPtr grid);
   void setNumberOfProcesses(int np);
protected:
   enum PartLevel {BUNDLE, PROCESS, THREAD};
   void collectData(Grid3DPtr grid, int nofSegments, PartLevel level);
   void buildMetisGraphLevelIntersected(Grid3DPtr grid, int nofSegments, PartLevel level);
   void buildMetisGraphLevelBased(Grid3DPtr grid, int nofSegments, PartLevel level);
   bool getPartitionCondition(Block3DPtr block, PartLevel level);
   void distributePartitionData(Grid3DPtr grid, PartLevel level);
   void clear();
   int  nofSegments;
   int numOfDirs;
   std::vector<int> blockID;
   std::vector<idx_t> parts;
   CommunicatorPtr comm;
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

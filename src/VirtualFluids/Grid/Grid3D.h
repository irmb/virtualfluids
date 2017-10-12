#ifndef GRID3D_H
#define GRID3D_H

#include <sstream>
#include <iostream>
#include <vector>
#include <set>

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/signals2/signal.hpp>
#include <boost/bind.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include "UbKeys.h"
#include "CoordinateTransformation3D.h"

#include <Communicator.h>
//#include <Observer.h>
//#include <FastSignal.hpp>

class Grid3D;
typedef boost::shared_ptr<Grid3D> Grid3DPtr;
typedef boost::weak_ptr<Grid3D>   Grid3DWeakPtr;

#include <Block3D.h>
#include <Block3DVisitor.h>
#include <Grid3DVisitor.h>
#include <Interactor3D.h>

//////////////////////////////////////////////////////////////////////////
class Grid3D : public boost::enable_shared_from_this<Grid3D>
{
public:
   typedef UbKeys::Key3<int>                  Block3DKey;
   typedef std::map< Block3DKey, Block3DPtr > Block3DMap;
   typedef std::map< int, Block3DPtr >    BlockIDMap;
   typedef std::vector<Block3DMap>        LevelSet;
   typedef std::vector<Interactor3DPtr>   Interactor3DSet;
   //typedef std::set<ObserverPtr>     GridObserversSet;

   typedef boost::signals2::signal<void (double)>  signal_t;
   typedef boost::signals2::connection  connection_t;

public:
   Grid3D();
   Grid3D(CommunicatorPtr comm);
   Grid3D(CommunicatorPtr comm, int blockNx1, int blockNx2, int blockNx3, int gridNx1, int gridNx2, int gridNx3);
   virtual ~Grid3D(){}
   //////////////////////////////////////////////////////////////////////////
   //blocks control
   void addBlock(Block3DPtr block);
   bool deleteBlock(Block3DPtr block);
   bool deleteBlock(int ix1, int ix2, int ix3, int level);
   void deleteBlocks(const std::vector<int>& ids);
   void replaceBlock(Block3DPtr block);
   Block3DPtr getBlock(int ix1, int ix2, int ix3, int level) const;
   Block3DPtr getBlock(int id) const;
   void getBlocksByCuboid(double minX1, double minX2, double minX3, 
                          double maxX1, double maxX2, double maxX3, 
                          std::vector<Block3DPtr>& blocks);
   void getBlocksByCuboid(int level, double minX1, double minX2, double minX3, 
                          double maxX1, double maxX2, double maxX3, 
                          std::vector<Block3DPtr>& blocks);
   //get blocks for level
   void getBlocks(int level, std::vector<Block3DPtr>& blockVector);
   //get blocks for level with current rank
   void getBlocks(int level, int rank, std::vector<Block3DPtr>& blockVector);
   //get only active or not active blocks 
   void getBlocks(int level, int rank, bool active, std::vector<Block3DPtr>& blockVector);
   int getNumberOfBlocks();
   int getNumberOfBlocks(int level);
   //const Block3DMap& getBlocks(int level);
   const BlockIDMap& getBlockIDs();
   void deleteBlockIDs();
   void updateDistributedBlocks(CommunicatorPtr comm);
   Block3DPtr getSuperBlock(Block3DPtr block);
   Block3DPtr getSuperBlock(int ix1, int ix2, int ix3, int level);
   void getSubBlocks(Block3DPtr block, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getSubBlocks(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blockVector);
   Block3DPtr getNeighborBlock(int dir, int ix1, int ix2, int ix3, int level) const;
   Block3DPtr getNeighborBlock(int dir, Block3DPtr block) const;
   bool expandBlock(int ix1, int ix2, int ix3, int level);
   Block3DPtr collapseBlock(int fix1, int fix2, int fix3, int flevel, int levelDepth);
   void getAllNeighbors(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getAllNeighbors(Block3DPtr block, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborBlocksForDirection(int dir, int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborBlocksForDirectionWithDirZero(int dir, int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);

   void getNeighborsZero(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsNorth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsSouth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsTop(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottom(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);

   void getNeighborsNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);

   void getNeighborsTopNorth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsTopSouth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsTopEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsTopWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);

   void getNeighborsBottomNorth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottomSouth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottomEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottomWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);

   void getNeighborsTopNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsTopNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsTopSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsTopSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottomNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottomNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottomSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   void getNeighborsBottomSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<Block3DPtr>& blocks);
   //////////////////////////////////////////////////////////////////////////
   //level control
   int getFinestInitializedLevel();
   int getCoarsestInitializedLevel();
   //////////////////////////////////////////////////////////////////////////
   void deleteConnectors();
   //////////////////////////////////////////////////////////////////////////
   //interactors control
   void addInteractor(Interactor3DPtr interactor);
   void addAndInitInteractor(Interactor3DPtr interactor, double timestep=0);
   Interactor3DSet getInteractors();
   //////////////////////////////////////////////////////////////////////////
   //visitors
   void accept(Block3DVisitor& blockVisitor);
   void accept(Grid3DVisitor& gridVisitor);
   void accept(Grid3DVisitorPtr gridVisitor);
   //////////////////////////////////////////////////////////////////////////
   //post processing
   connection_t connect(signal_t::slot_function_type subscriber);
   void disconnect(connection_t subscriber);
   void coProcess(double step);
   //////////////////////////////////////////////////////////////////////////
   //bundle and rank for distributed memory
   void setBundle(int bundle);
   int  getBundle() const;
   int  getRank() const;
   void setRank(int rank);
   //////////////////////////////////////////////////////////////////////////
   //periodic boundary
   bool isPeriodicX1() const;
   bool isPeriodicX2() const;
   bool isPeriodicX3() const;
   void setPeriodicX1(bool value);
   void setPeriodicX2(bool value);
   void setPeriodicX3(bool value);
   //////////////////////////////////////////////////////////////////////////
   //Topology
   UbTupleInt3 getBlockIndexes(double blockX1Coord, double blockX2Coord, double blockX3Coord)  const;
   UbTupleInt3 getBlockIndexes(double blockX1Coord, double blockX2Coord, double blockX3Coord, int level)  const;
   UbTupleDouble3 getBlockLengths(Block3DPtr block) const;
   UbTupleDouble6 getBlockOversize() const ;
   void setCoordinateTransformator(CoordinateTransformation3DPtr trafo);
   const CoordinateTransformation3DPtr getCoordinateTransformator() const ;
   void setDeltaX(double dx);
   void setDeltaX(double worldUnit, double gridUnit);
   double getDeltaX(int level) const;
   double getDeltaX(Block3DPtr block) const;
   UbTupleDouble3 getNodeOffset(Block3DPtr block) const ;
   UbTupleDouble3 getNodeCoordinates(Block3DPtr block, int ix1, int ix2, int ix3) const;
   UbTupleInt3 getNodeIndexes(Block3DPtr block, double nodeX1Coord, double nodeX2Coord, double nodeX3Coord) const;
   void setBlockNX(int nx1, int nx2, int nx3);
   UbTupleInt3 getBlockNX() const;
   UbTupleDouble3 getBlockWorldCoordinates(Block3DPtr block) const;
   UbTupleDouble3 getBlockWorldCoordinates(int blockX1Index, int blockX2Index, int blockX3Index, int level) const;
   void setNX1(int nx1);
   void setNX2(int nx2);
   void setNX3(int nx3);
   int  getNX1() const;
   int  getNX2() const;
   int  getNX3() const;
   void calcStartCoordinatesAndDelta(Block3DPtr block, double& worldX1, double& worldX2, double& worldX3, double& deltaX);
   void calcStartCoordinatesWithOutOverlap(Block3DPtr block, double& worldX1, double& worldX2, double& worldX3);
   //////////////////////////////////////////////////////////////////////////
   //LBM
   //double getDeltaT(Block3DPtr) const;
   //////////////////////////////////////////////////////////////////////////
   void setTimeStep(double step);
   double getTimeStep() const;
protected:
   void checkLevel(int level);
   bool hasLevel(int level) const;

   void fillExtentWithBlocks( UbTupleInt3 minInd, UbTupleInt3 maxInd );

   void getSubBlocksZero(int ix1, int ix2, int ix3, int level,std::vector<Block3DPtr>& blockVector, int levelDepth);

   void getSubBlocksEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksNorth(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksSouth(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksTop(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottom(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);

   void getSubBlocksSouthEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksSouthWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksNorthEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksNorthWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);

   void getSubBlocksTopEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksTopWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksTopNorth(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksTopSouth(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);

   void getSubBlocksBottomEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottomWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottomNorth(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottomSouth(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);

   void getSubBlocksTopNorthEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksTopNorthWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksTopSouthEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksTopSouthWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottomNorthEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottomNorthWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottomSouthEast(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);
   void getSubBlocksBottomSouthWest(int ix1, int ix2, int ix3, int level, std::vector<Block3DPtr>& blockVector, int levelDepth);

private:
   LevelSet levelSet;
   BlockIDMap blockIdMap;
   Interactor3DSet interactors;
   //GridObserversSet observers;

   signal_t sig;
   int rank;
   int bundle;
   
   bool periodicX1;
   bool periodicX2;
   bool periodicX3;

   int blockNx1;    
   int blockNx2;    
   int blockNx3; 

   int nx1;    
   int nx2;    
   int nx3;    

   CoordinateTransformation3DPtr trafo;
   double orgDeltaX;

   double timeStep;
   
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & nx1;
      ar & nx2;
      ar & nx3;
      ar & orgDeltaX;
      ar & trafo;
      ar & blockNx1;
      ar & blockNx2;
      ar & blockNx3;
      ar & rank;
      ar & bundle;
      ar & periodicX1;
      ar & periodicX2;
      ar & periodicX3;
      ar & levelSet;
      ar & blockIdMap;
      ar & timeStep;
      ar & interactors;
   }
};

#endif 

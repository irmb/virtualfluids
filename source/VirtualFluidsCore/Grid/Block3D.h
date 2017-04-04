#ifndef BLOCK3D_H
#define BLOCK3D_H

#include <sstream>
#include <iostream>

#include "UbMath.h"
#include "CoordinateTransformation3D.h"

#include "Block3DConnector.h"

#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/foreach.hpp>

#include <boost/shared_ptr.hpp>
class Block3D;
typedef boost::shared_ptr<Block3D> Block3DPtr;

#include "LBMKernel.h"

class Block3D
{
public:
   Block3D();
   Block3D(int x1, int x2, int x3, int level);
   virtual ~Block3D();
   bool operator==(const Block3D& src) const;
   bool operator!=(const Block3D& src) const;

   int getX1() const;
   int getX2() const;
   int getX3() const;

   void setActive(bool active);
   bool isActive()    const;
   bool isNotActive() const;

   void setKernel(LBMKernelPtr kernel);
   LBMKernelPtr getKernel() const;
   void deleteKernel();

   void setBundle(int bundle);
   int  getBundle() const;

   void setRank(int rank);
   int  getRank() const;

   void setLocalRank(int rank);
   int  getLocalRank() const;

   int  getGlobalID() const;
   void setGlobalID(int id);

   int  getLocalID() const;
   void setLocalID(int id);

   int  getPart() const;
   void setPart(int part);

   int  getLevel() const;
   void setLevel(int level);

   //Connector-Section
   void                 setConnector(Block3DConnectorPtr connector);
   Block3DConnectorPtr  getConnector(int dir) const;
   bool                 hasConnectors();
   void                 deleteConnectors();
   void pushBackSameLevelConnectors(  std::vector<Block3DConnectorPtr>& localSameLevelConnectors
                                    , std::vector<Block3DConnectorPtr>& remoteSameLevelConnectors );
   void pushBackLocalSameLevelConnectors( std::vector<Block3DConnectorPtr>& localSameLevelConnectors );
   void pushBackRemoteSameLevelConnectors( std::vector<Block3DConnectorPtr>& remoteSameLevelConnectors );
   void pushBackLocalInterpolationConnectorsCF( std::vector<Block3DConnectorPtr>& localInterpolationConnectors );
   void pushBackRemoteInterpolationConnectorsCF( std::vector<Block3DConnectorPtr>& remoteInterpolationConnectors );
   void pushBackLocalInterpolationConnectorsFC( std::vector<Block3DConnectorPtr>& localInterpolationConnectors );
   void pushBackRemoteInterpolationConnectorsFC( std::vector<Block3DConnectorPtr>& remoteInterpolationConnectors );
   void pushBackLocalSameLevelConnectors( std::vector<Block3DConnectorPtr>& localSameLevelConnectors, const int& dir);
   void pushBackRemoteSameLevelConnectors( std::vector<Block3DConnectorPtr>& remoteSameLevelConnectors, const int& dir );
   int getNumberOfLocalConnectors();
   int getNumberOfRemoteConnectors();
   int getNumberOfLocalConnectorsForSurfaces();
   int getNumberOfRemoteConnectorsForSurfaces();

   void setWeight(int rank, int weight);
   int  getWeight(int rank);
   void addWeightForAll(int weight);
   void addWeight(int rank, int weight);
   void clearWeight();
   int  getWeightSize();

   //interpolation
   bool hasInterpolationFlag();
   bool hasInterpolationFlag(int dir);
   void deleteInterpolationFlag();

   void setInterpolationFlagCF(int dir);
   int  getInterpolationFlagCF();
   bool hasInterpolationFlagCF(int dir);
   bool hasInterpolationFlagCF();

   void setInterpolationFlagFC(int dir);
   int  getInterpolationFlagFC();
   bool hasInterpolationFlagFC(int dir);
   bool hasInterpolationFlagFC();

   double getWorkLoad();

   std::string toString() ;

   static int getMaxGlobalID() { return counter; }
private:
  int   x1;
  int   x2;
  int   x3;

  bool active;

  int interpolationFlagCF;
  int interpolationFlagFC;

  LBMKernelPtr kernel;
  std::vector<Block3DConnectorPtr> connectors;
  std::map<int, int> weight;

  int bundle;
  int rank;
  int lrank;
  int globalID;
  int localID;
  int part;
  int level;
  static int counter;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & x1;
     ar & x2;
     ar & x3;
     ar & active;
     ar & bundle;
     ar & rank;
     ar & lrank;
     ar & part;
     ar & globalID;
     ar & localID;
     ar & level;
     ar & kernel;
     ar & interpolationFlagCF;
     ar & interpolationFlagFC;
     ar & counter;
  }
};

#endif  //BLOCK3D_H

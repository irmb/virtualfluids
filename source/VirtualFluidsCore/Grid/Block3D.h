#ifndef BLOCK3D_H
#define BLOCK3D_H

#include <PointerDefinitions.h>
#include <vector>
#include <map>

class Block3DConnector;
class LBMKernel;
class ILBMKernel;

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

   void setKernel(SPtr<LBMKernel> kernel);
   SPtr<ILBMKernel> getKernel() const;
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
   void                 setConnector(SPtr<Block3DConnector> connector);
   SPtr<Block3DConnector>  getConnector(int dir) const;
   bool                 hasConnectors();
   void                 deleteConnectors();
   void pushBackSameLevelConnectors(  std::vector<SPtr<Block3DConnector> >& localSameLevelConnectors
                                    , std::vector<SPtr<Block3DConnector> >& remoteSameLevelConnectors );
   void pushBackLocalSameLevelConnectors( std::vector<SPtr<Block3DConnector> >& localSameLevelConnectors );
   void pushBackRemoteSameLevelConnectors( std::vector<SPtr<Block3DConnector> >& remoteSameLevelConnectors );
   void pushBackLocalInterpolationConnectorsCF( std::vector<SPtr<Block3DConnector> >& localInterpolationConnectors );
   void pushBackRemoteInterpolationConnectorsCF( std::vector<SPtr<Block3DConnector> >& remoteInterpolationConnectors );
   void pushBackLocalInterpolationConnectorsFC( std::vector<SPtr<Block3DConnector> >& localInterpolationConnectors );
   void pushBackRemoteInterpolationConnectorsFC( std::vector<SPtr<Block3DConnector> >& remoteInterpolationConnectors );
   void pushBackLocalSameLevelConnectors( std::vector<SPtr<Block3DConnector> >& localSameLevelConnectors, const int& dir);
   void pushBackRemoteSameLevelConnectors( std::vector<SPtr<Block3DConnector> >& remoteSameLevelConnectors, const int& dir );
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

   int  getCollectionOfInterpolationFlagCF();
   void setCollectionOfInterpolationFlagCF(int flags);
   
   void setInterpolationFlagCF(int dir);
   bool hasInterpolationFlagCF(int dir);
   bool hasInterpolationFlagCF();

   int  getCollectionOfInterpolationFlagFC();
   void setCollectionOfInterpolationFlagFC(int flags);
   
   void setInterpolationFlagFC(int dir);
   bool hasInterpolationFlagFC(int dir);
   bool hasInterpolationFlagFC();

   double getWorkLoad();

   std::string toString() ;

   static int getMaxGlobalID() { return counter; }
   static void setMaxGlobalID(int c) { counter = 0; }

private:
  int   x1;
  int   x2;
  int   x3;

  bool active;

  int interpolationFlagCF;
  int interpolationFlagFC;

  SPtr<LBMKernel> kernel;
  std::vector<SPtr<Block3DConnector> > connectors;
  std::map<int, int> weight;

  int bundle;
  int rank;
  int lrank;
  int globalID;
  int localID;
  int part;
  int level;
  static int counter;

};

#endif  //BLOCK3D_H

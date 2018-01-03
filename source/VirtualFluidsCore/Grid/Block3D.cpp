#include "Block3D.h"

#include "Grid3DSystem.h"
#include "Block3DConnector.h"
#include "LBMKernel.h"


int Block3D::counter = 0;
//////////////////////////////////////////////////////////////////////////
Block3D::Block3D() : x1(0),x2(0),x3(0)
                     ,active(true)
                     ,globalID(-1)
                     ,rank(-1),part(-1)
                     ,interpolationFlagCF(0)
                     ,interpolationFlagFC(0)
                     ,level(-1)
                     ,bundle(-1)
{
}
//////////////////////////////////////////////////////////////////////////
Block3D::Block3D(int x1, int x2, int x3, int level)
               : x1(x1), x2(x2), x3(x3)
               ,active(true)
               ,globalID(-1)
               ,rank(0),part(0)
               ,interpolationFlagCF(0)
               ,interpolationFlagFC(0)
               ,level(level)
               ,bundle(0)
{
   globalID = counter++;
}
//////////////////////////////////////////////////////////////////////////
Block3D::~Block3D()
{
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::operator==(const Block3D& src) const
{
   return (x1==src.x1 && x2==src.x2 && x3==src.x3); 
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::operator!=(const Block3D& src) const
{
   return !(*this==src);
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getX1() const 
{ 
   return this->x1; 
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getX2() const 
{ 
   return this->x2; 
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getX3() const 
{ 
   return this->x3; 
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setActive(bool active) 
{ 
   this->active = active; 
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::isActive()    const           
{ 
   return this->active;   
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::isNotActive() const           
{ 
   return(!this->active); 
}
//////////////////////////////////////////////////////////////////////////
void  Block3D::setKernel(LBMKernelPtr kernel) 
{  
   this->kernel = kernel; 
}
//////////////////////////////////////////////////////////////////////////
ILBMKernelPtr Block3D::getKernel() const              
{  
   return this->kernel; 
}
//////////////////////////////////////////////////////////////////////////
void Block3D::deleteKernel()             
{  
   this->kernel = LBMKernelPtr(); 
}
//////////////////////////////////////////////////////////////////////////
int  Block3D::getBundle() const          
{ 
   return bundle;       
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setBundle(int bundle) 
{ 
   this->bundle = bundle; 
} 
//////////////////////////////////////////////////////////////////////////
void  Block3D::setRank(int rank) 
{  
   this->rank = rank; 
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getRank() const            
{  
   return this->rank; 
}
//////////////////////////////////////////////////////////////////////////
void  Block3D::setLocalRank(int rank) 
{  
   this->lrank = rank; 
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getLocalRank() const            
{  
   return this->lrank; 
}
//////////////////////////////////////////////////////////////////////////
int  Block3D::getGlobalID() const        
{ 
   return this->globalID; 
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setGlobalID(int id) 
{ 
   this->globalID = id; 
}
//////////////////////////////////////////////////////////////////////////
int  Block3D::getLocalID() const        
{ 
   return this->localID; 
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setLocalID(int id) 
{ 
   this->localID = id; 
}
//////////////////////////////////////////////////////////////////////////
int  Block3D::getPart() const        
{ 
   return this->part; 
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setPart(int part) 
{ 
   this->part = part; 
}
//////////////////////////////////////////////////////////////////////////
int  Block3D::getLevel() const        
{ 
   return this->level; 
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setLevel(int level) 
{ 
   this->level = level; 
}
//////////////////////////////////////////////////////////////////////////
Block3DConnectorPtr Block3D::getConnector(int dir) const
{ 
   for(Block3DConnectorPtr c : connectors)
   {
      if( c ) 
      {
            if(c->getSendDir() == dir) return c;
      }
   }
  return Block3DConnectorPtr();     
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setConnector(Block3DConnectorPtr connector)
{
   connectors.push_back(connector);
}
//////////////////////////////////////////////////////////////////////////
void Block3D::deleteConnectors()
{
   connectors.clear();
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::hasConnectors()
{
   for(Block3DConnectorPtr c : connectors)
      if( c ) return true;
   
   return false;
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackSameLevelConnectors(  std::vector<Block3DConnectorPtr>& localSameLevelConnectors
                                            , std::vector<Block3DConnectorPtr>& remoteSameLevelConnectors )
{
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isLocalConnector() && !connector->isInterpolationConnectorCF() && !connector->isInterpolationConnectorFC() ) 
            localSameLevelConnectors.push_back(this->connectors[i]);
         else                                
            remoteSameLevelConnectors.push_back(this->connectors[i]);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackLocalSameLevelConnectors( std::vector<Block3DConnectorPtr>& localSameLevelConnectors )
{
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isLocalConnector() && !connector->isInterpolationConnectorCF() && !connector->isInterpolationConnectorFC() ) 
            localSameLevelConnectors.push_back(this->connectors[i]);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackLocalSameLevelConnectors( std::vector<Block3DConnectorPtr>& localSameLevelConnectors, const int& dir)
{
   Block3DConnectorPtr connector = this->connectors[dir];
   if( this->connectors[dir] )
   {
      if( connector->isLocalConnector() && !connector->isInterpolationConnectorCF() && !connector->isInterpolationConnectorFC() ) 
         localSameLevelConnectors.push_back(this->connectors[dir]);
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackRemoteSameLevelConnectors( std::vector<Block3DConnectorPtr>& remoteSameLevelConnectors )
{
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isRemoteConnector() && !connector->isInterpolationConnectorCF() && !connector->isInterpolationConnectorFC() ) 
            remoteSameLevelConnectors.push_back(this->connectors[i]);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackRemoteSameLevelConnectors( std::vector<Block3DConnectorPtr>& remoteSameLevelConnectors, const int& dir )
{
   Block3DConnectorPtr connector = this->connectors[dir];
   if( this->connectors[dir] )
   {
      if( connector->isRemoteConnector() && !connector->isInterpolationConnectorCF() && !connector->isInterpolationConnectorFC() ) 
         remoteSameLevelConnectors.push_back(this->connectors[dir]);
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackLocalInterpolationConnectorsCF( std::vector<Block3DConnectorPtr>& localInterpolationConnectors )
{
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isLocalConnector() && connector->isInterpolationConnectorCF() )
            localInterpolationConnectors.push_back(this->connectors[i]);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackRemoteInterpolationConnectorsCF( std::vector<Block3DConnectorPtr>& remoteInterpolationConnectors )
{
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isRemoteConnector() && connector->isInterpolationConnectorCF() )
            remoteInterpolationConnectors.push_back(this->connectors[i]);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackLocalInterpolationConnectorsFC( std::vector<Block3DConnectorPtr>& localInterpolationConnectors )
{
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isLocalConnector() && connector->isInterpolationConnectorFC() )
            localInterpolationConnectors.push_back(this->connectors[i]);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::pushBackRemoteInterpolationConnectorsFC( std::vector<Block3DConnectorPtr>& remoteInterpolationConnectors )
{
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isRemoteConnector() && connector->isInterpolationConnectorFC() )
            remoteInterpolationConnectors.push_back(this->connectors[i]);
      }
   }
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getNumberOfLocalConnectors()
{
   int count = 0;
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isLocalConnector() ) count++;
      }
   }
   return count;
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getNumberOfRemoteConnectors()
{
   int count = 0;
   for(int i=0; i<(int)connectors.size(); i++)
   {
      Block3DConnectorPtr connector = this->connectors[i];
      if( this->connectors[i] )
      {
         if( connector->isRemoteConnector() ) count++;
      }
   }
   return count;
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getNumberOfLocalConnectorsForSurfaces()
{
   int count = 0;
   
   if(connectors.size() < 6)
      return count;

   for(int dir=0; dir<=5; dir++) //Hard coding. It works if you have 0...5 for E, N ... B 
   {
      Block3DConnectorPtr connector = this->connectors[dir];
      if( this->connectors[dir] )
      {
         if( connector->isLocalConnector() ) count++;
      }
   }
   return count;
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getNumberOfRemoteConnectorsForSurfaces()
{
   int count = 0;
   for(int dir=0; dir<=5; dir++) //Hard coding. It works if you have 0...5 for E, N ... B 
   {
      Block3DConnectorPtr connector = this->connectors[dir];
      if( this->connectors[dir] )
      {
         if( connector->isRemoteConnector() ) count++;
      }
   }
   return count;
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setInterpolationFlagCF(int dir)
{
   UbSystem::setBit(interpolationFlagCF, 1<<dir);
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getInterpolationFlagCF()
{
   return interpolationFlagCF;
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::hasInterpolationFlagCF(int dir)
{
   return UbSystem::bitCheck( interpolationFlagCF, 1<<dir );
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setInterpolationFlagFC(int dir)
{
   UbSystem::setBit(interpolationFlagFC, 1<<dir);
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getInterpolationFlagFC()
{
   return interpolationFlagFC;
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::hasInterpolationFlagFC(int dir)
{
   return UbSystem::bitCheck( interpolationFlagFC, 1<<dir );
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::hasInterpolationFlag()
{ 
   return(interpolationFlagCF!=0 || interpolationFlagFC!=0); 
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::hasInterpolationFlag(int direction)     
{ 
   return(hasInterpolationFlagCF(direction) || hasInterpolationFlagFC(direction)); 
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::hasInterpolationFlagCF()
{
   return(interpolationFlagCF!=0);
}
//////////////////////////////////////////////////////////////////////////
bool Block3D::hasInterpolationFlagFC()
{
   return(interpolationFlagFC!=0);
}
//////////////////////////////////////////////////////////////////////////
void Block3D::deleteInterpolationFlag()
{
   interpolationFlagFC = 0;
   interpolationFlagCF = 0;
}
//////////////////////////////////////////////////////////////////////////
double Block3D::getWorkLoad()
{
   double l = kernel->getCalculationTime();
   l *= static_cast<double>(1<<level);
   return l;
}
//////////////////////////////////////////////////////////////////////////
std::string Block3D::toString() 
{
   std::stringstream ss;
   ss<<"Block3D[(x1,x2,x3,level),";
   ss<<" ("<<this->x1<<", "<<this->x2<<", "<<this->x3<<", "<<this->level<<"), id=" << globalID; 
   ss<< ", active="<<this->active<< ", bundle="<<this->bundle<< ", rank="<<this->rank<<"]";
   ss<<" connectors:";
   for(std::size_t i=0; i<connectors.size(); i++)
      if( connectors[i] )
      {
         if(connectors[i]->isLocalConnector())
            ss <<"l."<< Grid3DSystem::getDirectionString(connectors[i]->getSendDir()) << ", ";
         if(connectors[i]->isRemoteConnector())
            ss <<"r."<< Grid3DSystem::getDirectionString(connectors[i]->getSendDir()) << ", ";
         if(connectors[i]->isInterpolationConnectorCF())
            ss <<"cf."<< Grid3DSystem::getDirectionString(connectors[i]->getSendDir()) << ", ";
         if(connectors[i]->isInterpolationConnectorFC())
            ss <<"fc."<< Grid3DSystem::getDirectionString(connectors[i]->getSendDir()) << ", ";
      }
   return ss.str();
}
//////////////////////////////////////////////////////////////////////////
void Block3D::setWeight( int rank, int weight )
{
   std::map<int, int>::iterator it;
   if((it = this->weight.find(rank)) != this->weight.end())
       it->second = weight;
   else
      this->weight.insert(std::make_pair(rank, weight));
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getWeight( int rank )
{
   std::map<int, int>::iterator it;
   if((it = this->weight.find(rank)) != this->weight.end())
      return it->second;
   else
      return 0;
}
//////////////////////////////////////////////////////////////////////////
void Block3D::addWeight( int rank, int weight )
{
   int weight_old = getWeight(rank);
   weight += weight_old;
   setWeight(rank, weight);
}
//////////////////////////////////////////////////////////////////////////
void Block3D::addWeightForAll( int weight )
{
   typedef std::map<int, int> wMap;
   for (wMap::value_type &w : this->weight)
   {
      w.second += weight;
   }
}
//////////////////////////////////////////////////////////////////////////
void Block3D::clearWeight()
{
   this->weight.clear();
}
//////////////////////////////////////////////////////////////////////////
int Block3D::getWeightSize()
{
   return static_cast<int>(this->weight.size());
}

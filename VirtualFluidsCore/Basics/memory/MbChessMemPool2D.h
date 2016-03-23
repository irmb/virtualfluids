//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef MBCHESSMEMPOOL2D_H
#define MBCHESSMEMPOOL2D_H

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <typeinfo>

#include <basics/utilities/UbException.h>


template <class TData, std::size_t cachSize>
class MbChessMemPool2D;

//////////////////////////////////////////////////////////////////////////
//class MbChessMap2DKey
//key zum Auffinden der ChessMem2DBlocks
class MbChessMap2DKey
{
public:
   //////////////////////////////////////////////////////////////////////////
   //Konstruktoren
   MbChessMap2DKey(): mVectorPos(0),mFirstCriteria(0),mSecondCriteria(0)
   {

   }
   /*==========================================================*/
   MbChessMap2DKey(std::size_t vectorPos, std::size_t firstCriteria, std::size_t secondCriteria)
      : mVectorPos(vectorPos), mFirstCriteria(firstCriteria), mSecondCriteria(secondCriteria)
   {
   }
   /*==========================================================*/
   MbChessMap2DKey& operator=(const MbChessMap2DKey& srcKey)
   {
      if(this == &srcKey ) return *this;

      mVectorPos      = srcKey.mVectorPos;
      mFirstCriteria  = srcKey.mFirstCriteria;
      mSecondCriteria = srcKey.mSecondCriteria;

      return *this;
   }

   //////////////////////////////////////////////////////////////////////////
   //global ueberladene Operatoren
   friend inline bool operator<(const MbChessMap2DKey& lhsKey,const MbChessMap2DKey& rhsKey)
   {
      if(lhsKey.mFirstCriteria  < rhsKey.mFirstCriteria ) return true;
      if(lhsKey.mFirstCriteria  > rhsKey.mFirstCriteria ) return false;
      if(lhsKey.mSecondCriteria < rhsKey.mSecondCriteria) return true;

      return false;
   }
   /*==========================================================*/
   friend inline bool operator==(const MbChessMap2DKey& lhsKey,const MbChessMap2DKey& rhsKey)
   {
      if(lhsKey.mVectorPos      != rhsKey.mVectorPos      ) return false;
      if(lhsKey.mFirstCriteria  != rhsKey.mFirstCriteria  ) return false;
      if(lhsKey.mSecondCriteria != rhsKey.mSecondCriteria ) return false;

      return true;
   }
   //ueberladene Operatoren
   friend inline bool operator!=(const MbChessMap2DKey& lhsKey,const MbChessMap2DKey& rhsKey)
   {
      return !(lhsKey==rhsKey);
   }
   //ueberladene Operatoren
   /*==========================================================*/
   friend inline std::ostream& operator << (std::ostream& os, const MbChessMap2DKey& key)
   {
      os<<"VectorPos,first-,second-,third Criteria) (";
      os<<key.mVectorPos<<","<<key.mFirstCriteria<<","<<key.mSecondCriteria<<")";
      return os;
   }

   //////////////////////////////////////////////////////////////////////////
   //public Methoden
   std::size_t getVectorPos() {return mVectorPos;}
private:
   //////////////////////////////////////////////////////////////////////////
   //private Member
   std::size_t mVectorPos;
   std::size_t mFirstCriteria;
   std::size_t mSecondCriteria;
};



template<class T,std::size_t cachSize>
class MbChessMem2DBlock
{
   friend class MbChessMemPool2D<T,cachSize>;
public:
   //////////////////////////////////////////////////////////////////////////
   //Konstruktoren
   MbChessMem2DBlock()
   {
      mUsedElements = 0;
      std::size_t arrayLength = mBlockWidth*mBlockWidth;
      //mDataElements = new T[arrayLength];
      mDataElements = operator new(arrayLength*sizeof(T));
      mFlagVector   = new bool[arrayLength];
      for(std::size_t i=0;i<arrayLength;i++) mFlagVector[i] = false;
   }
   //////////////////////////////////////////////////////////////////////////
   //Destruktor
   ~MbChessMem2DBlock()
   {
      //if(mDataElements) delete[] mDataElements;
      if(mDataElements) operator delete(mDataElements);
      if(mFlagVector)   delete[] mFlagVector;
   }

private:
   //////////////////////////////////////////////////////////////////////////
   //private Methoden
   void* getReference(std::size_t chessX1, std::size_t chessX2)
   {
      std::size_t arrayIndex = chessX2*mBlockWidth + chessX1;
      #ifdef _DEBUG
         if(arrayIndex>=mBlockWidth*mBlockWidth) UB_THROW( UbException(UB_EXARGS,"index out of range") );
      #endif

      if(mFlagVector[arrayIndex]==true) UB_THROW( UbException(UB_EXARGS,"memory already allocated!") );

      mUsedElements++;
      mFlagVector[arrayIndex]=true;

      return (void*)((T*)(mDataElements)+arrayIndex);//&(mDataElements[arrayIndex]);
   }
   /*==========================================================*/
   std::size_t freeReference(void* p)
   {
      //std::size_t arrayIndex = static_cast<std::size_t>(static_cast<T*>(p) - mDataElements); //mDataElements = &mDataElements[0]
      std::size_t arrayIndex = static_cast<std::size_t>(static_cast<T*>(p) - static_cast<T*>(mDataElements));

      #ifdef _DEBUG
         if(arrayIndex>=mBlockWidth*mBlockWidth) UB_THROW( UbException(UB_EXARGS,"index out of range") );
      #endif

      if(mFlagVector[arrayIndex]==false) UB_THROW( UbException(UB_EXARGS,"memory not allocated!") );

      mFlagVector[arrayIndex]=false;

      return --mUsedElements;
   }
   /*==========================================================*/
   std::size_t getNofUsedElements()   { return mUsedElements; }
   /*==========================================================*/
   void addPointerToTElementsToVector(std::vector<T*>& tdataVector)
   {
      std::size_t arrayLength = mBlockWidth*mBlockWidth;
      for(std::size_t arrayIndex=0;arrayIndex<arrayLength;arrayIndex++)
      {
         //if(mFlagVector[arrayIndex]) tdataVector.push_back(&mDataElements[arrayIndex]);
         if(mFlagVector[arrayIndex]) tdataVector.push_back(static_cast<T*>(mDataElements)+arrayIndex);
      }
   }
   /*==========================================================*/
   template<typename Pred>
   void addPointerToTElementsToVector(std::vector<T*>& tdataVector, Pred& pred)
   {
      std::size_t arrayLength = mBlockWidth*mBlockWidth;
      T* tmp;
      for(std::size_t arrayIndex=0;arrayIndex<arrayLength;arrayIndex++)
      {
         if(mFlagVector[arrayIndex])
         {
            //tmp = &mDataElements[arrayIndex];
            tmp = (static_cast<T*>(mDataElements))+arrayIndex;
            if( pred(*tmp) ) tdataVector.push_back(tmp);
         }
      }
   }
private:
   //////////////////////////////////////////////////////////////////////////
   //static Member
   static const std::size_t   mBlockWidth;

   //////////////////////////////////////////////////////////////////////////
   //private Member
   std::size_t   mUsedElements;
   //T*    mDataElements;
   void* mDataElements;
   bool* mFlagVector;

};

//////////////////////////////////////////////////////////////////////////
//class MbChessMemPool2D
//zum Verwalten von TData Elementen in einer Schabrett-artigen Struktur
//die ChessMemBloecke haben hier eine Groesse von ~cachSize
template <class TData, std::size_t cachSize>
class MbChessMemPool2D
{
private:
   //////////////////////////////////////////////////////////////////////////
   //protected static const Member
   const static std::size_t mCacheSize;

   //////////////////////////////////////////////////////////////////////////
   //protected Member
   static std::vector< std::map< MbChessMap2DKey , MbChessMem2DBlock< TData,cachSize >* > > mMapVector;
   static std::map< void*, MbChessMap2DKey > mPointerKeyMap;

   //////////////////////////////////////////////////////////////////////////
   //protected Konstrukoren
   MbChessMemPool2D() //protected, um max einmal vererbt werden zu koennen!!!
   {              //zudem kann man so keine elmente von TreeBasedMemPool erstellen

   }
   //////////////////////////////////////////////////////////////////////////
   //Destruktor
    ~MbChessMemPool2D()
   {
   }

public:

   //////////////////////////////////////////////////////////////////////////
   //static public Methoden
   static void* getReference(std::size_t level, std::size_t ix1, std::size_t ix2)
   {
      if(!MbChessMem2DBlock< TData,cachSize >::mBlockWidth)
      {
         std::stringstream ss;
         ss<<"TreeBasedMemPool() - InitialisationError\n";
         ss<<"\t size of StorageData ("<<typeid(TData).name()<<", "<<sizeof(TData)<<" byte)\n";
         ss<<"\t exceeds user-specifyed cache-zize ("<<mCacheSize<<" byte)\n";
         ss<<"\t cache-size has to be larger than data-size";
         UB_THROW( UbException(ss.str()) );
      }

      if( mMapVector.size()<=level ) mMapVector.resize(level+1);

      std::size_t chessX1 = ix1/(MbChessMem2DBlock<TData,cachSize>::mBlockWidth);
      std::size_t chessX2 = ix2/(MbChessMem2DBlock<TData,cachSize>::mBlockWidth);

      MbChessMap2DKey mapKey(level,chessX1,chessX2);

      typename std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* >::iterator pos = mMapVector[level].find(mapKey);

      MbChessMem2DBlock<TData,cachSize>* memBlock = NULL;

      if(pos==mMapVector[level].end())
      {
         memBlock = new MbChessMem2DBlock<TData,cachSize>;
         (mMapVector[level])[mapKey] = memBlock;
      }
      else memBlock = pos->second;

      std::size_t internalChessX1 = ix1%(MbChessMem2DBlock<TData,cachSize>::mBlockWidth);
      std::size_t internalChessX2 = ix2%(MbChessMem2DBlock<TData,cachSize>::mBlockWidth);

      void* p = memBlock->getReference(internalChessX1,internalChessX2);

      mPointerKeyMap[p]=mapKey;

      return p;
   }
   /*==========================================================*/
   static void freeReference(void *p)
   {
      typename std::map<void*,MbChessMap2DKey>::iterator posPointerKeyMap = mPointerKeyMap.find(p);

      if(posPointerKeyMap==mPointerKeyMap.end()) UB_THROW( UbException(UB_EXARGS,"pointer not in map") );

      MbChessMap2DKey mapKey = posPointerKeyMap->second;
      mPointerKeyMap.erase(posPointerKeyMap);

      typename std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* >::iterator posMemBlockMap;
      posMemBlockMap = mMapVector[mapKey.getVectorPos()].find(mapKey);

      if(posMemBlockMap == mMapVector[mapKey.getVectorPos()].end())
         UB_THROW( UbException(UB_EXARGS,"mapKey not in ChessMem2DBlockMap") );

      std::size_t leftElements = posMemBlockMap->second->freeReference(p);
      if(!leftElements)
      {
         MbChessMem2DBlock<TData,cachSize>* tmp = posMemBlockMap->second;
         mMapVector[mapKey.getVectorPos()].erase(posMemBlockMap);
         delete tmp;
      }
   }
   /*==========================================================*/
   static void fillVectorWithPointerToTDataElements(std::size_t level,std::vector<TData*>& tdataVector)
   {
      tdataVector.clear();

      if(level>=mMapVector.size()) return;
      typename std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* >::iterator pos;
      for(pos=mMapVector[level].begin();pos!=mMapVector[level].end();++pos)
      {
         pos->second->addPointerToTElementsToVector(tdataVector);
      }
   }
   /*==========================================================*/
   static void fillVectorWithPointerToTDataElements(std::vector<TData*>& tdataVector)
   {
      tdataVector.clear();

      typename std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* >::iterator pos;

      for(std::size_t vecIndex=0; vecIndex<mMapVector.size(); vecIndex++ )
      {
         for(pos=mMapVector[vecIndex].begin();pos!=mMapVector[vecIndex].end();++pos)
         {
            pos->second->addPointerToTElementsToVector(tdataVector);
         }
      }
   }
   /*==========================================================*/
   template<typename Pred>
   static void fillVectorWithPointerToTDataElements(std::size_t level,std::vector<TData*>& tdataVector, Pred pred)
   {
      tdataVector.clear();

      if(level>=mMapVector.size()) return;
      typename std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* >::iterator pos;
      for(pos=mMapVector[level].begin();pos!=mMapVector[level].end();++pos)
      {
         pos->second->addPointerToTElementsToVector(tdataVector,pred);
      }
   }
   /*==========================================================*/
   template<typename Pred>
   static void fillVectorWithPointerToTDataElements(std::vector<TData*>& tdataVector, Pred pred)
   {
      tdataVector.clear();

      typename std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* >::iterator pos;

      for(std::size_t vecIndex=0; vecIndex<mMapVector.size(); vecIndex++ )
      {
         for(pos=mMapVector[vecIndex].begin();pos!=mMapVector[vecIndex].end();++pos)
         {
            pos->second->addPointerToTElementsToVector(tdataVector,pred);
         }
      }
   }
   /*==========================================================*/
   static std::size_t getNumberOfChessMemoryBlocks()
   {
      std::size_t nofElements = 0;
      for(std::size_t i=0; i<mMapVector.size(); i++)
      {
         nofElements+=mMapVector[i].size();
      }
      return nofElements;
   }
   /*==========================================================*/
   static std::size_t getNumberOfChessMemoryBlocks(std::size_t level)
   {
      if(level<mMapVector.size() )return mMapVector[level].size();
      return 0;
   }
   /*==========================================================*/
   static std::size_t getNumberOfStoredDataElements()
   {
      return mPointerKeyMap.size();
   }
   /*==========================================================*/
   static std::size_t getNumberOfStoredDataElements(std::size_t level)
   {
      if(level<mMapVector.size() )
      {
         std::size_t nofElements = 0;
         typename std::map< MbChessMap2DKey , MbChessMem2DBlock< TData,cachSize >* >::iterator pos;

         for(pos=mMapVector[level].begin(); pos!=mMapVector[level].end(); ++pos)
         {
            nofElements += pos->second->getNofUsedElements();
         }
         return nofElements;
      }
      return 0;

   }
   /*==========================================================*/
   static std::string toString()
   {
      long double capaticityPerBlock   = (std::size_t)pow((double)MbChessMem2DBlock<TData,cachSize>::mBlockWidth,2.0);
      std::size_t storedElements       = MbChessMemPool2D<TData,cachSize>::getNumberOfStoredDataElements();
      std::size_t initialisedMemBlocks = MbChessMemPool2D<TData,cachSize>::getNumberOfChessMemoryBlocks();

      std::stringstream ss;
      ss<<std::endl;
      ss<<"****************** MbChessMemPool2D-Info (BEGIN) ******************"<<std::endl;
      ss<<"type of Storage-Data.................. : "<<typeid(TData).name()<<std::endl;
      ss<<"size of Storage-Data........... [bytes]: "<<sizeof(TData)<<std::endl;
      ss<<"specified cache-size........... [bytes]: "<<mCacheSize<<std::endl;
      ss<<"#elements per MbChessMem2DBlock [bytes]: "<<capaticityPerBlock<<std::endl;
      ss<<"mem per MbChessMem2DBlock...... [bytes]: "<<capaticityPerBlock*sizeof(TData)<<std::endl;
      ss<<"used cache-size[%]............. [bytes]: "<<capaticityPerBlock*sizeof(TData)/(double)mCacheSize*100<<std::endl;
      ss<<"\n";
      ss<<"#stored Elements  = "<<storedElements<<std::endl;
      ss<<"#ChessMem2DBlocks = "<<initialisedMemBlocks<<std::endl;
      ss<<std::endl;
      ss<<"level | #ChessMem2DBlocks | #stored Elements | used capaticity [%] \n";
      ss<<"----------------------------------------------------------------\n";
      for(std::size_t level=0;level<mMapVector.size();level++)
      {
         std::size_t nofStoredElements = MbChessMemPool2D<TData,cachSize>::getNumberOfStoredDataElements(level);
         std::size_t nofChessMem2DBlocks = MbChessMemPool2D<TData,cachSize>::getNumberOfChessMemoryBlocks(level);

         ss<<std::left<<" "<<std::setw(5)<<level<<"| "
            <<std::setw(16)<<nofChessMem2DBlocks<<"| "
            <<std::setw(17)<<nofStoredElements<<"| ";
         if(nofStoredElements)
            ss<<setw(15)<<nofStoredElements/(double)(capaticityPerBlock*nofChessMem2DBlocks)*100<<std::endl;
         else ss<<"-"<<std::endl;
      }
      ss<<std::endl;
      ss<<"called memory..... [bytes]: "<<storedElements*sizeof(TData)<<std::endl;
      ss<<"initialised memory [bytes]: "<<initialisedMemBlocks*capaticityPerBlock*sizeof(TData)<<std::endl;
      double denominator = (double)(initialisedMemBlocks*capaticityPerBlock*sizeof(TData));
      if(fabs(denominator)>1.E-13) ss<<"used.............. [%]    : "<<100.*storedElements*sizeof(TData)/denominator<<std::endl;
      else                         ss<<"used.............. [%]    : 0.0"<<std::endl;
      ss<<"****************** MbChessMemPool2D-Info (END)  *******************"<<std::endl;
      return ss.str();
   }
   /*==========================================================*/
   static void writeStatisticFiles(const std::string& filename)
   {
      //liefert Statistik ueber aufuellung der einzelnen bloecke (gesamt und pro level)
      //x-Achse: 0... max moegliche Anzahl von moeglichen Elementen pro MemBlock
      //y-Achse: Anzahl an Bloecken, die die Anzahl an Elementen beinhalten
      std::ofstream spreadingFile(((std::string)(filename+"_spreading.txt")).c_str());
      if(!spreadingFile) UB_THROW( UbException(UB_EXARGS,"couldn't open file") );

      std::size_t initialisedMemBlocks       =   MbChessMemPool2D<TData,cachSize>::getNumberOfChessMemoryBlocks();
      std::size_t maxNofDataElementsPerBlock =   MbChessMem2DBlock<TData,cachSize>::mBlockWidth
                                               * MbChessMem2DBlock<TData,cachSize>::mBlockWidth
                                               * MbChessMem2DBlock<TData,cachSize>::mBlockWidth;
      std::vector<std::size_t> spreading;
      spreading.resize(maxNofDataElementsPerBlock+1,0);
      std::vector< std::vector<std::size_t> > spreadingPerLevel;
      spreadingPerLevel.resize(mMapVector.size());

      typename std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* >::iterator pos;

      for(std::size_t level=0; level<mMapVector.size(); level++ )
      {
         spreadingPerLevel[level].resize(maxNofDataElementsPerBlock+1,0);
         for(pos=mMapVector[level].begin();pos!=mMapVector[level].end();++pos)
         {
            std::size_t number = pos->second->getNofUsedElements();
            spreading[number]++;
            spreadingPerLevel[level][number]++;
         }
      }
      spreadingFile<<"#BlockUsage nofBlocks(all Level) ";
      for(std::size_t level=0; level<mMapVector.size(); level++ )
         spreadingFile<<"nofBlockLevel"<<level<<" ";
      spreadingFile<<std::endl;

      for(std::size_t i=0; i<spreading.size(); i++)
      {
         spreadingFile<<i<<" "<<spreading[i];
         for(std::size_t level=0; level<mMapVector.size(); level++ )
            spreadingFile<<" "<<spreadingPerLevel[level][i];
         spreadingFile<<std::endl;
      }
      spreadingFile.flush();
      spreadingFile.close();
   }
   //////////////////////////////////////////////////////////////////////////
   //ueberladene operatoren
   void* operator new(size_t size, std::size_t level, std::size_t ix1, std::size_t ix2)
   {
      if(level<0) UB_THROW( UbException(UB_EXARGS,"level ist negativ!") );
      void *p = getReference(level,ix1,ix2);
      return p;
   }
   /*==========================================================*/
   void operator delete(void* p, std::size_t level, std::size_t ix1, std::size_t ix2)
   {
      //ACHTUNG: wenn man hier ne Exception schmeisst, dann gibts einen BoeSEN compilerFehler!!!
      //UB_THROW( UbException(UB_EXARGS,"Scheisse noch nicht gerafft, wie das geht!") );
      std::cerr<<"MbChessMemPool2D::delete(void* p, std::size_t level, std::size_t ix1, std::size_t ix2) - Scheisse noch nicht gerafft, wie das geht!\n";
   }

   /*==========================================================*/
   void operator delete(void* p)
   {
      freeReference(p);
   }

private:
   //////////////////////////////////////////////////////////////////////////
   //private statische Methoden
};

//statische Variablen initialisieren
template <class TData, std::size_t cachSize>
std::vector< std::map<MbChessMap2DKey,MbChessMem2DBlock<TData,cachSize>* > > MbChessMemPool2D<TData,cachSize>::mMapVector;

template <class TData, std::size_t cachSize>
std::map<void*,MbChessMap2DKey >  MbChessMemPool2D< TData, cachSize>::mPointerKeyMap;

template <class TData, std::size_t cachSize>
const std::size_t  MbChessMemPool2D<TData,cachSize>::mCacheSize=cachSize;

template <class TData,std::size_t cachSize>
const std::size_t  MbChessMem2DBlock<TData,cachSize>::mBlockWidth=static_cast<std::size_t>(pow(static_cast<double>(static_cast<std::size_t>(cachSize/sizeof(TData))),1./3.));

#endif

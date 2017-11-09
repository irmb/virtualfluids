//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef MBCHESSMEMPOOL3D_H
#define MBCHESSMEMPOOL3D_H

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include <basics/utilities/UbException.h>

template <class TData, std::size_t cachSize>
class MbChessMemPool3D;

//////////////////////////////////////////////////////////////////////////
//class MbChessMap3DKey
//key zum Auffinden der ChessMem3DBlocks
class MbChessMap3DKey
{
public:
   //////////////////////////////////////////////////////////////////////////
   //Konstruktoren
   MbChessMap3DKey(): mVectorPos(0),mFirstCriteria(0),mSecondCriteria(0),mThirdCriteria(0)
   {
   }
   /*==========================================================*/
   MbChessMap3DKey(std::size_t vectorPos,std::size_t firstCriteria,std::size_t secondCriteria,std::size_t thirdCriteria)
      : mVectorPos(vectorPos), mFirstCriteria(firstCriteria), mSecondCriteria(secondCriteria), mThirdCriteria(thirdCriteria)
   {
   }
   /*==========================================================*/
   MbChessMap3DKey& operator=(const MbChessMap3DKey& srcKey)
   {
      if(this == &srcKey ) return *this;

      mVectorPos      = srcKey.mVectorPos;
      mFirstCriteria  = srcKey.mFirstCriteria;
      mSecondCriteria = srcKey.mSecondCriteria;
      mThirdCriteria  = srcKey.mThirdCriteria;

      return *this;
   }

   //////////////////////////////////////////////////////////////////////////
   //global ueberladene Operatoren
   friend inline bool operator<(const MbChessMap3DKey& lhsKey,const MbChessMap3DKey& rhsKey)
   {
      if(lhsKey.mFirstCriteria  < rhsKey.mFirstCriteria ) return true;
      if(lhsKey.mFirstCriteria  > rhsKey.mFirstCriteria ) return false;
      if(lhsKey.mSecondCriteria < rhsKey.mSecondCriteria) return true;
      if(lhsKey.mSecondCriteria > rhsKey.mSecondCriteria) return false;
      if(lhsKey.mThirdCriteria  < rhsKey.mThirdCriteria ) return true;

      return false;
   }
   /*==========================================================*/
   friend inline bool operator==(const MbChessMap3DKey& lhsKey,const MbChessMap3DKey& rhsKey)
   {
      if(lhsKey.mVectorPos      != rhsKey.mVectorPos      ) return false;
      if(lhsKey.mFirstCriteria  != rhsKey.mFirstCriteria  ) return false;
      if(lhsKey.mSecondCriteria != rhsKey.mSecondCriteria ) return false;
      if(lhsKey.mThirdCriteria  != rhsKey.mThirdCriteria  ) return false;

      return true;
   }
   //ueberladene Operatoren
   friend inline bool operator!=(const MbChessMap3DKey& lhsKey,const MbChessMap3DKey& rhsKey)
   {
      return !(lhsKey==rhsKey);
   }
   //ueberladene Operatoren
   /*==========================================================*/
   friend inline std::ostream& operator << (std::ostream& os, const MbChessMap3DKey& key)
   {
      os<<"VectorPos,first-,second-,third Criteria) (";
      os<<key.mVectorPos<<","<<key.mFirstCriteria<<","<<key.mSecondCriteria<<","<<key.mThirdCriteria<<")";
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
   std::size_t mThirdCriteria;
};



template<class T,std::size_t cachSize>
class MbChessMem3DBlock
{
   friend class MbChessMemPool3D<T,cachSize>;
public:
   //////////////////////////////////////////////////////////////////////////
   //Konstruktoren
   MbChessMem3DBlock()
   {
      mUsedElements = 0;
      std::size_t arrayLength = mBlockWidth*mBlockWidth*mBlockWidth;
      //mDataElements = new T[arrayLength];
      mDataElements = operator new(arrayLength*sizeof(T));
      if(!mDataElements) throw UbException(UB_EXARGS,"out of memeory!");
      mFlagVector   = new bool[arrayLength];
      if(!mFlagVector) throw UbException(UB_EXARGS,"out of memeory!");
      for(std::size_t i=0;i<arrayLength;i++) mFlagVector[i] = false;
   }
   //////////////////////////////////////////////////////////////////////////
   //Destruktor
   ~MbChessMem3DBlock()
   {
      //if(mDataElements) delete[] mDataElements;
      if(mDataElements) operator delete(mDataElements);
      if(mFlagVector)   delete[] mFlagVector;
   }

private:
   //////////////////////////////////////////////////////////////////////////
   //private Methoden
   void* getReference(std::size_t chessX1, std::size_t chessX2, std::size_t chessX3)
   {
      //std::size_t arrayIndex = (chessX1*mBlockWidth+chessX2)*mBlockWidth + chessX3;
      std::size_t arrayIndex = (chessX3*mBlockWidth+chessX2)*mBlockWidth + chessX1;
      #ifdef _DEBUG
         if(arrayIndex>=mBlockWidth*mBlockWidth*mBlockWidth) throw UbException(UB_EXARGS,"index out of range");
      #endif

      if(mFlagVector[arrayIndex]==true) throw UbException(UB_EXARGS,"memory already allocated!");

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
        if(arrayIndex>=mBlockWidth*mBlockWidth*mBlockWidth) throw UbException(UB_EXARGS,"index out of range");
      #endif

      if(mFlagVector[arrayIndex]==false) throw UbException(UB_EXARGS,"memory not allocated!");

      mFlagVector[arrayIndex]=false;

      return --mUsedElements;
   }
   /*==========================================================*/
   std::size_t  getNofUsedElements()   { return mUsedElements; }
   /*==========================================================*/
   void addPointerToTElementsToVector(std::vector<T*>& tdataVector)
   {
      std::size_t arrayLength = mBlockWidth*mBlockWidth*mBlockWidth;
      for(std::size_t arrayIndex=0; arrayIndex<arrayLength; arrayIndex++)
      {
         //if(mFlagVector[arrayIndex]) tdataVector.push_back(&mDataElements[arrayIndex]);
         if(mFlagVector[arrayIndex]) tdataVector.push_back(static_cast<T*>(mDataElements)+arrayIndex);
      }
   }
   /*==========================================================*/
   template<typename Pred>
   void addPointerToTElementsToVector(std::vector<T*>& tdataVector,Pred& pred )
   {
      std::size_t arrayLength = mBlockWidth*mBlockWidth*mBlockWidth;
      T* tmp;
      for(std::size_t arrayIndex=0;arrayIndex<arrayLength;arrayIndex++)
      {
         if(mFlagVector[arrayIndex])
         {
            //tmp = &mDataElements[arrayIndex];
            tmp = static_cast<T*>(mDataElements)+arrayIndex;
            if( pred(*tmp) ) tdataVector.push_back(tmp);
         }
      }
   }
private:
   //////////////////////////////////////////////////////////////////////////
   //static Member
   static const std::size_t  mBlockWidth;

   //////////////////////////////////////////////////////////////////////////
   //private Member
   std::size_t mUsedElements;
   //T*    mDataElements;
   void* mDataElements;
   bool* mFlagVector;

};

//////////////////////////////////////////////////////////////////////////
//class MbChessMemPool3D
//zum Verwalten von TData Elementen in einer Schabrett-artigen Struktur
//die ChessMemBloecke haben hier eine Groesse von ~cachSize
template <class TData, std::size_t cachSize>
class MbChessMemPool3D
{
private:
   //////////////////////////////////////////////////////////////////////////
   //protected static const Member
   static const std::size_t mCacheSize;

   //////////////////////////////////////////////////////////////////////////
   //protected Member
   static std::vector< std::map< MbChessMap3DKey , MbChessMem3DBlock< TData,cachSize >* > > mMapVector;
   static std::map< void*, MbChessMap3DKey > mPointerKeyMap;

   //////////////////////////////////////////////////////////////////////////
   //protected Konstrukoren
   MbChessMemPool3D() //private, da NUR static erlaubt!!!
   {

   }
   //////////////////////////////////////////////////////////////////////////
   //Destruktor
   ~MbChessMemPool3D()
   {
   }

public:
   //////////////////////////////////////////////////////////////////////////
   //static public Methoden
   static void* getReference(std::size_t level, std::size_t ix1, std::size_t ix2, std::size_t ix3)
   {
      if(!MbChessMem3DBlock< TData,cachSize >::mBlockWidth)
      {
         std::stringstream ss;
         ss<<"TreeBasedMemPool() - InitialisationError\n";
         ss<<"\t size of StorageData ("<<typeid(TData).name()<<", "<<sizeof(TData)<<" byte)\n";
         ss<<"\t exceeds user-specifyed cache-zize ("<<mCacheSize<<" byte)\n";
         ss<<"\t cache-size has to be larger than data-size";
         throw UbException(UB_EXARGS,ss.str());
      }

      if( mMapVector.size()<=level ) mMapVector.resize(level+1);

      std::size_t chessX1 = ix1/(MbChessMem3DBlock<TData,cachSize>::mBlockWidth);
      std::size_t chessX2 = ix2/(MbChessMem3DBlock<TData,cachSize>::mBlockWidth);
      std::size_t chessX3 = ix3/(MbChessMem3DBlock<TData,cachSize>::mBlockWidth);

      MbChessMap3DKey mapKey(level,chessX1,chessX2,chessX3);

      typename std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* >::iterator pos = mMapVector[level].find(mapKey);

      MbChessMem3DBlock<TData,cachSize>* memBlock = NULL;

      if(pos==mMapVector[level].end())
      {
         memBlock = new MbChessMem3DBlock<TData,cachSize>;
         (mMapVector[level])[mapKey] = memBlock;
      }
      else memBlock = pos->second;

      std::size_t internalChessX1 = ix1%(MbChessMem3DBlock<TData,cachSize>::mBlockWidth);
      std::size_t internalChessX2 = ix2%(MbChessMem3DBlock<TData,cachSize>::mBlockWidth);
      std::size_t internalChessX3 = ix3%(MbChessMem3DBlock<TData,cachSize>::mBlockWidth);

      void* p = memBlock->getReference(internalChessX1,internalChessX2,internalChessX3);

      mPointerKeyMap[p]=mapKey;

      return p;
   }
   static void freeReference(void *p)
   {
      typename std::map<void*,MbChessMap3DKey>::iterator posPointerKeyMap = mPointerKeyMap.find(p);

      if(posPointerKeyMap==mPointerKeyMap.end()) throw UbException(UB_EXARGS,"pointer not in map");

      MbChessMap3DKey mapKey = posPointerKeyMap->second;
      mPointerKeyMap.erase(posPointerKeyMap);


      typename std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* >::iterator posMemBlockMap;
      posMemBlockMap = mMapVector[mapKey.getVectorPos()].find(mapKey);


      if(posMemBlockMap == mMapVector[mapKey.getVectorPos()].end())
         throw UbException(UB_EXARGS,"mapKey not in ChessMem3DBlockMap");

      std::size_t leftElements = posMemBlockMap->second->freeReference(p);
      if(!leftElements)
      {
         MbChessMem3DBlock<TData,cachSize>* tmp = posMemBlockMap->second;
         mMapVector[mapKey.getVectorPos()].erase(posMemBlockMap);
         try{ delete tmp; }
         catch(...){throw UbException(UB_EXARGS,"could not delete MbChessMem3DBlock");}
      }
   }
   /*==========================================================*/
   static void fillVectorWithPointerToTDataElements(std::size_t level,std::vector<TData*>& tdataVector)
   {
      tdataVector.clear();

      if(level>=mMapVector.size()) return;
      typename std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* >::iterator pos;
      for(pos=mMapVector[level].begin();pos!=mMapVector[level].end();++pos)
      {
         pos->second->addPointerToTElementsToVector(tdataVector);
      }
   }
   /*==========================================================*/
   static void fillVectorWithPointerToTDataElements(std::vector<TData*>& tdataVector)
   {
      tdataVector.clear();

      typename std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* >::iterator pos;

      for(std::size_t vecIndex=0; vecIndex<mMapVector.size(); vecIndex++ )
      {
         for(pos=mMapVector[vecIndex].begin();pos!=mMapVector[vecIndex].end();++pos)
         {
            pos->second->addPointerToTElementsToVector(tdataVector);
         }
      }
   }
   /*==========================================================*/
   template<class Pred>
   static void fillVectorWithPointerToTDataElements(std::size_t level,std::vector<TData*>& tdataVector, Pred pred)
   {
      tdataVector.clear();

      if(level>=mMapVector.size()) return;
      typename std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* >::iterator pos;
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

      typename std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* >::iterator pos;

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
      for(std::size_t i=0;i<mMapVector.size();i++)
      {
         nofElements+=mMapVector[i].size();
      }
      return nofElements;
   }
   /*==========================================================*/
   static std::size_t getNumberOfChessMemoryBlocks(std::size_t level)
   {
      if(level<mMapVector.size() ) return mMapVector[level].size();
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
         typename std::map< MbChessMap3DKey , MbChessMem3DBlock< TData,cachSize >* >::iterator pos;

         for(pos=mMapVector[level].begin(); pos!=mMapVector[level].end(); ++pos)
         {
            nofElements+= pos->second->getNofUsedElements();
         }
         return nofElements;
      }
      return 0;
   }
   /*==========================================================*/
   static std::string toString()
   {
      long double capaticityPerBlock   = pow((double)MbChessMem3DBlock<TData,cachSize>::mBlockWidth,3.0);
      std::size_t storedElements       = MbChessMemPool3D<TData,cachSize>::getNumberOfStoredDataElements();
      std::size_t initialisedMemBlocks = MbChessMemPool3D<TData,cachSize>::getNumberOfChessMemoryBlocks();

      std::stringstream ss;
      ss<<std::endl;
      ss<<"****************** MbChessMemPool3D-Info (BEGIN) ******************"<<std::endl;
      ss<<"type of Storage-Data.................. : "<<typeid(TData).name()<<std::endl;
      ss<<"size of Storage-Data........... [bytes]: "<<sizeof(TData)<<std::endl;
      ss<<"specified cache-size........... [bytes]: "<<mCacheSize<<std::endl;
      ss<<"#elements per MbChessMem3DBlock [bytes]: "<<capaticityPerBlock<<std::endl;
      ss<<"mem per MbChessMem3DBlock...... [bytes]: "<<capaticityPerBlock*sizeof(TData)<<std::endl;
      ss<<"used cache-size[%]............. [bytes]: "<<capaticityPerBlock*sizeof(TData)/(double)mCacheSize*100<<std::endl;
      ss<<"\n";
      ss<<"#stored Elements   = "<<storedElements<<std::endl;
      ss<<"#ChessMem3DBlocks  = "<<initialisedMemBlocks<<std::endl;
      ss<<std::endl;
      ss<<"level | #ChessMem3DBlocks | #stored Elements | used capaticity [%] \n";
      ss<<"----------------------------------------------------------------\n";
      for(std::size_t level=0;level<mMapVector.size();level++)
      {
         std::size_t nofStoredElements   = MbChessMemPool3D<TData,cachSize>::getNumberOfStoredDataElements(level);
         std::size_t nofChessMem3DBlocks = MbChessMemPool3D<TData,cachSize>::getNumberOfChessMemoryBlocks(level);

         ss<<std::left<<" "<<std::setw(5)<<level<<"| "
            <<std::setw(16)<<nofChessMem3DBlocks<<"| "
            <<std::setw(17)<<nofStoredElements<<"| ";
         if(nofStoredElements)
            ss<<std::setw(15)<<nofStoredElements/(double)(capaticityPerBlock*nofChessMem3DBlocks)*100<<std::endl;
         else ss<<"-"<<std::endl;
      }
      ss<<std::endl;
      ss<<"called memory..... [bytes]: "<<storedElements*sizeof(TData)<<std::endl;
      ss<<"initialised memory [bytes]: "<<initialisedMemBlocks*capaticityPerBlock*sizeof(TData)<<std::endl;
      double denominator = (double)(initialisedMemBlocks*capaticityPerBlock*sizeof(TData));
      if(fabs(denominator)>1.E-13) ss<<"used.............. [%]    : "<<100.*storedElements*sizeof(TData)/denominator<<std::endl;
      else                         ss<<"used.............. [%]    : 0.0"<<std::endl;
      ss<<"****************** MbChessMemPool3D-Info (END)  *******************"<<std::endl;
      return ss.str();
   }
   /*==========================================================*/
   static void writeStatisticFiles(const std::string& filename)
   {
      //liefert Statistik ueber aufuellung der einzelnen bloecke (gesamt und pro level)
      //x-Achse: 0... max moegliche Anzahl von moeglichen Elementen pro MemBlock
      //y-Achse: Anzahl an Bloecken, die die Anzahl an Elementen beinhalten
      std::ofstream spreadingFile(((std::string)(filename+"_spreading.txt")).c_str());
      if(!spreadingFile) throw UbException(UB_EXARGS,"couldn't open file");

      //std::size_t initialisedMemBlocks       =  MbChessMemPool3D<TData,cachSize>::getNumberOfChessMemoryBlocks();
      std::size_t maxNofDataElementsPerBlock =  MbChessMem3DBlock<TData,cachSize>::mBlockWidth
                                               *MbChessMem3DBlock<TData,cachSize>::mBlockWidth
                                               *MbChessMem3DBlock<TData,cachSize>::mBlockWidth;
      std::vector<std::size_t> spreading;
      spreading.resize(maxNofDataElementsPerBlock+1,0);
      std::vector< std::vector<std::size_t> > spreadingPerLevel;
      spreadingPerLevel.resize(mMapVector.size());

      typename std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* >::iterator pos;

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

      for(std::size_t i=0;i<spreading.size();i++)
      {
         spreadingFile<<i<<" "<<spreading[i];
         for(std::size_t level=0; level<mMapVector.size(); level++ )
            spreadingFile<<" "<<spreadingPerLevel[level][i];
         spreadingFile<<std::endl;
      }
      spreadingFile.flush();
      spreadingFile.close();
   }

   ////////////////////////////////////////////////////////////////////////////
   ////ueberladene operatoren
   //void* operator new(size_t size, std::size_t level, std::size_t ix1, std::size_t ix2, std::size_t ix3)
   //{
   //   if(level<0) throw UbException(UB_EXARGS,"level ist negativ!");
   //   void *p = getReference(level,ix1,ix2,ix3);
   //   return p;
   //}
   ///*==========================================================*/
   //void operator delete(void* p, std::size_t level, std::size_t ix1, std::size_t ix2, std::size_t ix3)
   //{
   //   //ACHTUNG: wenn man hier ne Exception schmeisst, dann gibts einen BoeSEN compilerFehler!!!
   //   //throw UbException(__FILE__, __LINE__, "MbChessMemPool3D::delete - Scheisse noch nicht gerafft, wie das geht!");
   //   cout<<"MbChessMemPool3D::delete(void* p, std::size_t level, std::size_t ix1, std::size_t ix2, std::size_t ix3) - Scheisse noch nicht gerafft, wie das geht!\n";
   //}

   ///*==========================================================*/
   //void operator delete(void* p)
   //{
   //   freeReference(p);
   //}

private:
   //////////////////////////////////////////////////////////////////////////
   //private statische Methoden
};


//statische Variablen initialisieren
template <class TData, std::size_t cachSize>
std::vector< std::map<MbChessMap3DKey,MbChessMem3DBlock<TData,cachSize>* > > MbChessMemPool3D<TData,cachSize>::mMapVector;

template <class TData, std::size_t cachSize>
std::map<void*,MbChessMap3DKey >  MbChessMemPool3D< TData, cachSize>::mPointerKeyMap;

template <class TData, std::size_t cachSize>
const std::size_t  MbChessMemPool3D<TData,cachSize>::mCacheSize=cachSize;

//template <class TData, std::size_t cachSize>
//const std::size_t  MbChessMemPool3D<TData,cachSize>::mNofElementsWidthMemBlock=static_cast<std::size_t>(pow(static_cast<double>(static_cast<std::size_t>(cachSize/sizeof(TData))),1./3.));

//template <class TData, std::size_t cachSize>
//const std::size_t  MbChessMemPool3D<TData,cachSize>::mNofElementsInMemBlock=static_cast<std::size_t>(pow(static_cast<double>(static_cast<std::size_t>(pow(static_cast<double>(static_cast<std::size_t>(cachSize/sizeof(TData))),1./3.))),3.0));

template <class TData,std::size_t cachSize>
const std::size_t  MbChessMem3DBlock<TData,cachSize>::mBlockWidth=static_cast<std::size_t>(std::pow(static_cast<double>(static_cast<std::size_t>(cachSize/sizeof(TData))),1./3.));

//template <class TData,std::size_t cachSize>
//const std::size_t  MbChessMem3DBlock<TData,cachSize>::mMaxElements=static_cast<std::size_t>(pow(static_cast<double>(static_cast<std::size_t>(pow(static_cast<double>(static_cast<std::size_t>(pow(static_cast<double>(static_cast<std::size_t>(cachSize/sizeof(TData))),1.0/3.0))),3.0))),3.0));

#endif

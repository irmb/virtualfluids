//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef BCARRAY_H
#define BCARRAY_H

#include <basics/container/CbArray3D.h>
#include <typeinfo>

#include <boost/serialization/serialization.hpp>

//IndexType must be int or long or short!!!
//no unsiged!

template< typename BCClass , class IndexType = int >
class BCArray3D
{
public:
   typedef typename boost::shared_ptr<BCClass> BCClassPtr;
   typedef std::size_t size_type;

public:
   //////////////////////////////////////////////////////////////////////////
   BCArray3D() {}
   //////////////////////////////////////////////////////////////////////////
   BCArray3D(const size_type& nx1, const size_type& nx2, const size_type& nx3)
   {
      bcindexmatrix.resize(nx1,nx2,nx3,UNDEFINED);
   }
   //////////////////////////////////////////////////////////////////////////
   BCArray3D(const size_type& nx1, const size_type& nx2, const size_type& nx3, const IndexType& val)
   {
      bcindexmatrix.resize(nx1,nx2,nx3,val);
   }
   //////////////////////////////////////////////////////////////////////////
   virtual ~BCArray3D() {}
   //////////////////////////////////////////////////////////////////////////
   inline size_type getNX1() const { return bcindexmatrix.getNX1(); }
   //////////////////////////////////////////////////////////////////////////
   inline size_type getNX2() const { return bcindexmatrix.getNX2(); }
   //////////////////////////////////////////////////////////////////////////
   inline size_type getNX3() const { return bcindexmatrix.getNX3(); }
   //////////////////////////////////////////////////////////////////////////
   void resize(const size_type& nx1, const size_type& nx2, const size_type& nx3)
   {
      bcindexmatrix.resize(nx1,nx2,nx3);
   }
   //////////////////////////////////////////////////////////////////////////
   void resize(const size_type& nx1, const size_type& nx2, const size_type& nx3, const IndexType& val)
   {
      bcindexmatrix.resize(nx1,nx2,nx3,val);
   }
   //////////////////////////////////////////////////////////////////////////
   bool validIndices(const size_type& x1, const size_type& x2, const size_type& x3)  const
   {
      if( x1 < 0 || x1 >= this->bcindexmatrix.getNX1()) return false;
      if( x2 < 0 || x2 >= this->bcindexmatrix.getNX2()) return false;
      if( x3 < 0 || x3 >= this->bcindexmatrix.getNX3()) return false;
      return true;
   }
   //////////////////////////////////////////////////////////////////////////
   inline bool hasBC(const size_type& x1, const size_type& x2, const size_type& x3)  const
   {
      return bcindexmatrix(x1,x2,x3)>=0;
   }
   //////////////////////////////////////////////////////////////////////////
   void setBC(const size_type& x1, const size_type& x2, const size_type& x3, BCClassPtr const& bc)
   {
      if( this->hasBC(x1,x2,x3) )
      {
         if( this->getBC(x1,x2,x3)==bc ) return;
         else                            this->deleteBC(x1,x2,x3);
      }

      //wenn keine frei gewordene BCs vorhanden
      if( indexContainer.empty() )
      {
         bcvector.push_back(bc);
         bcindexmatrix(x1,x2,x3) = (IndexType)bcvector.size()-1;
      }
      else
      {
         IndexType index = indexContainer.back();
         bcindexmatrix(x1,x2,x3) = index;
         bcvector[index] = bc;
         indexContainer.pop_back();
      }
   }
   //////////////////////////////////////////////////////////////////////////
   inline IndexType getBCVectorIndex(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      return bcindexmatrix(x1,x2,x3);
   }
   //////////////////////////////////////////////////////////////////////////
   inline const BCClassPtr getBC(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      IndexType index = bcindexmatrix(x1,x2,x3);
      if(index<0) return BCClassPtr(); //=> NULL Pointer

      return bcvector[index];
   }
   //////////////////////////////////////////////////////////////////////////
   inline BCClassPtr getBC(const size_type& x1, const size_type& x2, const size_type& x3)
   {
      IndexType index = bcindexmatrix(x1,x2,x3);
      if(index<0) return BCClassPtr(); //=> NULL Pointer

      return bcvector[index];
   }
   //////////////////////////////////////////////////////////////////////////
   void setSolid(const size_type& x1, const size_type& x2, const size_type& x3)
   {
      if( this->hasBC(x1,x2,x3) ) this->deleteBC(x1,x2,x3);
      bcindexmatrix(x1,x2,x3)=SOLID;
   }
   //////////////////////////////////////////////////////////////////////////
   inline bool isSolid(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      return bcindexmatrix(x1,x2,x3)==SOLID;
   }
   //////////////////////////////////////////////////////////////////////////
   void setFluid(const size_type& x1, const size_type& x2, const size_type& x3)
   {
      if( this->hasBC(x1,x2,x3) ) this->deleteBC(x1,x2,x3);
      bcindexmatrix(x1,x2,x3)=FLUID;
   }
   //////////////////////////////////////////////////////////////////////////
   //true : FLUID or BC
   //false: UNDEFINED or SOLID
   inline bool isFluid(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      int tmp = bcindexmatrix(x1,x2,x3);
      return (tmp==FLUID || tmp>=0);
   }
   //////////////////////////////////////////////////////////////////////////
   inline bool isFluidWithoutBC(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      return bcindexmatrix(x1,x2,x3)==FLUID;
   }
   //////////////////////////////////////////////////////////////////////////
   inline bool isUndefined(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      return bcindexmatrix(x1,x2,x3)==UNDEFINED;
   }
   //////////////////////////////////////////////////////////////////////////
   void setUndefined(const size_type& x1, const size_type& x2, const size_type& x3)
   {
      if( this->hasBC(x1,x2,x3) ) this->deleteBC(x1,x2,x3);
      bcindexmatrix(x1,x2,x3)=UNDEFINED;
   }
   //////////////////////////////////////////////////////////////////////////
   //inline bool isInterface(const size_type& x1, const size_type& x2, const size_type& x3) const
   //{
   //   return bcindexmatrix(x1,x2,x3)==INTERFACE;
   //}
   ////////////////////////////////////////////////////////////////////////////
   //void setInterface(const size_type& x1, const size_type& x2, const size_type& x3)
   //{
   //   if( this->hasBC(x1,x2,x3) ) this->deleteBC(x1,x2,x3);
   //   bcindexmatrix(x1,x2,x3)=INTERFACE;
   //}
   inline bool isInterfaceCF(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      return bcindexmatrix(x1, x2, x3)==INTERFACECF;
   }
   //////////////////////////////////////////////////////////////////////////
   void setInterfaceCF(const size_type& x1, const size_type& x2, const size_type& x3)
   {
      //if (this->hasBC(x1, x2, x3)) this->deleteBC(x1, x2, x3);
      bcindexmatrix(x1, x2, x3) = INTERFACECF;
   }
   //////////////////////////////////////////////////////////////////////////
   inline bool isInterfaceFC(const size_type& x1, const size_type& x2, const size_type& x3) const
   {
      return bcindexmatrix(x1, x2, x3)==INTERFACEFC;
   }
   //////////////////////////////////////////////////////////////////////////
   void setInterfaceFC(const size_type& x1, const size_type& x2, const size_type& x3)
   {
      //if (this->hasBC(x1, x2, x3)) this->deleteBC(x1, x2, x3);
      bcindexmatrix(x1, x2, x3) = INTERFACEFC;
   }
   //////////////////////////////////////////////////////////////////////////
   std::size_t getNumberOfSolidEntries() const 
   {
      const std::vector<IndexType>& data = bcindexmatrix.getDataVector();
      std::size_t counter = 0;
      for(std::size_t i=0; i<data.size(); i++)
         if(data[i]==SOLID) counter++;
      return counter;
   }
   //////////////////////////////////////////////////////////////////////////
   std::size_t getNumberOfFluidEntries() const 
   {
      const std::vector<IndexType>& data = bcindexmatrix.getDataVector();
      std::size_t counter = 0;
      for(std::size_t i=0; i<data.size(); i++)
      {
         const IndexType& tmp = data[i];
         if(tmp==FLUID || tmp>=0) counter++;
      }
      return counter;
   }
   //////////////////////////////////////////////////////////////////////////
   std::size_t getNumberOfFluidWithoutBCEntries() const
   {
      const std::vector<IndexType>& data = bcindexmatrix.getDataVector();
      std::size_t counter = 0;
      for(std::size_t i=0; i<data.size(); i++)
         if(data[i]==FLUID) counter++;
      return counter;
   }
   //////////////////////////////////////////////////////////////////////////
   std::size_t getNumberOfBCEntries() const
   {
      const std::vector<IndexType>& data = bcindexmatrix.getDataVector();
      std::size_t counter = 0;
      for(std::size_t i=0; i<data.size(); i++)
         if(data[i]>=0) counter++;
      return counter;
   }
   //////////////////////////////////////////////////////////////////////////
   std::size_t getNumberOfUndefinedEntries() const
   {
      const std::vector<IndexType>& data = bcindexmatrix.getDataVector();
      std::size_t counter = 0;
      for(std::size_t i=0; i<data.size(); i++)
         if(data[i]==UNDEFINED) counter++;
      return counter;
   }
   //////////////////////////////////////////////////////////////////////////
   std::size_t getBCVectorSize() const
   {
      return this->bcvector.size();
   }
   //////////////////////////////////////////////////////////////////////////
   std::string toString() const
   {
      std::size_t solidCounter = 0;
      std::size_t fluidCounter = 0;
      std::size_t bcCounter    = 0;
      std::size_t undefCounter = 0;

      for(int x1=0; x1<bcindexmatrix.getNX1(); x1++)
      {
         for(int x2=0; x2<bcindexmatrix.getNX2(); x2++)
         {
            for(int x3=0; x3<bcindexmatrix.getNX3(); x3++)
            {
               if(bcindexmatrix(x1,x2,x3)>=0             ) bcCounter++;
               else if(bcindexmatrix(x1,x2,x3)==FLUID    ) fluidCounter++;
               else if(bcindexmatrix(x1,x2,x3)==SOLID    ) solidCounter++;
               else if(bcindexmatrix(x1,x2,x3)==UNDEFINED) undefCounter++;
               else throw UbException(UB_EXARGS,"invalid matrixEntry");
            }
         }
      }

      std::size_t unrefEntriesInBcVector=0;
      for(std::size_t i=0; i<bcvector.size(); i++) if(!bcvector[i]) unrefEntriesInBcVector++;

      std::stringstream text;
      text<<"BCArray3D<"<<typeid(BCClass).name()<<","<<typeid(IndexType).name()<<">";
      text<<"[ entries: "<<bcindexmatrix.getNX1()<<"x"<<bcindexmatrix.getNX2();
      text<<"x"<<bcindexmatrix.getNX3()<<"=";
      text<<bcindexmatrix.getNX1()*bcindexmatrix.getNX2()*bcindexmatrix.getNX3()<<" ]:\n";
      text<<" - #fluid entries : "<<fluidCounter<<std::endl;
      text<<" - #bc    entries : "<<bcCounter<<std::endl;
      text<<" - #solid entries : "<<solidCounter<<std::endl;
      text<<" - #undef entries : "<<undefCounter<<std::endl;
      text<<" - bcvector-entries      : "<<bcvector.size()<<" (empty ones: "<<unrefEntriesInBcVector<<")\n";
      text<<" - indexContainer-entries: "<<indexContainer.size()<<std::endl;

      return text.str();
   }
//////////////////////////////////////////////////////////////////////////

   static const IndexType SOLID;     //definiton erfolgt ausserhalb!!!
   static const IndexType FLUID;     //definiton erfolgt ausserhalb!!!
   //static const IndexType INTERFACE; //definiton erfolgt ausserhalb!!!
   static const IndexType INTERFACECF; //definiton erfolgt ausserhalb!!!
   static const IndexType INTERFACEFC; //definiton erfolgt ausserhalb!!!
   static const IndexType UNDEFINED; //definiton erfolgt ausserhalb!!!

private:
   //////////////////////////////////////////////////////////////////////////
   void deleteBCAndSetType(const size_type& x1, const size_type& x2, const size_type& x3, const IndexType& type)
   {
      this->deleteBC(x1,x2,x3);

      //matrix neuen Typ zuweisen
      bcindexmatrix(x1,x2,x3) = type;
   }
   //////////////////////////////////////////////////////////////////////////
   void deleteBC(const size_type& x1,const size_type& x2, const size_type& x3)
   {
      //ueberpruefen, ob ueberhaupt BC vorhanden
      IndexType index = bcindexmatrix(x1,x2,x3);
      if(index<0) return;

      //frei gewordenen Index in den Indexcontainer schieben
      indexContainer.push_back(index);

      //element "loeschen"
      bcvector[index] = BCClassPtr();
   }

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & bcindexmatrix;
      ar & bcvector;
      ar & indexContainer;
   }
protected:
   //////////////////////////////////////////////////////////////////////////
   //CbArray3D<IndexType,IndexerX1X2X3> bcindexmatrix;  //-1 solid // -2 fluid -...
   CbArray3D<IndexType,IndexerX3X2X1> bcindexmatrix;
   std::vector<BCClassPtr> bcvector;
   std::vector<IndexType> indexContainer;
};

template< typename BCClass , class IndexType>
const IndexType BCArray3D<BCClass,IndexType>::SOLID  = -1;
template< typename BCClass , class IndexType>
const IndexType BCArray3D<BCClass,IndexType>::FLUID  = -2;
template< typename BCClass , class IndexType> 
const IndexType BCArray3D<BCClass,IndexType>::INTERFACECF = -3;
template< typename BCClass, class IndexType>
const IndexType BCArray3D<BCClass, IndexType>::INTERFACEFC = -4;
template< typename BCClass , class IndexType>
const IndexType BCArray3D<BCClass,IndexType>::UNDEFINED = -5;

#endif //D3Q19BCMATRIX_H

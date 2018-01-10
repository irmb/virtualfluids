//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef CBARRAY4D_H
#define CBARRAY4D_H

//////////////////////////////////////////////////////////////////////////
// 4D Array
// die Daten werden in einem Vector gehalten
//
// Ver 1.0
// Sept. 2006 muffmolch@gmx.de
// Ver 1.1
// Jul. 2006 - size_t + range check bei getIndex
// Ver 1.2
// Mrz. 2008 - typedefs, new index checks, NO_CB_RANGECHECK, no base class
//             assigmetcomparison between Matrices with different value_type and/or index-class
// Oct. 2008 - +isEmpty()
//
// Rangecheck aktiv, wenn:
// -debug  : not defined "NO_CB_RANGECHECK"
// -release: not defined "NO_CB_RANGECHECK" && defined "CB_RANGECHECK"
//////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbEqual.h>
#include <algorithm>
#include <typeinfo>
#include <boost/serialization/serialization.hpp>

#include <VirtualFluidsDefinitions.h>

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// IndexClasses

//IndexerX1X2X3X4:
//x4-reihen "liegen am stueck" im speicher
//optimaler schleifendurchlauf
//for(alle X1)
//  for(alle X2)
//    for(alle X3)
//      for(alle X4)
class VF_PUBLIC IndexerX1X2X3X4
{
public:
   typedef int size_type;
public:
   inline std::size_t getIndex( const size_type& x1 , const size_type& x2 , const size_type& x3 , const size_type& x4
                            , const size_type& nx1, const size_type& nx2, const size_type& nx3, const size_type& nx4 )  const
   {
      return nx4*(nx3*(nx2*x1+ x2)+x3)+x4 ;
   }
   inline std::size_t getStartIndexOfSortedArray(  const size_type& x1 , const size_type& x2 , const size_type& x3 , const size_type& x4
                                               , const size_type& nx1, const size_type& nx2, const size_type& nx3, const size_type& nx4 )  const
   {
      return  nx4*(nx3*(nx2*x1+ x2)+x3);
   }
};
//////////////////////////////////////////////////////////////////////////
// IndexClasses

//IndexerX4X3X2X1:
//x1-reihen "liegen am stueck" im speicher
//optimaler schleifendurchlauf
//for(alle X4)
//  for(alle X3)
//    for(alle X2)
//      for(alle X1)
class IndexerX4X3X2X1
{
public:
   typedef size_t size_type;
public:
   inline std::size_t getIndex( const size_type& x1 , const size_type& x2 , const size_type& x3 , const size_type& x4
      , const size_type& nx1, const size_type& nx2, const size_type& nx3, const size_type& nx4 )  const
   {
      return nx1*(nx2*(nx3*x4+ x3)+x2)+x1;
   }
   inline std::size_t getStartIndexOfSortedArray(  const size_type& x1 , const size_type& x2 , const size_type& x3 , const size_type& x4
      , const size_type& nx1, const size_type& nx2, const size_type& nx3, const size_type& nx4 )  const
   {
      return  nx1*(nx2*(nx3*x4+ x3)+x2);
   }
};
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// CbArray4D
template<typename T, typename IndexClass = IndexerX4X3X2X1>
class CbArray4D
{
public:
   typedef std::shared_ptr< CbArray4D <T,IndexClass> > CbArray4DPtr;

   typedef T                                                   value_type;
   typedef IndexClass                                          indexer_type;
   typedef typename IndexClass::size_type                      size_type;
   typedef typename std::vector< value_type >::reference       reference;
   typedef typename std::vector< value_type >::const_reference const_reference;
   typedef typename std::vector< value_type >::pointer         pointer;
   typedef typename std::vector< value_type >::const_pointer   const_pointer;

private:
   template< typename value_type2, typename IndexClass2 > friend class CbArray4D;

public:
   /*=======================================================================*/
   CbArray4D()
   {
      this->resize(0,0,0,0);
   }
   /*=======================================================================*/
   CbArray4D(const size_type& nx1,const size_type& nx2, const size_type& nx3, const size_type& nx4)
   {
      this->resize(nx1,nx2,nx3,nx4);
   }
   /*=======================================================================*/
   CbArray4D(const size_type& nx1,const size_type& nx2, const size_type& nx3, const size_type& nx4, const value_type& val)
   {
      this->resize(nx1,nx2,nx3,nx4,val);
   }
   /*=======================================================================*/
   CbArray4D(const size_type& uniformDimensionSize /*nx1=nx2=nx3=nx4*/)
   {
      this->resize(uniformDimensionSize,uniformDimensionSize,uniformDimensionSize,uniformDimensionSize);
   }
   /*=======================================================================*/
   //ubernimmt vector als daten vector! (erstellt KEINE kopie!!!, vec ist anschließend leer, da swap verwendet wird)
   CbArray4D(std::vector<value_type>& vec, const size_type& nx1,const size_type& nx2, const size_type& nx3, const size_type& nx4)
   {
      assert( (nx1*nx2*nx3*nx4)==vec.size() );
      this->data.swap(vec);
      this->resize(nx1,nx2,nx3,nx4);
   }
   /*=======================================================================*/
   CbArray4D(const CbArray4D& src)
      :  nx1(src.nx1)
       , nx2(src.nx2)
       , nx3(src.nx3)
       , nx4(src.nx4)
       , data(src.data)
   {
   }
   /*=======================================================================*/
   template< typename value_type2 >
   CbArray4D(const CbArray4D< value_type2 >& src)
      :  nx1(src.nx1)
       , nx2(src.nx2)
       , nx3(src.nx3)
       , nx4(src.nx4)
   {
      //Sourcedaten kopieren
      this->data.resize( src.data.size() );
      for(std::size_t i=0; i<data.size(); ++i)
         this->data[i] = src.data[i];
   }
   /*=======================================================================*/
   virtual ~CbArray4D()
   {
      //vector wird automatisch zerstoert
   }
   /*=======================================================================*/
   CbArray4D& operator= (const CbArray4D& rhs)
   {
      if(this == &rhs) return *this;

      this->nx1 = rhs.nx1;
      this->nx2 = rhs.nx2;
      this->nx3 = rhs.nx3;
      this->nx4 = rhs.nx4;

      //gespeicherte Datenelemente loeschen
      //Laenge anpassen
      this->data.resize(rhs.data.size());
      //gespeicherte Datenelemente loeschen
      this->data.clear();

      //Sourcedaten kopieren
      this->data = rhs.data;

      return *this;
   }
   /*=======================================================================*/
   //durch value_type2 kann man z.B. ein float Array einem double Array zuweisen!
   template< typename value_type2, typename IndexClass2 >
   CbArray4D& operator= (const CbArray4D< value_type2, IndexClass2 >& rhs)
   {
      this->nx1 = rhs.nx1;
      this->nx2 = rhs.nx2;
      this->nx3 = rhs.nx3;
      this->nx4 = rhs.nx4;

      //gespeicherte Datenelemente loeschen
      this->data.clear();
      //Laenge anpassen
      this->data.resize(rhs.data.size());

      //Sourcedaten kopieren (!! koennte anderen Indexer besitzen!!! -> operator() benutzen)
      for(int x1=0; x1<this->nx1; x1++)
         for(int x2=0; x2<this->nx2; x2++)
            for(int x3=0; x3<this->nx3; x3++)
               for(int x4=0; x4<this->nx4; x4++)
                  this->operator()(x1,x2,x3,x4) = static_cast< value_type >( rhs.operator()(x1,x2,x3,x4));

      return *this;
   }
   /*=======================================================================*/
   bool operator== (const CbArray4D& rhs) const
   {
      if( this == &rhs ) return true;

      if(   this->nx1!=rhs.nx1
         || this->nx2!=rhs.nx2
         || this->nx3!=rhs.nx3
         || this->nx4!=rhs.nx4
         || this->data.size() != rhs.data.size() )
      {
         return false;
      }

      return std::equal( this->data.begin(), this->data.end(), rhs.data.begin(), UbEqual<value_type, value_type >() );
   }
   /*=======================================================================*/
   template< typename value_type2, typename IndexClass2 >
   bool operator== (const CbArray4D< value_type2, IndexClass2 >& rhs) const
   {
      if( this->data.size() != rhs.data.size() ) return false;

      //Sourcedaten einzeln checken (!! koennte anderen Indexer besitzen!!! -> operator() benutzen)
      for(int x4=0; x4<this->nx4; x4++)
         for(int x3=0; x3<this->nx3; x3++)
            for(int x2=0; x2<this->nx2; x2++)
             for(int x1=0; x1<this->nx1; x1++)
               if( !isUbEqual(this->operator()(x1,x2,x3,x4), rhs.operator()(x1,x2,x3,x4)) )
                  return false;

      return true;
   }
   /*=======================================================================*/
   bool operator!= (const CbArray4D& rhs) const
   {
      return !(*this==rhs);
   }
   /*=======================================================================*/
   template< typename value_type2, typename IndexClass2 >
   bool operator!= (const CbArray4D< value_type2, IndexClass2 >& rhs) const
   {
      return !(*this==rhs);
   }
   /*=======================================================================*/
   reference operator() (const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4)
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif

      return this->data[indexer.getIndex(x1,x2,x3,x4,nx1,nx2,nx3,nx4)];
   }
   /*=======================================================================*/
   const_reference operator() (const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4)	const
   {
      #ifdef CbArray4D_RANGECHECKING
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif

      return this->data[indexer.getIndex(x1,x2,x3,x4,nx1,nx2,nx3,nx4)];
   }
   /*=======================================================================*/
   pointer getStartAdressOfSortedArray(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4)
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif

      return &this->data[indexer.getStartIndexOfSortedArray(x1,x2,x3,x4,nx1,nx2,nx3,nx4)];
   }
   /*=======================================================================*/
   const_pointer getStartAdressOfSortedArray(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4)  const
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif

      return &this->data[indexer.getStartIndexOfSortedArray(x1,x2,x3,x4,nx1,nx2,nx3,nx4)];
   }
   /*=======================================================================*/
   void setObject(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4, const value_type& value)
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif

      this->data[indexer.getIndex(x1,x2,x3,x4,nx1,nx2,nx3,nx4)] = value;
   }
   /*=======================================================================*/
   reference getObject(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4)
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif

      return this->data[indexer.getIndex(x1,x2,x3,x4,nx1,nx2,nx3,nx4)];
   }
   /*=======================================================================*/
   const_reference getObject(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4) const
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif
      return (*this)(x1,x2,x3,x4,nx1,nx2,nx3,nx4);
   }
   /*=======================================================================*/
   bool      isEmpty() const { return data.empty(); }
   size_type getNX1()  const { return this->nx1;    }
   size_type getNX2()  const { return this->nx2;    }
   size_type getNX3()  const { return this->nx3;    }
   size_type getNX4()  const { return this->nx4;    }
   /*=======================================================================*/
   void reset(const value_type& val)
   {
      std::fill( this->data.begin(), this->data.end(), val );
   }
   /*=======================================================================*/
   std::string toString() const
   {
      std::stringstream text;
      text<<std::setprecision(19);
      for(size_type x1=0; x1<this->nx1; x1++)
      {
      	for(size_type x2=0; x2<this->nx2; x2++)
      	{
         	for(size_type x3=0; x3<this->nx3; x3++)
         	{
               for(size_type x4=0; x4<this->nx4; x4++)
               {
                 	text<<(*this)(x1,x2,x3,x4)<<", ";
               }
               text<<std::endl;
         	}
         	text<<std::endl;
      	}
      	text<<std::endl<<std::endl;
      }

      return text.str();
   }
   /*=======================================================================*/
   std::string getInfo() const
   {
      std::stringstream text;
      text<<"CbArray4D< storageType="<<typeid(T).name()<<", indexer="<<typeid(IndexClass).name()<<" >";
      text<<"( nx1="<<this->nx1<<", nx2="<<this->nx2<<", nx3="<<this->nx3<<", nx4="<<this->nx4<<")";
      return text.str();
   }
   /*=======================================================================*/
   void resize(const size_type& uniformDimensionSize) { this->resize(uniformDimensionSize,uniformDimensionSize,uniformDimensionSize); }
   /*=======================================================================*/
   void resize(const size_type& nx1, const size_type& nx2, const size_type& nx3, const size_type& nx4)
   {
      this->nx1 = nx1;
      this->nx2 = nx2;
      this->nx3 = nx3;
      this->nx4 = nx4;
      this->data.resize(nx1*nx2*nx3*nx4);
   }
   /*=======================================================================*/
   void resize(const size_type& nx1, const size_type& nx2, const size_type& nx3, const size_type& nx4, const value_type& val)
   {
      this->nx1 = nx1;
      this->nx2 = nx2;
      this->nx3 = nx3;
      this->nx4 = nx4;
      this->data.resize(nx1*nx2*nx3*nx4,val);
   }
   /*=======================================================================*/
   std::vector< value_type >& getDataVector() { return this->data; }
   /*=======================================================================*/
   const std::vector< value_type >& getDataVector() const { return this->data; }
   /*=======================================================================*/
   inline std::size_t getDataVectorIndex(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4) const
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if( !this->indicesInRange(x1,x2,x3,x4) )
            UB_THROW( UbException(UB_EXARGS,getExceptionErrorString(x1,x2,x3,x4)) );
      #endif

      return indexer.getIndex(x1,x2,x3,x4,nx1,nx2,nx3,nx4);
   }

#ifdef CAB_RCF
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & nx1;
      ar & nx2;
      ar & nx3;
      ar & nx4;
      ar & data;
   }
#endif //CAB_RCF

protected:
   /*=======================================================================*/
   //success -> true
   //else    -> false
   inline bool indicesInRange(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4) const
   {
      if(   x1 < 0 || x1 >= this->nx1
         || x2 < 0 || x2 >= this->nx2
         || x3 < 0 || x3 >= this->nx3
         || x4 < 0 || x4 >= this->nx4 )
      {
         return false;
      }
      return true;
   }
   /*=======================================================================*/
   std::string getExceptionErrorString(const size_type& x1, const size_type& x2, const size_type& x3, const size_type& x4) const
   {
      std::stringstream out("index out of range - ");
      out<<"("<<x1<<","<<x2<<","<<x3<<","<<x4<<") not in ("<<nx1<<","<<nx2<<","<<nx3<<","<<nx4<<")";
      return out.str();
   }
   /*=======================================================================*/

protected:
   size_type    nx1;
   size_type    nx2;
   size_type    nx3;
   size_type    nx4;
   indexer_type indexer;
   std::vector< value_type > data;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & nx1;
      ar & nx2;
      ar & nx3;
      ar & nx4;
      ar & data;
   }
};

#endif //CBARRAY4D_H

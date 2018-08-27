//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef CBVECTOR_H
#define CBVECTOR_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
   #include <RCF/ByteBuffer.hpp>
#endif

#include <vector>
#include <algorithm> //for std::swap
#include <typeinfo>  //for typeid
#include <PointerDefinitions.h>    //for memcopy

#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbEqual.h>

/*=========================================================================*/
/*  CbVector                                                               */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 08.11.07
@version 1.1 - 09.02.08
@version 1.2 - 23.04.08 - swap added
@version 1.3 - 08.05.08 - boosting up serialization performance!
*/

/*
usage: ...
Da es Voraussetzun bei doeser Klasse war, dass lediglich der Typ als
template-parameter miteingeht, muss der allcocator eine abstrakte klasse sein
ansonsten hätte sich hier der allokator als zweites argument
wie beim STL vector angeboten, womit man auch keinen pointer speichern muesste.
Im letzteren Fall würde aber jeweils ein bestimmeter Klassentyp in Abhaengigkeit
des allokators zur compilezeit erzeugt. Problem wir wollen ein und denselben
typ benutzen und nur der allokator innerhalb der klasse soll sich unterscheiden
//
// Rangecheck aktiv, wenn:
// -debug  : not defined "NO_CB_RANGECHECK"
// -release: not defined "NO_CB_RANGECHECK" && defined "CB_RANGECHECK"
*/

template< typename T > class CbVectorAllocator;
template< typename T > class CbVectorAllocatorStd;

//////////////////////////////////////////////////////////////////////////
template< typename T >
class CbVector
{
public:
   typedef T           value_type;
   typedef value_type* pointer;
   typedef std::size_t size_type;

   friend class CbVectorAllocator<value_type>; //um auf ptrData und dataSize zugreifen zu koennen!

public:
   /*==========================================================*/
   CbVector( CbVectorAllocator<value_type>* const& allocator = new CbVectorAllocatorStd<value_type> )
      :  ptrData(NULL)
       , dataSize(0)
       , allocator(allocator)
   {
      this->allocator->alloc(*this,0,value_type());
   }
   /*==========================================================*/
   CbVector( const size_type size, CbVectorAllocator<value_type>* const& allocator = new CbVectorAllocatorStd<value_type>, const value_type& value=value_type() )
      :  ptrData(NULL)
       , dataSize(0)
       , allocator(allocator)
   {
      this->allocator->alloc(*this,size,value);
   }
   /*==========================================================*/
   virtual ~CbVector()
   {
      if(allocator)
      {
         this->allocator->dealloc(*this);
         delete allocator;
         allocator=NULL;
      }
   }
   /*=======================================================================*/
   CbVector& operator= (const CbVector& src)
   {
      if(this == &src) return *this;

      //gespeicherte Datenelemente loeschen
      //Laenge anpassen
      this->allocator->resize(*this, src.size());

      //gespeicherte Datenelemente kopieren
      if( !src.empty() ) 
      {
         memcpy( (char*)ptrData, (char*)&src[0], src.size()*sizeof(value_type) ); 
         //for(size_type i=0; i<src.size(); i++)
         //   (*this)[i] = src[i];
      }

      return *this;
   }
   /*=======================================================================*/
   CbVector& operator= (const std::vector< value_type >& src)
   {
      //gespeicherte Datenelemente loeschen
      //Laenge anpassen
      this->allocator->resize(*this, src.size());

      //gespeicherte Datenelemente kopieren
      if( !src.empty() ) 
      {
         memcpy( (char*)ptrData, (char*)&src[0], src.size()*sizeof(value_type) ); 
         //for(size_type i=0; i<src.size(); i++)
         //   (*this)[i] = src[i];
      }
      
      return *this;
   }
   /*=======================================================================*/
   bool operator== (const CbVector& rhs) const
   {
      if( this           == &rhs         ) return true;
      if( this->dataSize != rhs.dataSize ) return false;

      for(size_type i=0; i<rhs.size(); i++)
         if( !isUbEqual( this->operator[](i), rhs.operator[](i) ) )
            return false;

      return true;
   }
   /*==========================================================*/
   void setAllocator( CbVectorAllocator<value_type>* const& allocator )
   {
      if(this->allocator)
      {
         if(this->allocator==allocator) return;
         this->allocator->dealloc(*this);
         delete this->allocator;
      }
      this->allocator = allocator;
      this->allocator->alloc(*this,0);
   }
   /*==========================================================*/
   size_type size() const { return dataSize; }
   /*==========================================================*/
   bool empty() const { return dataSize==0; }
   /*==========================================================*/
   bool resize(const size_type& dataSize)
   {
      return allocator->resize(*this, dataSize);
   }
   /*==========================================================*/
   bool resize(const size_type& dataSize, const value_type& value)
   {
      return allocator->resize(*this, dataSize, value);
   }
   /*==========================================================*/
   void swap(CbVector& rhs)
   {
      if( this == &rhs ) return;

      std::swap( this->ptrData  , rhs.ptrData   );
      std::swap( this->dataSize , rhs.dataSize  );
      std::swap( this->allocator, rhs.allocator );
   }
   /*==========================================================*/
   value_type& operator[](const size_type& i)
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if(i>=dataSize) 
            UB_THROW( UbException(UB_EXARGS,"T="+(std::string)typeid(*this).name()+UbSystem::toString(i)+" out of range (size="+UbSystem::toString(dataSize)+")") );
      #endif // _DEBUG

      return ptrData[i];
   }
   /*==========================================================*/
   const value_type& operator[](const size_type& i) const
   {
      #if !defined(NO_CB_RANGECHECK) && ( defined(_DEBUG) || defined(CB_RANGECHECK) )
         if(i>=dataSize) 
            UB_THROW( UbException(UB_EXARGS,"T="+(std::string)typeid(*this).name()+UbSystem::toString(i)+" out of range (size="+UbSystem::toString(dataSize)+")") );
      #endif // _DEBUG

      return ptrData[i];
   }
   /*==========================================================*/
   CbVectorAllocator<value_type>* getAllocator() const { return allocator; }
   /*==========================================================*/
   #ifdef CAB_RCF
      template<typename Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
         if( ArchiveTools::isWriting(ar) )
         {
            ar & allocator;
            ar & dataSize; //!!!erst hier

            //old:
            //for(size_type i=0; i<dataSize; i++)
            // ar & ptrData[i];

            //new and boosting to the sky:
            RCF::ByteBuffer byteBuffer( (char*) &ptrData[0], dataSize*sizeof(value_type) );
            ar & byteBuffer;
         }
         else
         {
            CbVectorAllocator<value_type>* tmpCbVectorAllocator(NULL);
            size_type tmpInteger;
            ar & tmpCbVectorAllocator;
            ar & tmpInteger;
            this->setAllocator(tmpCbVectorAllocator);
            allocator->resize(*this,tmpInteger);

            //old:
            //for(size_type i=0; i<dataSize; i++)
            // ar & ptrData[i];

            //new and boosting to the sky:
            RCF::ByteBuffer byteBuffer;
            ar & byteBuffer;
            memcpy( (char*)ptrData, byteBuffer.getPtr(), byteBuffer.getLength() ); 
         }
      }
   #endif //CAB_RCF

private:
   value_type* ptrData;
   size_type   dataSize;
   CbVectorAllocator<value_type>* allocator;
   CbVector<value_type>(const CbVector<value_type>& src);
   //CbVector<value_type>& operator=(const CbVector<value_type>& src);
};

//////////////////////////////////////////////////////////////////////////
// CbVectorAllocator-Interface
//////////////////////////////////////////////////////////////////////////
template< typename T >
class CbVectorAllocator
{
public:
   typedef typename CbVector<T>::value_type          value_type;
   typedef typename CbVector<value_type>::size_type  size_type;

public:
   CbVectorAllocator() {}
   virtual ~CbVectorAllocator() {}

   virtual bool alloc(CbVector< value_type >& vec, const size_type& dataSize, const value_type& value=value_type()) = 0;
   virtual bool resize(CbVector< value_type >& vec, const size_type& dataSize, const value_type& value=value_type()) = 0;
   virtual bool dealloc(CbVector< value_type >& vec) = 0;

#ifdef CAB_RCF
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
   }
#endif //CAB_RCF

protected:
   //folgende Methoden ersparen eine friend Deklaierung aller moeglichen Allocatoren
   //denn durch diese beiden Methoden haben sie exklusive Zugriffsrechte!
   //**********************************************************************************//
   inline value_type*& ptrDataOf( CbVector< value_type >& vec )
   {
      if( vec.getAllocator()!=this ) UB_THROW( UbException(UB_EXARGS,"allocator is not member of vec!") );
      return vec.ptrData;
   }
   //**********************************************************************************//
   inline size_type& dataSizeOf( CbVector< value_type >& vec )
   {
      if( vec.getAllocator()!=this ) UB_THROW( UbException(UB_EXARGS,"allocator is not member of vec!") );
      return vec.dataSize;
   }
};

#ifdef RCF_USE_SF_SERIALIZATION
SF_NO_CTOR(CbVectorAllocator<double>);
SF_NO_CTOR(CbVectorAllocator<float>);
#endif //RCF_USE_SF_SERIALIZATION


//////////////////////////////////////////////////////////////////////////
// CbVectorAllocatorStd
//////////////////////////////////////////////////////////////////////////
template< typename T >
class CbVectorAllocatorStd : public CbVectorAllocator<T>
{
public:
   //typedefs wiederholen, da Basisklasse = template -> "Dependent-Base"-Problem
   typedef typename CbVector<T>::value_type          value_type;
   typedef typename CbVector<value_type>::size_type  size_type;

public:
   CbVectorAllocatorStd() : CbVectorAllocator<value_type>()
   {

   }
   /*==========================================================*/
   bool alloc(CbVector< value_type >& src, const size_type& dataSize, const value_type& value=value_type())
   {
      return this->resize(src,dataSize,value);
   }
   /*==========================================================*/
   bool resize(CbVector< value_type >& vec, const size_type& dataSize, const value_type& value=value_type())
   {
      if( CbVectorAllocatorStd< value_type >::dataSizeOf(vec) == dataSize) return false;

      //new array
      value_type* new_data = new value_type[dataSize];
      //copy existing data to array
      if( this->ptrDataOf(vec) )
      {
         for(size_type i=0; (i<vec.size() && i<dataSize); ++i) new_data[i] = CbVectorAllocatorStd< value_type >::ptrDataOf(vec)[i];
         delete[] this->ptrDataOf(vec);
      }
      this->ptrDataOf(vec) = new_data;
      //new value for new items
      for(size_type i=this->dataSizeOf(vec); i<dataSize; ++i) this->ptrDataOf(vec)[i] = value;
      //assign new dataSize
      this->dataSizeOf(vec) = dataSize;

      return true;
   }
   /*==========================================================*/
   bool dealloc(CbVector< value_type >& vec)
   {
      if( this->ptrDataOf(vec) )
      {
         delete[] this->ptrDataOf(vec);
         this->ptrDataOf(vec) = NULL;
      }
      this->dataSizeOf(vec) = 0;
      return true;
   }
   /*==========================================================*/
   #ifdef CAB_RCF
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
         serializeParent< CbVectorAllocator<value_type> >(ar, *this);
      }
   #endif //CAB_RCF

private:
};


#ifdef RCF_USE_SF_SERIALIZATION
   UB_AUTO_RUN_NAMED(   SF::registerType< CbVectorAllocatorStd<double> >(" CbVectorAllocatorStd<double> ")       , SF_CbVectorAllocatorStd_double );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< CbVectorAllocator<double>, CbVectorAllocatorStd<double> >() ), SF_CbVectorAllocatorStd_double_BD1 );

   UB_AUTO_RUN_NAMED(   SF::registerType< CbVectorAllocatorStd<float> >(" CbVectorAllocatorStd<float> "  )       , SF_CbVectorAllocatorStd_float  );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< CbVectorAllocator<float> , CbVectorAllocatorStd<float> >()  ), SF_CbVectorAllocatorStd_float_BD2 );
#endif //RCF_USE_SF_SERIALIZATION

#endif //CBVECTOR_H

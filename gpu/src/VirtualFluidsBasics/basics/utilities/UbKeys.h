//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBKEYS_H
#define UBKEYS_H

#include <iostream>


#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

/*=========================================================================*/
/*  UbKeys                                                             */
/*                                                                         */
/**
namespace for global Keys (e.g. for STL-maps)
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 08.08.07
*/ 

/*
usage: ...
*/

namespace UbKeys
{
   //nested class
   template< typename T1, typename T2 = T1 >
   class Key2 
   {
   public:
      //////////////////////////////////////////////////////////////////////////
      //Konstruktoren
      Key2(const T1& t1, const T2& t2) 
         : t1(t1), t2(t2)
      { 
      }
      /*==========================================================*/
      Key2& operator=(const Key2& srcKey) 
      {
         if(this == &srcKey ) return *this;

         t1    = srcKey.t1;
         t2    = srcKey.t2;

         return *this;
      }
      /*==========================================================*/
      T1 getT1() const { return t1; }
      T2 getT2() const { return t2; }

      //////////////////////////////////////////////////////////////////////////
      //global ueberladene Operatoren
      friend inline bool operator<(const Key2& lhsKey,const Key2& rhsKey)  
      {
         if( lhsKey.t1 < rhsKey.t1 ) return true;
         if( lhsKey.t1 > rhsKey.t1 ) return false;
         if( lhsKey.t2 < rhsKey.t2 ) return true;

         return false;
      }
      /*==========================================================*/
      friend inline bool operator==(const Key2& lhsKey, const Key2& rhsKey)  
      {
         if(lhsKey.t1 != rhsKey.t1 ) return false;
         if(lhsKey.t2 != rhsKey.t2 ) return false;

         return true;
      }
      //ueberladene Operatoren
      friend inline bool operator!=(const Key2& lhsKey, const Key2& rhsKey)  
      {
         return !(lhsKey == rhsKey); 
      }
      //ueberladene Operatoren
      /*==========================================================*/
      friend inline std::ostream& operator << (std::ostream& os, const Key2& key) 
      {
         os<<"Key2<"<<typeid(T1).name()<<","<<typeid(T2).name()<<">,("<<key.t1<<","<<key.t2<<")";
         return os;
      }
      /*==========================================================*/
      #ifdef CAB_RCF
	      template<class Archive>
	      void serialize(Archive & ar, const unsigned int version)
	      {
		      ar & t1;
		      ar & t2;
	      }
      #endif //CAB_RCF

   private:
      //////////////////////////////////////////////////////////////////////////
      //private Member
      T1 t1;
      T2 t2;

   };

   //////////////////////////////////////////////////////////////////////////
   //
   //////////////////////////////////////////////////////////////////////////
   template< typename T1, typename T2 = T1, typename T3 = T1 >
   class Key3 
   {
   public:
      //////////////////////////////////////////////////////////////////////////
      //Konstruktoren
      Key3() : t1(0), t2(0), t3(0)
      {

      }
      Key3(const T1& t1, const T2& t2, const T3& t3) 
         : t1(t1), t2(t2), t3(t3)
      { 
      }
      /*==========================================================*/
      T1 getT1() const { return t1; }
      T2 getT2() const { return t2; }
      T3 getT3() const { return t3; }
      /*==========================================================*/
      Key3& operator=(const Key3& srcKey) 
      {
         if(this == &srcKey ) return *this;

         t1    = srcKey.t1;
         t2    = srcKey.t2;
         t3    = srcKey.t3;

         return *this;
      }

      //////////////////////////////////////////////////////////////////////////
      //global ueberladene Operatoren
      friend inline bool operator<(const Key3& lhsKey,const Key3& rhsKey)  
      {
         if( lhsKey.t1 < rhsKey.t1 ) return true;
         if( lhsKey.t1 > rhsKey.t1 ) return false;
         if( lhsKey.t2 < rhsKey.t2 ) return true;
         if( lhsKey.t2 > rhsKey.t2 ) return false;
         if( lhsKey.t3 < rhsKey.t3 ) return true;

         return false;
      }
      /*==========================================================*/
      friend inline bool operator==(const Key3& lhsKey,const Key3& rhsKey)  
      {
         if(lhsKey.t1 != rhsKey.t1 ) return false;
         if(lhsKey.t2 != rhsKey.t2 ) return false;
         if(lhsKey.t3 != rhsKey.t3 ) return false;

         return true;
      }
      /*==========================================================*/
      //ueberladene Operatoren
      friend inline bool operator!=(const Key3& lhsKey, const Key3& rhsKey) 
      {
         return !(lhsKey == rhsKey); 
      }

      //ueberladene Operatoren
      /*==========================================================*/
      friend inline std::ostream& operator << (std::ostream& os, const Key3& key) 
      {
         os<<"Key3<"<<typeid(T1).name()<<","<<typeid(T2).name()<<","<<typeid(T3).name();
         os<<">,("<<key.t1<<","<<key.t2<<","<<key.t3<<")";
         return os;
      }
      /*==========================================================*/
      #ifdef CAB_RCF
	      template<class Archive>
	      void serialize(Archive & ar, const unsigned int version)
	      {
		      ar & t1;
		      ar & t2;
            ar & t3;
         }
      #endif //CAB_RCF

   private:
      //////////////////////////////////////////////////////////////////////////
      //private Member
      T1 t1;
      T2 t2;
      T3 t3;

   };

   //////////////////////////////////////////////////////////////////////////
   //
   //////////////////////////////////////////////////////////////////////////
   template< typename T1, typename T2 = T1, typename T3 = T1, typename T4 = T1 >
   class Key4 
   {
   public:
      //////////////////////////////////////////////////////////////////////////
      //Konstruktoren
      Key4(const T1& t1, const T2& t2, const T3& t3, const T4& t4) 
         : t1(t1), t2(t2), t3(t3), t4(t4)
      { 
      }
      /*==========================================================*/
      T1 getT1() const { return t1; }
      T2 getT2() const { return t2; }
      T3 getT3() const { return t3; }
      T4 getT4() const { return t4; }
      /*==========================================================*/
      Key4& operator=(const Key4& srcKey) 
      {
         if(this == &srcKey ) return *this;

         t1    = srcKey.t1;
         t2    = srcKey.t2;
         t3    = srcKey.t3;
         t4    = srcKey.t4;

         return *this;
      }
      //////////////////////////////////////////////////////////////////////////
      //global ueberladene Operatoren
      friend inline bool operator<(const Key4& lhsKey,const Key4& rhsKey)  
      {
         if( lhsKey.t1 < rhsKey.t1 ) return true;
         if( lhsKey.t1 > rhsKey.t1 ) return false;
         if( lhsKey.t2 < rhsKey.t2 ) return true;
         if( lhsKey.t2 > rhsKey.t2 ) return false;
         if( lhsKey.t3 < rhsKey.t3 ) return true;
         if( lhsKey.t3 > rhsKey.t3 ) return false;
         if( lhsKey.t4 < rhsKey.t4 ) return true;

         return false;
      }
      /*==========================================================*/
      friend inline bool operator==(const Key4& lhsKey,const Key4& rhsKey)  
      {
         if(lhsKey.t1 != rhsKey.t1 ) return false;
         if(lhsKey.t2 != rhsKey.t2 ) return false;
         if(lhsKey.t3 != rhsKey.t3 ) return false;
         if(lhsKey.t4 != rhsKey.t4 ) return false;

         return true;
      }

      //ueberladene Operatoren
      friend inline bool operator!=(const Key4& lhsKey, const Key4& rhsKey) 
      {
         return !(lhsKey == rhsKey); 
      }
      //ueberladene Operatoren
      /*==========================================================*/
      friend inline std::ostream& operator << (std::ostream& os, const Key4& key) 
      {
         os<<"Key4<"<<typeid(T1).name()<<","<<typeid(T2).name()<<","<<typeid(T3).name()<<","<<typeid(T4).name();
         os<<">,("<<key.t1<<","<<key.t2<<","<<key.t3<<","<<key.t4<<")";
         return os;
      }
      /*==========================================================*/
      #ifdef CAB_RCF
	      template<class Archive>
	      void serialize(Archive & ar, const unsigned int version)
	      {
		      ar & t1;
		      ar & t2;
            ar & t3;
            ar & t4;
         }
      #endif //CAB_RCF

   private:
      //////////////////////////////////////////////////////////////////////////
      //private Member
      T1 t1;
      T2 t2;
      T3 t3;
      T4 t4;

   };
}

#endif //UBKEYS_H

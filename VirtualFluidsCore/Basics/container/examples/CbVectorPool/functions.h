#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>

#include <basics/utilities/UbTuple.h>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbTiming.h>

#include <basics/memory/MbSmartPtr.h>

#include <basics/container/CbVector.h>
#include <basics/container/CbVectorPool.h>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

typedef long double value_type;
typedef MbSmartPtr<CbVector< value_type > > CbVectorPtr;
typedef MbSmartPtr<vector< value_type > >   StlVectorPtr;

/*==========================================================*/
template<typename T>
inline void setValues(vector<T>& stlvec, CbVector<T>& cbvec, CbVector<T>& cbpoolvec)
{
   if(stlvec.size() != cbvec.size() || stlvec.size() != cbpoolvec.size() )
   {
      cerr<<"sizes:"<<endl;
      cerr<<"stlvec... = "<<(int)stlvec.size()<<endl;
      cerr<<"cbvec.... = "<<(int)cbvec.size()<<endl;
      cerr<<"cbpoolvec = "<<(int)cbpoolvec.size()<<endl;
      throw UB_THROW( UbException("setValues - sizeCheck failed") );
   }
   static value_type stlVal    = 1;
   static value_type cbVal     = 1;
   static value_type cbPoolVal = 1;

   for(size_t i=0; i<cbvec.size(); i++) stlvec[i]    = stlVal   ++;
   for(size_t i=0; i<cbvec.size(); i++) cbvec[i]     = cbVal    ++;
   for(size_t i=0; i<cbvec.size(); i++) cbpoolvec[i] = cbPoolVal++;
}
/*==========================================================*/
template<typename T>
inline void setValues(vector< StlVectorPtr >& stlvecs, vector< CbVectorPtr >& cbvecs, vector< CbVectorPtr >& cbpoolvecs)
{
   if(stlvecs.size() != cbvecs.size() || stlvecs.size() != cbpoolvecs.size() )
   {
      cerr<<"sizes:"<<endl;
      cerr<<"stlvec... = "<<(int)stlvecs.size()<<endl;
      cerr<<"cbvec.... = "<<(int)cbvecs.size()<<endl;
      cerr<<"cbpoolvec = "<<(int)cbpoolvecs.size()<<endl;
      throw UB_THROW( UbException("setValues glob - sizeCheck failed") );
   }

   for(size_t i=0; i<cbvecs.size(); i++)
      setValues(*stlvecs[i],*cbvecs[i],*cbpoolvecs[i]);
}
/*==========================================================*/
template<typename T>
inline void resize(vector<T>& stlvec, CbVector<T>& cbvec, CbVector<T>& cbpoolvec, std::size_t size, const T& val)
{
   stlvec.resize(size,val);
   cbvec.resize(size,val);
   cbpoolvec.resize(size,val);
}
/*==========================================================*/
template<typename T>
inline void resize(vector< StlVectorPtr >& stlvecs, vector< CbVectorPtr >& cbvecs, vector< CbVectorPtr >& cbpoolvecs, std::size_t size, const value_type& val, bool timed=false)
{
   if(stlvecs.size() != cbvecs.size() || stlvecs.size() != cbpoolvecs.size() )
   {
      cerr<<"sizes:"<<endl;
      cerr<<"stlvec... = "<<(int)stlvecs.size()<<endl;
      cerr<<"cbvec.... = "<<(int)cbvecs.size()<<endl;
      cerr<<"cbpoolvec = "<<(int)cbpoolvecs.size()<<endl;
      throw UB_THROW( UbException("resize glob - sizeCheck failed") );
   }

   if(timed)
   {
      UbTimer timer;
      timer.start(); for(size_t i=0; i<cbvecs.size(); i++) stlvecs[i]->resize(size,val);    if(timed) cout<<"stl-resize    in "<<timer.stop()<<"s"<<endl;
      timer.start(); for(size_t i=0; i<cbvecs.size(); i++) cbvecs[i]->resize(size,val);     if(timed) cout<<"cbStd-resize  in "<<timer.stop()<<"s"<<endl;
      timer.start(); for(size_t i=0; i<cbvecs.size(); i++) cbpoolvecs[i]->resize(size,val); if(timed) cout<<"cbPool-resize in "<<timer.stop()<<"s"<<endl;
   }
   else
   {
      for(size_t i=0; i<cbvecs.size(); i++)
         resize(*stlvecs[i],*cbvecs[i],*cbpoolvecs[i],size,val);
   }
}
/*==========================================================*/
inline void createVecs(size_t number, int size,vector< StlVectorPtr >& stlvecs, vector< CbVectorPtr >& cbvecs, vector< CbVectorPtr >& cbpoolvecs, CbVectorPool<value_type>*& pool, bool timed=false)
{
   UbTimer timer;
   timer.start(); for(size_t i=0; i<number; i++) stlvecs.push_back(StlVectorPtr(new vector<value_type>(size)));                                                  if(timed) cout<<"stl-createVecs    in "<<timer.stop()<<"s"<<endl;
   timer.start(); for(size_t i=0; i<number; i++) cbvecs.push_back(CbVectorPtr(new CbVector<value_type>(size)));                                                  if(timed) cout<<"cbStd-createVecs  in "<<timer.stop()<<"s"<<endl;
   timer.start(); for(size_t i=0; i<number; i++) cbpoolvecs.push_back(CbVectorPtr(new CbVector<value_type>(size,new CbVectorAllocatorPool<value_type>(pool))));  if(timed) cout<<"cbPool-createVecs in "<<timer.stop()<<"s"<<endl;

   for(size_t i=0; i<cbvecs.size(); i++) setValues(*stlvecs.back(),*cbvecs.back(),*cbpoolvecs.back());
}
/*==========================================================*/
inline void createVecs(size_t number, size_t size, const value_type& val,vector< StlVectorPtr >& stlvecs, vector< CbVectorPtr >& cbvecs, vector< CbVectorPtr >& cbpoolvecs, CbVectorPool<value_type>*& pool, bool timed=false)
{
   UbTimer timer;
   timer.start(); for(size_t i=0; i<number; i++) stlvecs.push_back(StlVectorPtr(new vector<value_type>(size,val)));                                                  if(timed) cout<<"stl-createVecs    in "<<timer.stop()<<"s"<<endl;
   timer.start(); for(size_t i=0; i<number; i++) cbvecs.push_back(CbVectorPtr(new CbVector<value_type>(size,new CbVectorAllocatorStd<value_type>(),val)));           if(timed) cout<<"cbStd-createVecs  in "<<timer.stop()<<"s"<<endl;
   timer.start(); for(size_t i=0; i<number; i++) cbpoolvecs.push_back(CbVectorPtr(new CbVector<value_type>(size,new CbVectorAllocatorPool<value_type>(pool),val)));  if(timed) cout<<"cbPool-createVecs in "<<timer.stop()<<"s"<<endl;
}
/*==========================================================*/
template<typename T>
inline void equalCheck(vector<T>& stlvec, CbVector<T>& cbvec, CbVector<T>& cbpoolvec)
{
   if(stlvec.size() != cbvec.size() || stlvec.size() != cbpoolvec.size() )
   {
      cerr<<"sizes:"<<endl;
      cerr<<"stlvec... = "<<(int)stlvec.size()<<endl;
      cerr<<"cbvec.... = "<<(int)cbvec.size()<<endl;
      cerr<<"cbpoolvec = "<<(int)cbpoolvec.size()<<endl;
      throw UB_THROW( UbException("equalCheck - sizeCheck failed") );
   }

   bool check=true;
   for(size_t i=0; i<cbvec.size(); i++)
      if(stlvec[i] != cbvec[i] || stlvec[i] != cbpoolvec[i]  )
         check=false;

   if(!check)
   {
      cerr<<"\nstl - "; for(size_t i=0; i<cbvec.size(); i++) cout<<stlvec[i]<<" ";    cout<<endl;
      cerr<<  "cbv - "; for(size_t i=0; i<cbvec.size(); i++) cout<<cbvec[i]<<" ";     cout<<endl;
      cerr<<  "cbp - "; for(size_t i=0; i<cbvec.size(); i++) cout<<cbpoolvec[i]<<" "; cout<<endl;
      throw UB_THROW( UbException("equalCheck - equalCheck failed") );
   }
}
/*==========================================================*/
template<typename T>
void equalCheck(vector< StlVectorPtr >& stlvecs, vector< CbVectorPtr >& cbvecs, vector< CbVectorPtr >& cbpoolvecs)
{
   if(stlvecs.size() != cbvecs.size() || stlvecs.size() != cbpoolvecs.size() )
   {
      cerr<<"sizes:"<<endl;
      cerr<<"stlvec... = "<<(int)stlvecs.size()<<endl;
      cerr<<"cbvec.... = "<<(int)cbvecs.size()<<endl;
      cerr<<"cbpoolvec = "<<(int)cbpoolvecs.size()<<endl;
      throw UB_THROW( UbException("equalCheck - sizeCheck failed") );
   }

   for(size_t i=0; i<cbvecs.size(); i++)
   {
      //cout<<"equalCheck i="<<i<<"/"<<cbvecs.size()-1;
      equalCheck(*stlvecs[i],*cbvecs[i],*cbpoolvecs[i]);
      //cout<<" passed"<<endl;
   }
}
/*==========================================================*/
void accessCheck(int times,vector< StlVectorPtr >& stlvecs, vector< CbVectorPtr >& cbvecs, vector< CbVectorPtr >& cbpoolvecs)
{
   UbTimer timer;
   timer.start();
   for(size_t i=0; i<stlvecs.size(); i++)
   {
      vector<value_type>& vec = *stlvecs[i];
      for(int m=0; m<times; m++)
         for(vector<value_type>::size_type k=0; k<vec.size(); k++) vec[k] = k;
   }
   cout<<"stl-accessCheck       in "<<timer.stop()<<"s"<<endl;
   timer.start();
   for(size_t i=0; i<cbvecs.size(); i++)
   {
      CbVector<value_type>& vec = *cbvecs[i];
      for(int m=0; m<times; m++)
         for(vector<value_type>::size_type k=0; k<vec.size(); k++) vec[k] = k;
   }
   cout<<"cbStd-accessCheck     in "<<timer.stop()<<"s"<<endl;
   timer.start();
   for(size_t i=0; i<cbpoolvecs.size(); i++)
   {
      CbVector<value_type>& vec = *cbpoolvecs[i];
      for(int m=0; m<times; m++)
         for(vector<value_type>::size_type k=0; k<vec.size(); k++) vec[k] = k;
   }
   cout<<"cbPool-accessCheck    in "<<timer.stop()<<"s"<<endl;
}

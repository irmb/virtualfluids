#ifndef RBAGGREGATION_H
#define RBAGGREGATION_H

#include <sstream>
#include <iostream>
#include <map>

using namespace std;

template <class T1, class T2>
class RbAggregation 
{
private:
   string name;
   std::multimap<T1,T2> obj1Map;
   std::multimap<T2,T1> obj2Map;
                       
public:
   RbAggregation(string name)
   {
      this->name = name;
   }
   /*=========================================================================*/
   void insertPair(T1& to1, T2& to2)
   {
      obj1Map.insert(pair<T1,T2>(to1,to2));
      obj2Map.insert(pair<T2,T1>(to2,to1));
   }     
   /*=========================================================================*/
   int countObj2forObj1(T1& to1)
   {                                                                
      return((int)obj1Map.count(to1));
   }

   /*=========================================================================*/
   int countObj1forObj2(T2& to2)
   {
      return((int)obj2Map.count(to2));
   }
   /*=========================================================================*/
   vector<T2> getObj2vectorForObj1(T1& to1)
   {
      vector<T2> obj2vector;
      unsigned number = (unsigned)obj1Map.count(to1);
      typedef std::multimap<T1, T2>::iterator obj1MapIterator = obj1Map.find(to1);
      for(unsigned u =0; u<number; u++) 
      {
         obj2vector.push_back(obj1MapIterator->second);
         obj1MapIterator++;
      }
      return obj2vector;
   }
   ///*=========================================================================*/
   vector<T1>  getObj1vectorForObj2(T2& to2)
   {
      vector<T1> obj1vector;
      unsigned number = (unsigned)obj2Map.count(to2);
      typedef std::multimap<T2, T1>::iterator obj2MapIterator = obj2Map.find(to2);
      for(unsigned u =0; u<number; u++) 
      {
         obj1vector.push_back(obj2MapIterator->second);
         obj2MapIterator++;
      }
      return obj1vector;
   }
};
/*=========================================================================*/
#endif



/**
* @file PriorityQueueDecompositor.h
* @brief Priority Queue for threads decomposition.
* @author Kostyantyn Kucher
* @date 06/06/2011
*/
#ifndef PRIORITYQUEUEDECOMPOSITOR_H
#define PRIORITYQUEUEDECOMPOSITOR_H

#include <algorithm>
#include <vector>
#include <map>

struct sortMinMax {
   bool operator() (int i,int j) const { return (i<j);}
};

struct sortMaxMin {
   bool operator() (int i,int j) const { return (i>j);}
};

template <class T>
class PriorityQueueDecompositor
{
public:
   PriorityQueueDecompositor(const std::vector<T>& objcts, const std::vector<int>& weights, const int& numberOfParts)
   {
      for (int i = 0; i < (int)objcts.size(); i++)
      {
         objects.insert(std::pair<int, T>(weights[i], objcts[i]));
      }
      for (int i = 0; i < numberOfParts; i++)
      {
         std::vector<T> part;
         parts.insert(std::pair<int,std::vector<T> >(0, part));
      }
   }
   virtual ~PriorityQueueDecompositor()
   {

   }
   void getDecomposition(std::vector< std::vector<T> >& prts)
   {
      for( itOb=objects.begin() ; itOb != objects.end(); itOb++)
      {
         itP = parts.begin();
         int weight = (*itP).first;
         std::vector<T> obj = (*itP).second;
         parts.erase(itP);
         weight += (*itOb).first;
         obj.push_back((*itOb).second);
         parts.insert(std::pair<int,std::vector<T> >(weight, obj));
      }

      for( itP=parts.begin() ; itP != parts.end(); itP++)
      {
         prts.push_back((*itP).second);
      }
   }
protected:
private:
   std::multimap<int, T, sortMaxMin> objects;
   typename std::multimap<int, T, sortMaxMin>::iterator itOb;
   std::multimap<int, std::vector<T>, sortMinMax> parts;
   typename std::multimap<int, std::vector<T>, sortMinMax>::iterator itP;
};

#endif

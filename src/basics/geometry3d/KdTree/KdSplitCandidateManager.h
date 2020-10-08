#ifndef KDSPLITCANDIDATEMANAGER_H
#define KDSPLITCANDIDATEMANAGER_H

#include <geometry3d/KdTree/KdSplitCandidate.h>

#include <map>
#include <vector>
#include <algorithm>


namespace Kd
{
   template< typename T >
   class SplitCandidateManager  
   {
   public:
      SplitCandidateManager()
          
      {

      }
      /* ======================================================================================= */
      SplitCandidate<T>& operator[] (const std::size_t& i)
      { 
         #ifdef DEBUG
            return splitCandidatesVec.at(i);
         #else
            return splitCandidatesVec[i];
         #endif  
      }
      /* ======================================================================================= */
      typename std::vector< SplitCandidate< T > >::size_type size()
      { 
         return splitCandidatesVec.size();
      }
      /* ======================================================================================= */
      void add(const T& pos, const int& axis, const int& starting, const int& ending, const int& np)
      {
         typename std::map<T, SplitCandidate<T> >::iterator it = splitCandidates.find(pos); 
         if ( it != splitCandidates.end() )   //split candidate is already available -> increase parameter (starting, ending and np)
         {
            SplitCandidate<T>& sc = it->second;
            sc.np       += np;
            sc.starting += starting;
            sc.ending   += ending;
         } 
         else // split candidate is not available -> add new split candidate
         {
            this->splitCandidates[pos] = SplitCandidate<T>(axis, pos, starting, ending, np);
         }
      }
      /* ======================================================================================= */
      void createSortedArray()
      {
         splitCandidatesVec.clear();
         typename std::map<T, SplitCandidate<T> >::iterator it;
         for( it=splitCandidates.begin(); it!=splitCandidates.end(); ++it)
            splitCandidatesVec.push_back(it->second);
         splitCandidates.clear();
         std::sort(splitCandidatesVec.begin(), splitCandidatesVec.end(), std::less< SplitCandidate<T> >() );
      }
      /* ======================================================================================= */

   public:
      int objects_starting_outside_left{0};
      int objects_fully_outside_node{0};

   private:
      std::map<T, SplitCandidate<T> > splitCandidates;
      std::vector< SplitCandidate<T> >    splitCandidatesVec;
   };
}

#endif //KDSPLITCANDIDATEMANAGER_H

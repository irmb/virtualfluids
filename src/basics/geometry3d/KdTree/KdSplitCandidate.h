#ifndef KDSPLITCANDIDATE_H
#define KDSPLITCANDIDATE_H

#include <basics/utilities/UbMath.h>

namespace Kd
{
   template< typename T >
   class SplitCandidate   
   {
   public:
      SplitCandidate() 
         :  
           position(0.0)
          , 
           Cn(0.0)
           
      {

      }
      /* ======================================================================================= */
      SplitCandidate(const int& axis, const T& position, const int& starting, const int& ending, const int& insidePlane)
         : np_left(false)
         , np_right(false)
         , axis(axis)
         , position(position)
         , starting(starting)
         , ending(ending)
         , np(insidePlane)
         , Cn(0.0) 
         , nr(0)  
         , nl(0)  
         , isValid(true)
      {
      }
      /* ======================================================================================= */
      bool operator!() const
      {
         return isValid; 
      }
      /* ======================================================================================= */
      friend inline bool operator< (const SplitCandidate& lhs, const SplitCandidate& rhs)  
      {
         return  lhs.position < rhs.position;
      }
      /* ======================================================================================= */

   public:
      int     axis{0};
      T       Cn;
      T       position;
      int     nl{0};
      int     nr{0};
      int     np;
      int     starting{0};
      int     ending{0};
      bool    np_left{false};
      bool    np_right{false};
      bool    isValid{false};
   };
}

#endif //KDSPLITCANDIDATE_H

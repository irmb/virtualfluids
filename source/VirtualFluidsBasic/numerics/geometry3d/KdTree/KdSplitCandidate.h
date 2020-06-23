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
         :  axis(0)
          , position(0.0)
          , starting(0)
          , ending(0)
          , np_left(false)
          , np_right(false)
          , Cn(0.0)
          , nr(0)
          , nl(0)  
          , isValid(false)
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
      int     axis;
      T       Cn;
      T       position;
      int     nl;
      int     nr;
      int     np;
      int     starting;
      int     ending;
      bool    np_left;
      bool    np_right;
      bool    isValid;
   };
}

#endif //KDSPLITCANDIDATE_H

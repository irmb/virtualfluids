#ifndef KDRAY_H
#define KDRAY_H

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbMath.h>


namespace Kd
{
   /*
   * Ray class, for use with the optimized ray-box intersection test
   * described in:
   *
   *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
   *      "An Efficient and Robust Ray-Box Intersection Algorithm"
   *      Journal of graphics tools, 10(1):49-54, 2005
   * 
   */
   template< typename T>
   class Ray 
   {
   public:
      Ray(  const T& originX   , const T& originY   , const T& originZ
          , const T& directionX, const T& directionY, const T& directionZ ) 
      {
         this->originX        = originX;
         this->originY        = originY;
         this->originZ        = originZ;

         //normierung (fuer ray-triangle-intersection)
         T oneOverLength = T(1.0/std::sqrt(  directionX*directionX 
                                            + directionY*directionY 
                                            + directionZ*directionZ ) );

         this->directionX     = directionX*oneOverLength;
         this->directionY     = directionY*oneOverLength;
         this->directionZ     = directionZ*oneOverLength;

         this->inv_directionX = T(1.0/this->directionX);   //ACHTUNG: BEWUSST KEINE ==0 Abfrage
         this->inv_directionY = T(1.0/this->directionY);   //Alg verwendet exlitzit INF
         this->inv_directionZ = T(1.0/this->directionZ);
         
         if(this->inv_directionX < 0.0) this->signX = 1;
         else                           this->signX = 0;
         if(this->inv_directionY < 0.0) this->signY = 1;
         else                           this->signY = 0;
         if(this->inv_directionZ < 0.0) this->signZ = 1;
         else                           this->signZ = 0;
      }

      T   originX;
      T   originY;
      T   originZ;
        
      T   directionX;
      T   directionY;
      T   directionZ;
        
      T   inv_directionX;
      T   inv_directionY;
      T   inv_directionZ;
      
      int signX;
      int signY;
      int signZ;
   };
} //namespace Kd

#endif //KDRAY_H

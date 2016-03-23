//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBRANDOM_H 
#define UBRANDOM_H 

#include <cstdlib> 
#include <ctime> 
#include <cassert> 
#include <cmath> 

/*=========================================================================*/
/*  UbRandom                                                             */
/*                                                                         */
/**
generates a random number
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 04.10.2007
*/ 
/*
usage: 
   int main() 
   { 
      char* hand[] = {"Schere", "Stein", "Papier"}; 
      for (unsigned u = 0; u < 20; u++) 
      { 
         cout << hand[UbRandom::rand(0, 2, 1)] << endl; 
      } 

      return 0; 
   } 
*/

class UbRandom 
{ 
private: 
   UbRandom() { std::srand(static_cast<int>(std::time(NULL)));  } 

public: 
   //returns arbitrary int value element of [min ; max]
   static inline int rand(const int& min, const int& max) 
   { 
      static UbRandom dummy; 
      assert(max - min < RAND_MAX); 
      return ( min + std::rand() % (max - min + 1) ); 
   } 
   //returns arbitrary float value element of "( (max - min) / gran ) * [min ; max]"
   //with other words: val = min+n*(max-min)/gran, n=0..gran-1
   static inline double rand(const double& min, const double& max, const double& gran)
   {
      static UbRandom dummy; 
      return (min + std::floor( std::rand() / (1.0 + RAND_MAX) * gran)* (max - min) / gran);
   }
}; 

#endif //UBRANDOM_H

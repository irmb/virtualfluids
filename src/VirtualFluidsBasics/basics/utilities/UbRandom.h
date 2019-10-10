//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file UbRandom.h
//! \ingroup utilities
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBRANDOM_H 
#define UBRANDOM_H 

#include <cstdlib> 
#include <ctime> 
#include <cassert> 
#include <cmath> 

//////////////////////////////////////////////////////////////////////////
//!
//! \brief A class generates a random number
//! \details
//! usage: 
//! \code
//!    int main() 
//!    { 
//!       char* hand[] = {"Schere", "Stein", "Papier"}; 
//!       for (unsigned u = 0; u < 20; u++) 
//!       { 
//!          cout << hand[UbRandom::rand(0, 2, 1)] << endl; 
//!       } 
//! 
//!       return 0; 
//!    } 
//! \endcode
//////////////////////////////////////////////////////////////////////////

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

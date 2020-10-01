//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef MCTYPES_H
#define MCTYPES_H

//extension by CAB
#include <vector>
#include <string>
#include <fstream>

#include <basics/utilities/UbException.h>
#include <basics/container/CbArray3D.h>
#include <basics/container/CbArray4D.h>

namespace McCubes
{ 

   #if !defined(WIN32) || defined(__CYGWIN__)
   #pragma interface
   #endif // WIN32

   //_____________________________________________________________________________
   // types
   /** unsigned char alias */
   typedef unsigned char   uchar;
   /** signed char alias */
   typedef   signed char   schar;
   /** isovalue alias */
   typedef          double real;

} //namespace McCubes

#endif //MCTYPES_H

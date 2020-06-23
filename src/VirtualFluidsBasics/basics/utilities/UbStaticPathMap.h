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
//! \file UbStaticPathMap.h
//! \ingroup utilities
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBSTATICPATHMAP_H
#define UBSTATICPATHMAP_H

#include <iostream>
#include <string>
#include <map>

#include <basics/utilities/UbSystem.h>

//////////////////////////////////////////////////////////////////////////
//!
//! \brief A class stores pathnames for pathIDs
//! \details (e.g. on different processes different paths with same pathID)
//! adding an path autom. changes "\" to "/" and removed last "/" if exists
//!
//////////////////////////////////////////////////////////////////////////

class UbStaticPathMap
{
   typedef std::map< std::string, std::string > PathMap;
public:
   static const std::string GLOBAL;
public:

   static std::string addAndMakePath(const std::string& id, const std::string& path)
   {
      std::string tmpPath = UbStaticPathMap::addPath(id,path);
      if( !tmpPath.empty() ) UbSystem::makeDirectory(tmpPath,20);
      return tmpPath;
   }
   static std::string addPath(const std::string& id, const std::string& path)
   {
      std::string tmpPath = UbSystem::replaceInString(path,"\\","/");
      if(tmpPath.rfind("/") == tmpPath.size()-1) tmpPath.resize(tmpPath.size()-1);
      pathMap[id] = tmpPath;   
      return tmpPath;
   }
   static std::string getPath(const std::string& id)
   {
      PathMap::iterator it = pathMap.find(id);
      if(it == pathMap.end()) return "";
      return it->second;
   }
   static void removePath(const std::string& id)
   {
      pathMap.erase(id);
   }

protected:
   static PathMap pathMap;

private:
   UbStaticPathMap() {}
   UbStaticPathMap(const UbStaticPathMap&) {}
};

#endif //UBSTATICPATHMAP_H

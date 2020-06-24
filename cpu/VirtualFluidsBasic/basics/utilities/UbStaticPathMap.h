//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBSTATICPATHMAP_H
#define UBSTATICPATHMAP_H

#include <iostream>
#include <string>
#include <map>

#include <basics/utilities/UbSystem.h>

/*=========================================================================*/
/*  UbStaticPathMap                                                             */
/*                                                                         */
/**
...
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 12.10.2007
*/ 

/*
stores pathnames for pathIDs (e.g. on different processes different paths with same pathID)
adding an path autom. changes "\" to "/" and removed last "/" if exists

*/

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

#include "RestartPostprocessor.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

//#include "BoostSerializationClassExportHelper.h"

RestartPostprocessor::RestartPostprocessor(RestartObject *object, ArchiveType type)
{
	archiveType = type;
	obj = object;
}
//////////////////////////////////////////////////////////////////////////
RestartPostprocessor::~RestartPostprocessor()
{
}
//////////////////////////////////////////////////////////////////////////
void RestartPostprocessor::doCheckPoint( std::string fname, int step, int myID)
{
   std::string filename = fname + "_Restart_" + UbSystem::toString(myID) + "_" +  UbSystem::toString(step);

   if (archiveType == TXT)
   {
      saveTxtArchive(filename + ".txt");
   } 
   else if(archiveType == BINARY)
   {
      saveBinArchive(filename + ".bin");
   }
}
//////////////////////////////////////////////////////////////////////////
RestartObject* RestartPostprocessor::restart( std::string fname, int step, int myID)
{
   std::string filename = fname + "_Restart_" + UbSystem::toString(myID) + "_" +  UbSystem::toString(step);

   if (archiveType == TXT)
   {
      loadTxtArchive(filename + ".txt");
   } 
   else if(archiveType == BINARY)
   {
      loadBinArchive(filename + ".bin");
   }

   return obj;
}
//////////////////////////////////////////////////////////////////////////
void RestartPostprocessor::saveTxtArchive(std::string filename)
{
   std::ofstream file(filename.c_str()); 
   if(!file)
   { 
      file.clear(); 
      std::string path = UbSystem::getPathFromString(filename);
      if(path.size()>0){ UbSystem::makeDirectory(path); file.open(filename.c_str());}
      if(!file) throw UbException(UB_EXARGS,"couldn't open file "+filename);
   }
   boost::archive::text_oarchive oa(file);
   oa.register_type<RestartObject>();
   oa << obj;
}
//////////////////////////////////////////////////////////////////////////
void RestartPostprocessor::loadTxtArchive( std::string filename )
{
   std::ifstream file(filename.c_str()); 
   //if (!file.is_open()) UB_THROW( UbException(UB_EXARGS,"Can not open check point file \"" + filename + "\""));
   if (!file.is_open()) throw UbException(UB_EXARGS,"couldn't open check point file "+filename);
   printf("the check point file is open...");
   boost::archive::text_iarchive ia(file);
   printf("boost archive is allocated...");
   ia.register_type<RestartObject>();
   printf("register type is set...");
   ia >> obj; 
   printf("obj is filled...");
}
//////////////////////////////////////////////////////////////////////////
void RestartPostprocessor::saveBinArchive( std::string filename )
{
   //impotent for binary archive add std::ios::binary
   std::ofstream file(filename.c_str(), std::ios::binary); 
   if(!file)
   { 
      file.clear(); 
      std::string path = UbSystem::getPathFromString(filename);
      if(path.size()>0){ UbSystem::makeDirectory(path); file.open(filename.c_str());}
      if(!file) throw UbException(UB_EXARGS,"couldn't open file "+filename);
   }
   boost::archive::binary_oarchive oa(file);
   oa.register_type<RestartObject>();
   oa << obj;
}
//////////////////////////////////////////////////////////////////////////
void RestartPostprocessor::loadBinArchive( std::string filename )
{
   //impotent for binary archive add std::ios::binary
   std::ifstream file(filename.c_str(), std::ios::binary); 
   //if (!file.is_open()) UB_THROW( UbException(UB_EXARGS,"Can not open check point file \"" + filename + "\""));
   boost::archive::binary_iarchive ia(file);
   ia.register_type<RestartObject>();
   ia >> obj;
}


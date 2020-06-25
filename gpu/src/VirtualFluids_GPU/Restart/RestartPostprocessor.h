#ifndef Restarter_H
#define Restarter_H

#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include "Restart/RestartObject.h"
#include "basics/utilities/UbSystem.h"

#include <boost/shared_ptr.hpp>
class RestartPostprocessor;
typedef boost::shared_ptr<RestartPostprocessor> RestartPostprocessorPtr;

class RestartPostprocessor
{
public:
   enum ArchiveType {TXT, BINARY};
public:
   RestartPostprocessor(RestartObject *object, ArchiveType typetype = BINARY);
   ~RestartPostprocessor();
   void doCheckPoint( std::string fname, int step, int myID);
   RestartObject* restart( std::string fname, int step, int myID);
protected:
   void saveTxtArchive(std::string filename);
   void loadTxtArchive(std::string filename);
   void saveBinArchive(std::string filename);
   void loadBinArchive(std::string filename);
private:
   ArchiveType archiveType;
   RestartObject *obj;
};

#endif

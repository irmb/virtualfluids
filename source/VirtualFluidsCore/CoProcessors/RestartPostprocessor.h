#ifndef Restarter_H
#define Restarter_H

#include "Postprocessor.h"
#include "Grid3DVisitor.h"
#include "Block3DVisitor.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class RestartPostprocessor;
typedef boost::shared_ptr<RestartPostprocessor> RestartPostprocessorPtr;

class RestartPostprocessor : public Postprocessor
{
public:
   enum ArchiveType {TXT, BINARY};
public:
   RestartPostprocessor(Grid3DPtr& grid, UbSchedulerPtr s, CommunicatorPtr comm, const std::string& path, ArchiveType typetype = BINARY);
   RestartPostprocessor(Grid3DPtr& grid, UbSchedulerPtr s, CommunicatorPtr comm, const std::string& path, int restartStep, ArchiveType typetype = BINARY);
   ~RestartPostprocessor();
   void update(double step);
   void addPostprocessor(PostprocessorPtr p);
   PostprocessorPtr getPostprocessor(int index);
   std::vector<PostprocessorPtr> getPostprocessors();
   void addGridVisitor(Grid3DVisitorPtr v);
   void addBlockVisitor(Block3DVisitorPtr v);
   void doCheckPoint(int step);
   Grid3DPtr restart();
protected:
   void acceptGridVisitors();
   void acceptBlockVisitors();
   void saveTxtArchive(std::string filename);
   void loadTxtArchive(std::string filename);
   void saveBinArchive(std::string filename);
   void loadBinArchive(std::string filename);
   void writeMetafile(int step);
   int readMetafile();
   void checkMetafile();
private:
   std::vector<PostprocessorPtr> postprocessors;
   std::vector<Grid3DVisitorPtr> gridVisitors;
   std::vector<Block3DVisitorPtr> blockVisitors;
   Grid3DPtr grid;
   std::string path;
   ArchiveType archiveType;
   CommunicatorPtr comm;
   int restartStep;
   std::string metafile;
};

#endif

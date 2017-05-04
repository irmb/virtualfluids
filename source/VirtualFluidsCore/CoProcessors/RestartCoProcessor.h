#ifndef Restarter_H
#define Restarter_H

#include "CoProcessor.h"
#include "Grid3DVisitor.h"
#include "Block3DVisitor.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class RestartCoProcessor;
typedef boost::shared_ptr<RestartCoProcessor> RestartCoProcessorPtr;

class RestartCoProcessor : public CoProcessor
{
public:
   enum ArchiveType {TXT, BINARY};
public:
   RestartCoProcessor(Grid3DPtr& grid, UbSchedulerPtr s, CommunicatorPtr comm, const std::string& path, ArchiveType type = BINARY);
   RestartCoProcessor(Grid3DPtr& grid, UbSchedulerPtr s, CommunicatorPtr comm, const std::string& path, int restartStep, ArchiveType type = BINARY);
   ~RestartCoProcessor();
   void process(double step);
   void addCoProcessor(CoProcessorPtr p);
   CoProcessorPtr getCoProcessor(int index);
   std::vector<CoProcessorPtr> getCoProcessors();
   void addGridVisitor(Grid3DVisitorPtr v);
   void addBlockVisitor(Block3DVisitorPtr v);
   void doCheckPoint(int step);
   Grid3DPtr restart();
   void writeDistributedGrid(Grid3DPtr grid, int numberOfProcesses);
   void setArchiveType(ArchiveType type);
   ArchiveType getArchiveType();
protected:
   void acceptGridVisitors();
   void acceptBlockVisitors();
   void saveTxtArchive(std::string filename, Grid3DPtr grid);
   void loadTxtArchive(std::string filename);
   void saveBinArchive(std::string filename, Grid3DPtr grid);
   void loadBinArchive(std::string filename);
   void writeMetafile(int step);
   int readMetafile();
   void checkMetafile();
private:
   std::vector<CoProcessorPtr> CoProcessors;
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

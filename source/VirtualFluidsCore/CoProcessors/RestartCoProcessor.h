#ifndef RESTARTER_H
#define RESTARTER_H

#include <memory>
#include <string>
#include <vector>

#include "CoProcessor.h"


class Grid3D;
class Grid3DVisitor;
class Block3DVisitor;
class Communicator;
class UbScheduler;

class RestartCoProcessor;
typedef std::shared_ptr<RestartCoProcessor> RestartCoProcessorPtr;

class RestartCoProcessor : public CoProcessor
{
public:
   enum ArchiveType {TXT, BINARY};
public:
   RestartCoProcessor(std::shared_ptr<Grid3D>& grid, std::shared_ptr<UbScheduler> s, std::shared_ptr<Communicator> comm, const std::string& path, ArchiveType type = BINARY);
   RestartCoProcessor(std::shared_ptr<Grid3D>& grid, std::shared_ptr<UbScheduler> s, std::shared_ptr<Communicator> comm, const std::string& path, int restartStep, ArchiveType type = BINARY);
   ~RestartCoProcessor();
   void process(double step);
   void addCoProcessor(CoProcessorPtr p);
   CoProcessorPtr getCoProcessor(int index);
   std::vector<CoProcessorPtr> getCoProcessors();
   void addGridVisitor(std::shared_ptr<Grid3DVisitor> v);
   void addBlockVisitor(std::shared_ptr<Block3DVisitor> v);
   void doCheckPoint(int step);
   std::shared_ptr<Grid3D> restart();
   void writeDistributedGrid(std::shared_ptr<Grid3D> grid, int numberOfProcesses);
   void setArchiveType(ArchiveType type);
   ArchiveType getArchiveType();

protected:
   void acceptGridVisitors();
   void acceptBlockVisitors();
   void saveTxtArchive(std::string filename, std::shared_ptr<Grid3D> grid);
   void loadTxtArchive(std::string filename);
   void saveBinArchive(std::string filename, std::shared_ptr<Grid3D> grid);
   void loadBinArchive(std::string filename);
   void writeMetafile(int step);
   int readMetafile();
   void checkMetafile();

private:
   std::vector<CoProcessorPtr> CoProcessors;
   std::vector<std::shared_ptr<Grid3DVisitor> > gridVisitors;
   std::vector<std::shared_ptr<Block3DVisitor>> blockVisitors;
   std::shared_ptr<Grid3D> grid;
   std::string path;
   ArchiveType archiveType;
   std::shared_ptr<Communicator> comm;
   int restartStep;
   std::string metafile;
};

#endif

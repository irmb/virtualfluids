/**
* @file ZoltanPartitioningPatchVisitor.h
* @brief Visitor class wich apply Zoltan library partitioning.
* @author Kostyantyn Kucher
* @date 10.06.2011
*/

#ifndef ZoltanPartitioningGridVisitor_H
#define ZoltanPartitioningGridVisitor_H

#if defined VF_ZOLTAN && defined VF_MPI

#include "Grid3DVisitor.h"
#include "Communicator.h"
#include "ZoltanPartitioner.h"

class ZoltanPartitioningGridVisitor : public Grid3DVisitor
{
public:
   ZoltanPartitioningGridVisitor(CommunicatorPtr comm, int numOfDirs, int numOfLocalParts = 1);
   ~ZoltanPartitioningGridVisitor();
   void visit(Grid3DPtr grid);
protected:
   void collectData(Grid3DPtr grid);
   void repartGrid(Grid3DPtr grid, ZoltanPartitioner& zp);
private:
   CommunicatorPtr comm;
   int numOfDirs;
   int numOfLocalParts;
   ZoltanGraph *graph;
   std::vector<int> vertexGID;
   std::vector<int> numEdges;
   std::vector<int> nborGID;
   std::vector<int> nborProc;
};

#endif
#endif

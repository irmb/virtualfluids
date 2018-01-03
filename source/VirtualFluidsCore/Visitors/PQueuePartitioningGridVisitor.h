/**
* @file PQueuePartitioningPatchVisitor.h
* @brief Visitor class which apply Priority Queue for threads decomposition.
* @author Kostyantyn Kucher
* @date 06.06.2011
*/
#ifndef PQUEUEPARTITIONINGPATCHVISITOR_H 
#define PQUEUEPARTITIONINGPATCHVISITOR_H

#include <memory>

#include "Grid3DVisitor.h"

class Grid3D;

class PQueuePartitioningGridVisitor : public Grid3DVisitor
{
public:
   PQueuePartitioningGridVisitor(int numOfParts);

   void visit(std::shared_ptr<Grid3D> grid) override;

private:
   int numOfParts;
};

#endif

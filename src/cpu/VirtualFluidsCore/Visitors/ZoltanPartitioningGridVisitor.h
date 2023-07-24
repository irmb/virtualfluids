/**
 * @file ZoltanPartitioningPatchVisitor.h
 * @brief Visitor class wich apply Zoltan library partitioning.
 * @author Kostyantyn Kucher
 * @date 10.06.2011
 */

#ifndef ZoltanPartitioningGridVisitor_H
#define ZoltanPartitioningGridVisitor_H

#if defined VF_ZOLTAN && defined VF_MPI

#include <parallel/Communicator.h>
#include "Grid3DVisitor.h"
#include "ZoltanPartitioner.h"

class ZoltanPartitioningGridVisitor : public Grid3DVisitor
{
public:
    ZoltanPartitioningGridVisitor(std::shared_ptr<vf::parallel::Communicator> comm, int numOfDirs, int numOfLocalParts = 1);
    ~ZoltanPartitioningGridVisitor();
    void visit(SPtr<Grid3D> grid);

protected:
    void collectData(SPtr<Grid3D> grid);
    void repartGrid(SPtr<Grid3D> grid, ZoltanPartitioner &zp);

private:
    std::shared_ptr<vf::parallel::Communicator> comm;
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

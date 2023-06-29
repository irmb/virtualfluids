/**
 * @file MetisPartitioner.h
 * @brief Class use METIS library for graph-based partitioning.
 * @author Kostyantyn Kucher
 * @date 22.09.2011
 */

#ifndef METISPARTITIONER_H
#define METISPARTITIONER_H

#if defined VF_METIS

#include "basics/utilities/UbSystem.h"
#include "metis.h"
#include <PointerDefinitions.h>
#include <string>
#include <vector>

class MetisPartitioner
{
public:
    enum PartType { RECURSIVE, KWAY };

public:
    MetisPartitioner();
    idx_t *getMetisOptions();
    void setMetisOptions(int option, idx_t value);
    int partition(int nofParts, PartType ptype);

public:
    std::vector<idx_t> xadj; // adjncy offset of nodes
    //(size = n+1, n=nofNodes)
    std::vector<idx_t> adjncy; // array that stores the adjacency lists of nodes
    //(size = m*2, m= nofEdged, factor 2 because edge A->B AND B->A has to be stored)
    std::vector<idx_t> vwgt;   // vertex weights (size=n*ncon, ncon=nofNodeWeightsPerNode)
    std::vector<idx_t> adjwgt; // array that stores the weights of the adjacency lists
    // (size=2*m)
    idx_t *vsize; // array that stores the computation weights per node
    // (size=n)

    real_t *tpwgts; // holds the wished fratcion of segment i, e.g. tpwgts={0.2, 0.2, 0.6}
    // -> patisions 0 and one will get 20% of the weight each and part 3 60%!
    // (size=nofPartitions)  sum of tpwgts must be 1.0

    real_t *
        ubvec; // This is an array of size ncon that specifies the allowed load imbalance tolerance for each constraint.
    // For the ith partition and jth constraint the allowed weight is the ubvec[j]*tpwgts[i*ncon+j] fraction
    // of the jths constraint total weight. The load imbalances must be greater than 1.0.
    // A NULL value can be passed indicating that the load imbalance tolerance for each constraint should
    // be 1.001 (for ncon=1) or 1.01 (for ncon1).

    std::vector<idx_t>
        part; // This is a vector of size n that upon successful completion stores the partition vector of the graph.
              // The numbering of this vector starts from 0
private:
    idx_t options[METIS_NOPTIONS];
};

#endif

#endif

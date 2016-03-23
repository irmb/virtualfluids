/**
* @file ZoltanPartitioner.h
* @brief Class use Zoltan library for graph-based partitioning.
* @author Kostyantyn Kucher
* @date 10.06.2011
*/

#ifndef ZOLTANPARTITIONER_H
#define ZOLTANPARTITIONER_H

#if defined VF_ZOLTAN && defined VF_MPI

#include "zoltan.h"
#include <vector>
#include <string>

/* Structure to hold graph data */

struct ZoltanGraph{
   int numLocalVertices;        // total vertices in in this partition
   std::vector<int> vvertexGID; // global ID of each of my vertices
   std::vector<int> vnumEdges;  // number of Edges 
   std::vector<int> vnborGID;   // global ID of neighbors
   std::vector<int> vnborProc;  // process owning each nbor in nborGID
};

struct Zoltan_Output{
   ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
   int *importProcs, *importToPart, *exportProcs, *exportToPart;
   int changes, numGidEntries, numLidEntries, numImport, numExport;
};

class ZoltanPartitioner
{
public:
   ZoltanPartitioner(MPI_Comm comm , int rank, int numberOfLocalParts);
   virtual ~ZoltanPartitioner();
   void partition();
   ZoltanGraph* getGraphData();
   void setLB_APPROACH(std::string lb_approach);
   void setNumberOfLocalParts(int numberOfLocalParts);
   void getExportData(std::vector<int>& exportGlobalGids, std::vector<int>& exportToPart, std::vector<int>& exportProcs);
   bool areChanges();

protected:
   static int get_number_of_vertices(void *data, int *ierr);
   static void get_vertex_list(void *data, int sizeGID, int sizeLID,
                  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
   static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                     int num_obj,
                     ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                     int *numEdges, int *ierr);
   static void get_edge_list(void *data, int sizeGID, int sizeLID,
               int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
               int *num_edges,
               ZOLTAN_ID_PTR nborGID, int *nborProc,
               int wgt_dim, float *ewgts, int *ierr);

private:
   MPI_Comm comm;
   int rank;
   int numberOfLocalParts;
   struct Zoltan_Struct *zz;
   std::string lb_approach;
   ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
   int *importProcs, *importToPart, *exportProcs, *exportToPart;
   int changes, numGidEntries, numLidEntries, numImport, numExport;
   ZoltanGraph graph;
};

#endif

#endif 

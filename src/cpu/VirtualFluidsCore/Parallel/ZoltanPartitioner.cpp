#if defined VF_ZOLTAN && defined VF_MPI

#include "ZoltanPartitioner.h"
#include <iostream>
#include <stdlib.h>

#include "UbSystem.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
ZoltanPartitioner::ZoltanPartitioner(MPI_Comm comm, int rank, int numberOfLocalParts): 
                                    comm(comm), rank(rank), numberOfLocalParts(numberOfLocalParts),lb_approach("PARTITION")
{
   int rc;
   float ver;

   rc = Zoltan_Initialize(0, NULL, &ver);

   if (rc != ZOLTAN_OK){
      cout<<"Sorry, Zoltan can't be initialized\n"<<endl;
      MPI_Finalize();
      exit(0);
   } 
   /******************************************************************
   ** Create a Zoltan library structure for this instance of load
   ** balancing.  Set the parameters and query functions that will
   ** govern the library's calculation.  See the Zoltan User's
   ** Guide for the definition of these and many other parameters.
   ******************************************************************/

   zz = Zoltan_Create(comm);
}
//////////////////////////////////////////////////////////////////////////
ZoltanPartitioner::~ZoltanPartitioner()
{
  Zoltan_Destroy(&zz);
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioner::partition()
{
   //General parameters
   Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
   Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
   Zoltan_Set_Param(zz, "LB_APPROACH", lb_approach.c_str());
   Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
   Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
   Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
   string nparts(UbSystem::toString<int>(numberOfLocalParts));
   Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS", nparts.c_str());

   /* Query functions - defined in simpleQueries.h */

   Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &graph);
   Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &graph);
   Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &graph);
   Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &graph);
   
   /******************************************************************
   ** Zoltan can now partition the graph.
   ** In this case, we assume the number of partitions is
   ** equal to the number of processes.  Process rank 0 will own
   ** partition 0, process rank 1 will own partition 1, and so on.
   ******************************************************************/

   int rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
      &changes,        /* 1 if partitioning was changed, 0 otherwise */ 
      &numGidEntries,  /* Number of integers used for a global ID */
      &numLidEntries,  /* Number of integers used for a local ID */
      &numImport,      /* Number of vertices to be sent to me */
      &importGlobalGids,  /* Global IDs of vertices to be sent to me */
      &importLocalGids,   /* Local IDs of vertices to be sent to me */
      &importProcs,    /* Process rank for source of each incoming vertex */
      &importToPart,   /* New partition for each incoming vertex */
      &numExport,      /* Number of vertices I must send to other processes*/
      &exportGlobalGids,  /* Global IDs of the vertices I must send */
      &exportLocalGids,   /* Local IDs of the vertices I must send */
      &exportProcs,    /* Process to which I send each of the vertices */
      &exportToPart);  /* Partition to which each vertex will belong */



   if (rc != ZOLTAN_OK){
      cout << "Partitioning failed on process " << rank <<"\n" << endl;
      MPI_Finalize();
      Zoltan_Destroy(&zz);
      exit(0);
   }
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioner::setLB_APPROACH(std::string lb_approach)
{
   this->lb_approach = lb_approach;
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioner::setNumberOfLocalParts(int numberOfLocalParts)
{
   this->numberOfLocalParts = numberOfLocalParts;
}
//////////////////////////////////////////////////////////////////////////
// Application defined query functions //
//////////////////////////////////////////////////////////////////////////
int ZoltanPartitioner::get_number_of_vertices(void *data, int *ierr)
{
   ZoltanGraph *graph = (ZoltanGraph *)data;
   *ierr = ZOLTAN_OK;
   return graph->numLocalVertices;
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioner::get_vertex_list(void *data, int sizeGID, int sizeLID,
                                        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                        int wgt_dim, float *obj_wgts, int *ierr)
{
   ZoltanGraph *graph = (ZoltanGraph *)data;
   *ierr = ZOLTAN_OK;

   /* In this case, return the IDs of our vertices, but no weights.
   * Zoltan will assume equally weighted vertices.
   */

   for (int i=0; i<graph->numLocalVertices; i++){
      globalID[i] = graph->vvertexGID[i];
      localID[i] = i;
   }
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioner::get_num_edges_list(void *data, int sizeGID, int sizeLID,
                                           int num_obj,
                                           ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                           int *numEdges, int *ierr)
{
   ZoltanGraph *graph = (ZoltanGraph *)data;

   if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->numLocalVertices)){
      *ierr = ZOLTAN_FATAL;
      return;
   }

   for (int i=0;  i < num_obj ; i++){
      numEdges[i] = graph->vnumEdges[i];
   }

   *ierr = ZOLTAN_OK;
   return;
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioner::get_edge_list(void *data, int sizeGID, int sizeLID,
                                      int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                      int *num_edges,
                                      ZOLTAN_ID_PTR nborGID, int *nborProc,
                                      int wgt_dim, float *ewgts, int *ierr)
{
   int *nextNbor, *nextProc;

   ZoltanGraph *graph = (ZoltanGraph *)data;
   *ierr = ZOLTAN_OK;

   if ( (sizeGID != 1) || (sizeLID != 1) || 
      (num_obj != graph->numLocalVertices)||
      (wgt_dim != 0)){
         *ierr = ZOLTAN_FATAL;
         return;
   }

   nextNbor = (int *)nborGID;
   nextProc = nborProc;
   
   int n=0;
   for (int i=0; i < num_obj; i++){

      /*
      * In this case, we are not setting edge weights.  Zoltan will
      * set each edge to weight 1.0.
      */

      for (int j=0; j < num_edges[i]; j++){
         nborGID[n] = graph->vnborGID[n];
         nborProc[n] = graph->vnborProc[n];
         n++;
      }
   }
   return;
}
//////////////////////////////////////////////////////////////////////////
ZoltanGraph* ZoltanPartitioner::getGraphData()
{
   return &graph;
}
//////////////////////////////////////////////////////////////////////////
void ZoltanPartitioner::getExportData(vector<int>& exportGlobalGids, vector<int>& exportToPart, vector<int>& exportProcs)
{
   for (int i = 0; i < this->numExport; i++)
   {
      exportGlobalGids.push_back(static_cast<int> (this->exportGlobalGids[i]));
      exportToPart.push_back(this->exportToPart[i]);
      exportProcs.push_back(this->exportProcs[i]);
   }
}
//////////////////////////////////////////////////////////////////////////
 bool ZoltanPartitioner::areChanges()
 {
     return static_cast<bool>(this->changes);
 }
#endif

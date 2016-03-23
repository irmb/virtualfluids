#if defined VF_METIS

#include "MetisPartitioner.h"


MetisPartitioner::MetisPartitioner()
{
   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_NUMBERING] = 0;
}
//////////////////////////////////////////////////////////////////////////
MetisPartitioner::~MetisPartitioner()
{

}
//////////////////////////////////////////////////////////////////////////
idx_t* MetisPartitioner::getMetisOptions()
{
   return options;
}
//////////////////////////////////////////////////////////////////////////
int MetisPartitioner::partition(int nofParts, MetisPartitioner::PartType ptype)
{
   int rc;
   idx_t nvtxs = (idx_t)xadj.size()-1;  // number of nodes
   idx_t ncon = (idx_t)vwgt.size()/nvtxs; // number Of node constraints;

   part.resize(nvtxs);
   int edgecutCount = 0;

   switch (ptype)
   {
   case MetisPartitioner::RECURSIVE: 
      if     ( nofParts <  1 ) UB_THROW( UbException(UB_EXARGS,"invalid nofParts<1") );
      else if (nofParts == 1) { part.resize(nvtxs, 0); return 0; }
      //else if( nofParts >  8 ) UBLOG(logWARNING, "MetisPartitioner::Recursive: !!!Warning!!!  best for nofParts<=8 --> Kway is maybe a better option");
      
      rc = METIS_PartGraphRecursive(&nvtxs, &ncon, &xadj[0], &adjncy[0],
                                    &vwgt[0], &vsize[0], &adjwgt[0], &nofParts, 
                                    &tpwgts[0], &ubvec[0], options, &edgecutCount, &part[0]);
   	break;
   case MetisPartitioner::KWAY: 
      if     ( nofParts <  1 ) UB_THROW( UbException(UB_EXARGS,"invalid nofParts<1") );
      else if (nofParts == 1) { part.resize(nvtxs, 0); return 0; }
      //else if( nofParts <  9 ) UBLOG(logWARNING, "MetisPartitioner::Kway: !!!Warning!!!  best for nofParts>8 --> Recursive is maybe a better option");

      rc = METIS_PartGraphKway(&nvtxs, &ncon, &xadj[0], &adjncy[0],
                                &vwgt[0], &vsize[0], &adjwgt[0], &nofParts,
                                &tpwgts[0], &ubvec[0], options, &edgecutCount, &part[0]);
      break;
   }

   switch (rc)
   {
   case METIS_ERROR_INPUT:
      throw UbException(UB_EXARGS,"METIS: input error");
   	break;
   case METIS_ERROR_MEMORY:
      throw UbException(UB_EXARGS,"METIS: it could not allocate the required memory");
      break;
   case METIS_ERROR:
      throw UbException(UB_EXARGS,"METIS: error");
      break;
   }

   return edgecutCount;
}
#endif

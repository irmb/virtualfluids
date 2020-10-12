#if defined VF_METIS

#include "MetisPartitioner.h"


MetisPartitioner::MetisPartitioner()
{
   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_NUMBERING] = 0;
   vsize = NULL;
   tpwgts = NULL;
   ubvec = NULL;

   //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
   ////options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

   //options[METIS_OPTION_CTYPE]  = METIS_CTYPE_SHEM;
   //options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
}
//////////////////////////////////////////////////////////////////////////
MetisPartitioner::~MetisPartitioner()
= default;
//////////////////////////////////////////////////////////////////////////
idx_t* MetisPartitioner::getMetisOptions()
{
   return options;
}
void MetisPartitioner::setMetisOptions(int option, idx_t value)
{
   options[option] = value;
}
//////////////////////////////////////////////////////////////////////////
int MetisPartitioner::partition(int nofParts, MetisPartitioner::PartType ptype)
{
   int rc;
   idx_t nvtxs = (idx_t)xadj.size()-1;  // number of nodes
   idx_t ncon = (idx_t)vwgt.size()/nvtxs; // number Of node constraints;
   part.resize(nvtxs);
   idx_t edgecutCount = 0;
   idx_t nofPartsMetis = (idx_t)nofParts;

   switch (ptype)
   {
   case MetisPartitioner::RECURSIVE: 
      if     ( nofParts <  1 ) UB_THROW( UbException(UB_EXARGS,"invalid nofParts<1") );
      else if (nofParts == 1) { part.resize(nvtxs, 0); return 0; }
      //else if( nofParts >  8 ) UBLOG(logWARNING, "MetisPartitioner::Recursive: !!!Warning!!!  best for nofParts<=8 --> Kway is maybe a better option");
      
      rc = METIS_PartGraphRecursive(&nvtxs, &ncon, &xadj[0], &adjncy[0],
                                    &vwgt[0], vsize, &adjwgt[0], &nofPartsMetis, 
                                    tpwgts, ubvec, options, &edgecutCount, &part[0]);
   	break;
   case MetisPartitioner::KWAY: 
      if     ( nofParts <  1 ) UB_THROW( UbException(UB_EXARGS,"invalid nofParts<1") );
      else if (nofParts == 1) { part.resize(nvtxs, 0); return 0; }
      //else if( nofParts <  9 ) UBLOG(logWARNING, "MetisPartitioner::Kway: !!!Warning!!!  best for nofParts>8 --> Recursive is maybe a better option");

      rc = METIS_PartGraphKway(&nvtxs, &ncon, &xadj[0], &adjncy[0],
                                &vwgt[0], vsize, &adjwgt[0], &nofPartsMetis,
                                tpwgts, ubvec, options, &edgecutCount, &part[0]);
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

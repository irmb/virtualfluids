#ifndef CoarsenCrossAndInsideGbObjectBlockVisitor_H
#define CoarsenCrossAndInsideGbObjectBlockVisitor_H

#include "Block3DVisitor.h"
#include <numerics/geometry3d/GbObject3D.h>

//! \brief Refine blocks on base of bounding box which is defined with <i>geoObject</i>
//! \details The class uses a geometry object for define a bounding box. Inside and across this bounding box will be grid on block basis refinement.
//! \author K. Kutscher
class CoarsenCrossAndInsideGbObjectBlockVisitor : public Block3DVisitor
{
public:
   //! A default constructor
   CoarsenCrossAndInsideGbObjectBlockVisitor();
   //! A constructor
   //! \param geoObject a smart pointer to bounding box
   //! \param refineLevel an integer for refine on this level
   CoarsenCrossAndInsideGbObjectBlockVisitor(GbObject3DPtr geoObject, int fineLevel, int coarseLevel);
   virtual ~CoarsenCrossAndInsideGbObjectBlockVisitor();
   void visit(Grid3DPtr grid, Block3DPtr block);
   //////////////////////////////////////////////////////////////////////////
protected:
   GbObject3DPtr geoObject;
   bool notActive;
   int coarseLevel;
};

#endif 

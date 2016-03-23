#ifndef RefineCrossAndInsideGbObjectBlockVisitor_H
#define RefineCrossAndInsideGbObjectBlockVisitor_H

#include "Block3DVisitor.h"
#include <numerics/geometry3d/GbObject3D.h>

//! \brief Refine blocks on base of bounding box which is defined with <i>geoObject</i>
//! \details The class uses a geometry object for define a bounding box. Inside and across this bounding box will be grid on block basis refinement.
//! \author K. Kucher
class RefineCrossAndInsideGbObjectBlockVisitor : public Block3DVisitor
{
public:
   //! A default constructor
   RefineCrossAndInsideGbObjectBlockVisitor();
   //! A constructor
   //! \param geoObject a smart pointer to bounding box
   //! \param refineLevel an integer for refine on this level
   RefineCrossAndInsideGbObjectBlockVisitor(GbObject3DPtr geoObject, int refineLevel);
   virtual ~RefineCrossAndInsideGbObjectBlockVisitor();
   void visit(Grid3DPtr grid, Block3DPtr block);
   //////////////////////////////////////////////////////////////////////////
protected:
   GbObject3DPtr geoObject;
   bool notActive;
};

#endif 

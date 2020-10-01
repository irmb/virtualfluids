#ifndef CoarsenCrossAndInsideGbObjectBlockVisitor_H
#define CoarsenCrossAndInsideGbObjectBlockVisitor_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class GbObject3D;
class Block3D;
class Grid3D;

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
   CoarsenCrossAndInsideGbObjectBlockVisitor(SPtr<GbObject3D> geoObject, int fineLevel, int coarseLevel);
   virtual ~CoarsenCrossAndInsideGbObjectBlockVisitor();
      void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
   //////////////////////////////////////////////////////////////////////////
protected:
    SPtr<GbObject3D> geoObject;
   bool notActive;
   int coarseLevel;
};

#endif 

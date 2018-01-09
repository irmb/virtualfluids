#ifndef RefineCrossAndInsideGbObjectBlockVisitor_H
#define RefineCrossAndInsideGbObjectBlockVisitor_H

#include <vector>
#include <memory>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class GbObject3D;

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
   RefineCrossAndInsideGbObjectBlockVisitor(std::shared_ptr<GbObject3D> geoObject, int refineLevel);
   virtual ~RefineCrossAndInsideGbObjectBlockVisitor();

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;
protected:
    std::shared_ptr<GbObject3D> geoObject;
   bool notActive;
};

#endif 

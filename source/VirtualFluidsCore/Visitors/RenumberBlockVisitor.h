/**
* @file RenumberBlockVisitor.h
* @brief Visitor class which renumber blocks.
* @author Konstantin Kutscher
* @date 06.06.2011
*/

#ifndef RenumberBlockVisitor_h
#define RenumberBlockVisitor_h

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;

//! \brief  Visitor class which renumber blocks.
//! \details Visitor class which renumber blocks.            
//! \author  Konstantin Kutscher 
class RenumberBlockVisitor : public Block3DVisitor
{
public:
   RenumberBlockVisitor();

   virtual ~RenumberBlockVisitor() {}

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

private:
   static int counter;
};

#endif

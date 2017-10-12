/**
* @file RenumberBlockVisitor.h
* @brief Visitor class which renumber blocks.
* @author Konstantin Kutscher
* @date 06.06.2011
*/

#ifndef RenumberBlockVisitor_h
#define RenumberBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel.h"

//! \brief  Visitor class which renumber blocks.
//! \details Visitor class which renumber blocks.            
//! \author  Konstantin Kutscher 
class RenumberBlockVisitor : public Block3DVisitor
{
public:
   RenumberBlockVisitor();

   virtual ~RenumberBlockVisitor() {}

   virtual void visit(Grid3DPtr grid, Block3DPtr block);

private:
   static int counter;
};

#endif

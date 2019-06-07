/**
* @file RenumberGridVisitor.h
* @brief Visitor class which renumber blocks.
* @author Konstantin Kutscher
* @date 06.06.2011
*/

#ifndef RenumberGridVisitor_h
#define RenumberGridVisitor_h

#include "Grid3DVisitor.h"
#include "Communicator.h"

class Grid3D;

//! \brief  Visitor class which renumber blocks.
//! \details Visitor class which renumber blocks.            
//! \author  Konstantin Kutscher 
class RenumberGridVisitor : public Grid3DVisitor
{
public:
   RenumberGridVisitor(SPtr<Communicator> com);

   virtual ~RenumberGridVisitor() {}

   void visit(SPtr<Grid3D> grid) override;

private:
    SPtr<Communicator> comm;
//   static int counter;
};

#endif

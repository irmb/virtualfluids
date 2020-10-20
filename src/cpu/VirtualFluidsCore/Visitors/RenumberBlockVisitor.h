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

    ~RenumberBlockVisitor() override = default;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    static int counter;
};

#endif

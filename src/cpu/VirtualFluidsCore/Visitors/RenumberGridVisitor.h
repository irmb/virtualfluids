/**
 * @file RenumberGridVisitor.h
 * @brief Visitor class which renumber blocks.
 * @author Konstantin Kutscher
 * @date 06.06.2011
 */

#ifndef RenumberGridVisitor_h
#define RenumberGridVisitor_h

#include <mpi/Communicator.h>
#include "Grid3DVisitor.h"

class Grid3D;

//! \brief  Visitor class which renumber blocks in order: rank->level.
//! \details Visitor class which renumber blocks.
//! \author  Konstantin Kutscher
class RenumberGridVisitor : public Grid3DVisitor
{
public:
    RenumberGridVisitor(std::shared_ptr<vf::mpi::Communicator> com);

    ~RenumberGridVisitor() override = default;

    void visit(SPtr<Grid3D> grid) override;

private:
    std::shared_ptr<vf::mpi::Communicator> comm;
    //   static int counter;
};

#endif

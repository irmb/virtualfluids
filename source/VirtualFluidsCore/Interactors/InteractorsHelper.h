#ifndef SolidBlocksHelper_h 
#define SolidBlocksHelper_h

#include <vector>
#include <PointerDefinitions.h>
#include <SetSolidOrBoundaryBlockVisitor.h>

class Interactor3D;
class Block3D;
class Grid3D;
class Grid3DVisitor;
class SetSolidOrBoundaryBlockVisitor;

class InteractorsHelper
{
public:
   InteractorsHelper(SPtr<Grid3D> grid, SPtr<Grid3DVisitor> visitor);
   ~InteractorsHelper();

   void addInteractor(SPtr<Interactor3D> interactor);
   void selectBlocks();
   void setBC();
   void sendDomainDecompositionVisitor() const;

protected:
   void deleteSolidBlocks();
   void setBlocks(const SPtr<Interactor3D> interactor, SetSolidOrBoundaryBlockVisitor::BlockType type) const;
   void setBcBlocks();

private:
   void updateGrid();

   std::vector<SPtr<Interactor3D> > interactors;
   SPtr<Grid3D> grid;
   std::vector<SPtr<Block3D> > solidBlocks;
   SPtr<Grid3DVisitor> visitor;
};

#endif

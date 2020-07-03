#ifndef SolidBlocksHelper_h 
#define SolidBlocksHelper_h

#include <vector>
#include <PointerDefinitions.h>


class Interactor3D;
class Block3D;
class Grid3D;
class Grid3DVisitor;

class InteractorsHelper
{
public:
   InteractorsHelper(SPtr<Grid3D> grid, SPtr<Grid3DVisitor> visitor, bool deleteBlocks=true);
   ~InteractorsHelper();

   void addInteractor(SPtr<Interactor3D> interactor);
   void selectBlocks();
   void setBC();
   void sendDomainDecompositionVisitor() const;

protected:
   void deleteSolidBlocks();
   void setBcBlocks();

private:
   void updateGrid();

   std::vector<SPtr<Interactor3D> > interactors;
   SPtr<Grid3D> grid;
   std::vector<SPtr<Block3D> > solidBlocks;
   SPtr<Grid3DVisitor> visitor;
   bool deleteBlocks;
};

#endif

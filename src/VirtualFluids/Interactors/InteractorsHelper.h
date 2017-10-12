#ifndef SolidBlocksHelper_h 
#define SolidBlocksHelper_h

#include <Grid3D.h>
#include <Communicator.h>
#include <Interactor3D.h>

class InteractorsHelper
{
public:
   InteractorsHelper(Grid3DPtr grid, Grid3DVisitorPtr visitor);
   ~InteractorsHelper();
   void addInteractor(Interactor3DPtr interactor);
   void selectBlocks();
   void setBC();
protected:
   void deleteSolidBlocks();
   void setTransBlocks();
private:
   void updateGrid();
   std::vector<Interactor3DPtr> interactors;
   Grid3DPtr grid;
   std::vector<Block3DPtr> solidBlocks;
   Grid3DVisitorPtr visitor;
};

#endif

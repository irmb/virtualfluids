#ifndef SolidBlocksHelper_h 
#define SolidBlocksHelper_h

#include <vector>
#include <memory>

class Interactor3D;
class Block3D;
class Grid3D;
class Grid3DVisitor;
enum class BlockType;

class InteractorsHelper
{
public:
   InteractorsHelper(std::shared_ptr<Grid3D> grid, std::shared_ptr<Grid3DVisitor> visitor);
   ~InteractorsHelper();

   void addInteractor(std::shared_ptr<Interactor3D> interactor);
   void selectBlocks();
   void setBC();
    void sendDomainDecompositionVisitor() const;

protected:
   void deleteSolidBlocks();
    void setBlocks(const std::shared_ptr<Interactor3D> interactor, BlockType type) const;
    void setBcBlocks();

private:
   void updateGrid();

   std::vector<std::shared_ptr<Interactor3D> > interactors;
   std::shared_ptr<Grid3D> grid;
   std::vector<std::shared_ptr<Block3D> > solidBlocks;
   std::shared_ptr<Grid3DVisitor> visitor;
};

#endif

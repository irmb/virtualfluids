#ifndef RefineInterGbObjectsVisirtor_H
#define RefineInterGbObjectsVisirtor_H

#include <vector>
#include <memory>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class GbObject3D;

//////////////////////////////////////////////////////////////////////////
class RefineInterGbObjectsBlockVisitor : public Block3DVisitor
{
public:
   RefineInterGbObjectsBlockVisitor();
   RefineInterGbObjectsBlockVisitor(std::shared_ptr<GbObject3D> includeGbObject3D, std::shared_ptr<GbObject3D> excludeGbObject3D, int startlevel, int stoplevel);
   RefineInterGbObjectsBlockVisitor(std::vector<std::shared_ptr<GbObject3D> > includeGbObjects3D, std::vector<std::shared_ptr<GbObject3D> > excludeGbObjects3D, int startlevel, int stoplevel);
   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

private:
   std::vector<std::shared_ptr<GbObject3D> > includeGbObjects3D;
   std::vector<std::shared_ptr<GbObject3D> > excludeGbObjects3D;
};

#endif //RefineInterGbObjectsVisirtor_H

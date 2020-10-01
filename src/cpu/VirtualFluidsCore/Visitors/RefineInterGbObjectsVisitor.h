#ifndef RefineInterGbObjectsVisirtor_H
#define RefineInterGbObjectsVisirtor_H

#include <vector>
#include <PointerDefinitions.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class GbObject3D;

//////////////////////////////////////////////////////////////////////////
class RefineInterGbObjectsBlockVisitor : public Block3DVisitor
{
public:
   RefineInterGbObjectsBlockVisitor();
   RefineInterGbObjectsBlockVisitor(SPtr<GbObject3D> includeGbObject3D, SPtr<GbObject3D> excludeGbObject3D, int startlevel, int stoplevel);
   RefineInterGbObjectsBlockVisitor(std::vector<SPtr<GbObject3D> > includeGbObjects3D, std::vector<SPtr<GbObject3D> > excludeGbObjects3D, int startlevel, int stoplevel);
   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
   std::vector<SPtr<GbObject3D> > includeGbObjects3D;
   std::vector<SPtr<GbObject3D> > excludeGbObjects3D;
};

#endif //RefineInterGbObjectsVisirtor_H

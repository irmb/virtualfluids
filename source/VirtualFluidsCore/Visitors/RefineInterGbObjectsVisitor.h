#ifndef RefineInterGbObjectsVisirtor_H
#define RefineInterGbObjectsVisirtor_H

#include <vector>

#include "Block3DVisitor.h"
#include <numerics/geometry3d/GbObject3D.h>

//////////////////////////////////////////////////////////////////////////
class RefineInterGbObjectsBlockVisitor : public Block3DVisitor
{
public:
   RefineInterGbObjectsBlockVisitor();
   RefineInterGbObjectsBlockVisitor(GbObject3DPtr includeGbObject3D, GbObject3DPtr excludeGbObject3D, int startlevel, int stoplevel);
   RefineInterGbObjectsBlockVisitor(std::vector<GbObject3DPtr> includeGbObjects3D, std::vector<GbObject3DPtr> excludeGbObjects3D, int startlevel, int stoplevel);
   void visit(Grid3DPtr grid, Block3DPtr block);

private:
   std::vector<GbObject3DPtr> includeGbObjects3D;
   std::vector<GbObject3DPtr> excludeGbObjects3D;
};

#endif //RefineInterGbObjectsVisirtor_H

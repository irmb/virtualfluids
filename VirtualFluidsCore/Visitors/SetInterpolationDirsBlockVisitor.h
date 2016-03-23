#ifndef SetInterpolationDirsBlockVisitor_h
#define SetInterpolationDirsBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel3D.h"


class SetInterpolationDirsBlockVisitor : public Block3DVisitor
{
public:
   SetInterpolationDirsBlockVisitor(std::vector<int>& dirs);

   virtual ~SetInterpolationDirsBlockVisitor() {}

   virtual void visit(Grid3DPtr grid, Block3DPtr block);

private:
   std::vector<int> dirs;
   void checkFlagDir(Grid3DPtr grid, int dir1, int dir2, bool &flagDirection, int ix1, int ix2, int ix3, int level);
   void checkFlagDir(Grid3DPtr grid, int dir1, int dir2, int dir3, bool &flagDirection, int ix1, int ix2, int ix3, int level);
};

#endif

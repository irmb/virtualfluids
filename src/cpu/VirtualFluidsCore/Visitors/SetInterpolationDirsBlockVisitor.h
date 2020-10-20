#ifndef SetInterpolationDirsBlockVisitor_h
#define SetInterpolationDirsBlockVisitor_h

#include <PointerDefinitions.h>
#include <vector>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;

class SetInterpolationDirsBlockVisitor : public Block3DVisitor
{
public:
    SetInterpolationDirsBlockVisitor(std::vector<int> &dirs);

    ~SetInterpolationDirsBlockVisitor() override = default;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    std::vector<int> dirs;
    void checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, bool &flagDirection, int ix1, int ix2, int ix3, int level);
    void checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, int dir3, bool &flagDirection, int ix1, int ix2, int ix3,
                      int level);
};

#endif

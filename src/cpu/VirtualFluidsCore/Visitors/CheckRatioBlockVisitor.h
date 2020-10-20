#ifndef CheckRatioBlockVisitor_H
#define CheckRatioBlockVisitor_H

#include <PointerDefinitions.h>
#include <string>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;

class CheckRatioBlockVisitor : public Block3DVisitor
{
public:
    CheckRatioBlockVisitor(int levelDepth, bool includeNotActiveBlocks = true);

    ~CheckRatioBlockVisitor() override = default;

    bool getState();
    void resetState();
    std::string getStateString();

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    int levelDepth;
    bool includeNotActiveBlocks;
    bool state;
    SPtr<Block3D> falseBlock;
};

#endif // OverlapBlockVisitor_H

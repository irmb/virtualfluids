#ifndef InitDistributionsFromFileBlockVisitor_h__
#define InitDistributionsFromFileBlockVisitor_h__

#include "Block3DVisitor.h"

class InitDistributionsFromFileBlockVisitor : public Block3DVisitor
{
public:
   InitDistributionsFromFileBlockVisitor(LBMReal nu, LBMReal rho, std::string file);
   ~InitDistributionsFromFileBlockVisitor();
   void visit(Grid3DPtr grid, Block3DPtr block);
private:
   CbArray4D<LBMReal, IndexerX4X3X2X1> matrix;
   enum Velocity { Vx1, Vx2, Vx3 };
   LBMReal nu;
   LBMReal rho;
};
#endif // InitDistributionsFromFileBlockVisitor_h__



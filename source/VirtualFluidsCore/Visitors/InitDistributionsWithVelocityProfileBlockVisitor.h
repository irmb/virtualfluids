//#ifndef InitDistributionsWithVelocityProfileBlockVisitor_h__
//#define InitDistributionsWithVelocityProfileBlockVisitor_h__
//
//#include "Block3DVisitor.h"
//
//class InitDistributionsWithVelocityProfileBlockVisitor : public Block3DVisitor
//{
//public:
//   InitDistributionsWithVelocityProfileBlockVisitor(LBMReal nu, LBMReal rho, double dx, std::string file);
//   ~InitDistributionsWithVelocityProfileBlockVisitor();
//   void visit(Grid3DPtr grid, Block3DPtr block);
//private:
//   std::vector<double> z;
//   std::vector<double> vx1;
//   std::vector<double> vx2;
//   std::vector<double> vx3;
//
//   LBMReal nu;
//   LBMReal rho;
//
//   Grid3DPtr oldGrid;
//};
//
//#endif // InitDistributionsWithVelocityProfileBlockVisitor_h__

#ifndef D3Q27INTRPOLATIOPROCESSOR_H_
#define D3Q27INTRPOLATIOPROCESSOR_H_

#include "InterpolationProcessor.h"
#include "LBMSystem.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"
#include "BCArray3D.h"


struct D3Q27ICell
{
   LBMReal TSW[27]; 
   LBMReal TNW[27];
   LBMReal TNE[27];
   LBMReal TSE[27];
   LBMReal BSW[27];
   LBMReal BNW[27];
   LBMReal BNE[27];
   LBMReal BSE[27];
};

class InterpolationProcessor;
typedef boost::shared_ptr<InterpolationProcessor> InterpolationProcessorPtr;

#include "InterpolationHelper.h"

class InterpolationProcessor
{
public:
   InterpolationProcessor();
   virtual ~InterpolationProcessor();
   virtual InterpolationProcessorPtr clone() = 0;
   virtual void setOmegas(LBMReal omegaC, LBMReal omegaF) = 0;
   virtual void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF) = 0;
   virtual void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, LBMReal xoff, LBMReal yoff, LBMReal zoff) = 0;
   virtual void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC) = 0; 
   virtual void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC, LBMReal xoff, LBMReal yoff, LBMReal zoff) = 0; 

   static void readICell(DistributionArray3DPtr f, D3Q27ICell& icell, int x1, int x2, int x3);
   static void writeICell(DistributionArray3DPtr f, const D3Q27ICell& icell, int x1, int x2, int x3);
   static void writeICellInv(DistributionArray3DPtr f, const D3Q27ICell& icell, int x1, int x2, int x3);
   static void writeINode(DistributionArray3DPtr f, const LBMReal* const inode, int x1, int x2, int x3);
   static void writeINodeInv(DistributionArray3DPtr f, const LBMReal* const inode, int x1, int x2, int x3);
   static bool iCellHasSolid(const BCArray3DPtr bcArray, int x1, int x2, int x3);
   static int  iCellHowManySolids(const BCArray3DPtr bcArray, int x1, int x2, int x3);

   bool findNeighborICell(const BCArray3DPtr bcArray, DistributionArray3DPtr f, 
                          D3Q27ICell& icell, int maxX1, int maxX2, int maxX3, 
                          int x1, int x2, int x3, LBMReal& xoff, LBMReal& yoff, LBMReal& zoff);

protected:
   virtual void calcInterpolatedCoefficiets(const D3Q27ICell& icell, LBMReal omega, LBMReal eps_new){}
   virtual void calcInterpolatedNodeFC(LBMReal* f, LBMReal omega){}
   virtual void calcInterpolatedVelocity(LBMReal x, LBMReal y, LBMReal z,LBMReal& vx1, LBMReal& vx2, LBMReal& vx3){}
   virtual void calcInterpolatedShearStress(LBMReal x, LBMReal y, LBMReal z,LBMReal& tauxx, LBMReal& tauyy, LBMReal& tauzz,LBMReal& tauxy, LBMReal& tauxz, LBMReal& tauyz){}
   virtual void setOffsets(LBMReal xoff, LBMReal yoff, LBMReal zoff){}
   friend class InterpolationHelper;
private:
   bool inRange(int x1, int x2, int x3);
   int m_maxX1, m_maxX2, m_maxX3;
};

//////////////////////////////////////////////////////////////////////////
inline bool InterpolationProcessor::inRange(int x1, int x2, int x3)
{
   return x1 >= 0 && x1 < m_maxX1 && x2 >= 0 && x2 < m_maxX2 && x3 >= 0 && x3 < m_maxX3;
}

#endif

#ifndef InitDistributionsBlockVisitor_H
#define InitDistributionsBlockVisitor_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "D3Q27System.h"

#include <muParser.h>

/*================================================================================*/
/*  D3Q27ETInitDistributionsBlockVisitor                                             */
/*                                                                                */
/**
more flexible way to initialize flow area
you can define functions to calculate macroscopic values for feq 
!!! x1,x2,x3 are automatically defined via this adapter and are the real world
vertex coordinates !!!

if function is invalid an UbException with detailed information is thrown

<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 19.04.08
*/ 

//! \details example:<BR>
//! D3Q27InitDistributionsBlockVisitor init(1.0,0.0,0.0,0.0);<BR>
//! Bem.: rho=0.0 bei inkompressibel<BR>
//! init.setVx1("0.01*x2");<BR>
//! init.setVx1("0.01*x2^2");<BR>
//! patch.adaptByPatchCriterion(init);

class Grid3D;
class Block3D;


class InitDistributionsBlockVisitor : public Block3DVisitor
{
public:
   typedef std::numeric_limits<LBMReal> D3Q27RealLim;

public:
   InitDistributionsBlockVisitor();
   //////////////////////////////////////////////////////////////////////////
   //automatic vars are: x1,x2, x3
   //ussage example: setVx1("x1*0.01+x2*0.003")
   //////////////////////////////////////////////////////////////////////////
   void setVx1( const mu::Parser& parser);
   void setVx2( const mu::Parser& parser);
   void setVx3( const mu::Parser& parser);
   void setRho( const mu::Parser& parser);

   void setVx1( const std::string& muParserString);
   void setVx2( const std::string& muParserString);
   void setVx3( const std::string& muParserString);
   void setRho( const std::string& muParserString);
   //////////////////////////////////////////////////////////////////////////
   void setVx1( LBMReal vx1 );
   void setVx2( LBMReal vx2 );
   void setVx3( LBMReal vx3 );
   void setRho( LBMReal rho );

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

protected:
   void checkFunction(mu::Parser fct);

private:
   mu::Parser muVx1;
   mu::Parser muVx2;
   mu::Parser muVx3;
   mu::Parser muRho;
};

#endif //D3Q27INITDISTRIBUTIONSPATCHVISITOR_H

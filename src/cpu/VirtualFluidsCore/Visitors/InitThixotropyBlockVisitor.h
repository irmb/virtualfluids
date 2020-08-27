#ifndef InitThixotropyBlockVisitor_H
#define InitThixotropyBlockVisitor_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "D3Q27System.h"

#include <MuParser/include/muParser.h>

/*================================================================================*/
/*  D3Q27ETInitThixotropyBlockVisitor                                             */
/*                                                                                */




class InitThixotropyBlockVisitor : public Block3DVisitor
{
public:
	typedef std::numeric_limits<LBMReal> D3Q27RealLim;

public:
	InitThixotropyBlockVisitor();
	//D3Q27ETInitThixotropyBlockVisitor(LBMReal rho, LBMReal vx1=0.0, LBMReal vx2=0.0, LBMReal vx3=0.0);
	//! Constructor
	//! \param nu - viscosity
	//! \param rho - density
	//! \param vx1 - velocity in x
	//! \param vx2 - velocity in y
	//! \param vx3 - velocity in z
	//InitThixotropyBlockVisitor(LBMReal lambda /*LBMReal nu, LBMReal D, LBMReal rho, LBMReal vx1 = 0.0, LBMReal vx2 = 0.0, LBMReal vx3 = 0.0, LBMReal c=0.0, LBMReal f1 = 0.0, LBMReal f2 = 0.0, LBMReal f3 = 0.0*/);
	//////////////////////////////////////////////////////////////////////////
	//automatic vars are: x1,x2, x3
	//ussage example: setVx1("x1*0.01+x2*0.003")
	//////////////////////////////////////////////////////////////////////////
	//void setVx1(const mu::Parser& parser);
	//void setVx2(const mu::Parser& parser);
	//void setVx3(const mu::Parser& parser);
	//void setRho(const mu::Parser& parser);

	//void setVx1(const std::string& muParserString);
	//void setVx2(const std::string& muParserString);
	//void setVx3(const std::string& muParserString);
	//void setRho(const std::string& muParserString);
	////////////////////////////////////////////////////////////////////////////
	//void setVx1(LBMReal vx1);
	//void setVx2(LBMReal vx2);
	//void setVx3(LBMReal vx3);
	//void setRho(LBMReal rho);
	//void setNu(LBMReal nu);

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//void setf1(const mu::Parser& parser);
	//void setf2(const mu::Parser& parser);
	//void setf3(const mu::Parser& parser);
	void setLambda(const mu::Parser& parser);

	//void setf1(const std::string& muParserString);
	//void setf2(const std::string& muParserString);
	//void setf3(const std::string& muParserString);
	void setLambda(const std::string& muParserString);
	//////////////////////////////////////////////////////////////////////////
	//void setf1(LBMReal f1);
	//void setf2(LBMReal f2);
	//void setf3(LBMReal f3);
	void setLambda(LBMReal lambda);
	//void setD(LBMReal D);

	//void initialize(double* f, double x1, double x2, double x3, double vx1, double vx2, double vx3, double rho, UbTupleDouble3 coords, double dx, double o, bool NSE);

	void visit(const SPtr<Grid3D> grid, SPtr<Block3D> block);

protected:
	void checkFunction(mu::Parser fct);
	typedef void(*CalcFeqsFct)(LBMReal* const& /*feq[27]*/, const LBMReal& /*(d)rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);

private:
	mu::Parser muVx1;
	mu::Parser muVx2;
	mu::Parser muVx3;
	//mu::Parser muRho;
	//LBMReal nu;

	//mu::Parser muf1;
	//mu::Parser muf2;
	//mu::Parser muf3;
	mu::Parser muLambda;
	//LBMReal D;
};

#endif //D3Q27INITDISTRIBUTIONSPATCHVISITOR_H

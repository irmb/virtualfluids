#ifndef BasicLBMKernell_h__
#define BasicLBMKernel_h__

#include "LBMKernel.h"


class BasicLBMKernel : public LBMKernel
{
public:
	BasicLBMKernel();
	virtual ~BasicLBMKernel(void);
	virtual void calculate(int step);
	virtual SPtr<LBMKernel> clone() = 0;
protected:
	virtual void initData(){}
	virtual void nodeCollision(int step, int x1, int x2, int x3) {}
	int minX1;
	int minX2;
	int minX3;
	int maxX1;
	int maxX2;
	int maxX3;
};

#endif // BasicLBMKernel_h__
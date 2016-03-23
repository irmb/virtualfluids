/*
 *  D3Q27PressureDifferencePostprocessor.h
 *
 *  Created on: 28.12.2010
 *  Author: kucher
 */

#ifndef D3Q27PRESSUREDIFFERENCEPOSTPROCESSOR_H
#define D3Q27PRESSUREDIFFERENCEPOSTPROCESSOR_H

#include "Postprocessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

class D3Q27PressureDifferencePostprocessor: public Postprocessor {
public:
	D3Q27PressureDifferencePostprocessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& path, 
                                   D3Q27IntegrateValuesHelperPtr h1, D3Q27IntegrateValuesHelperPtr h2, 
                                   LBMReal rhoReal, LBMReal uReal, LBMReal uLB,
                                   /*const LBMUnitConverterPtr conv,*/ CommunicatorPtr comm);
	virtual ~D3Q27PressureDifferencePostprocessor();
	void update(double step);
protected:
	D3Q27IntegrateValuesHelperPtr h1, h2;
   std::string path;
   LBMUnitConverterPtr conv;
	void collectPostprocessData(double step);
   CommunicatorPtr comm;
   LBMReal factor1; //= (1/3)*rhoReal*(uReal/uLB)^2 for calculation pReal = rhoLB * (1/3)*rhoReal*(uReal/uLB)^2, rhoReal and uReal in SI
   LBMReal factor2; //= rhoReal*(uReal/uLB)^2       for calculation pReal = press * rhoReal*(uReal/uLB)^2,       rhoReal and uReal in SI
};


#endif /* D3Q27RHODIFFERENCEPOSTPROCESSOR_H_ */

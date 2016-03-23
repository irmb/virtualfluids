/*
 *  D3Q27PressureDifferencePostprocessor.h
 *
 *  Created on: 28.12.2010
 *  Author: kucher
 */

#ifndef D3Q27PRESSUREDIFFERENCEPOSTPROCESSOR_H
#define D3Q27PRESSUREDIFFERENCEPOSTPROCESSOR_H

#include "CoProcessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

class PressureDifferenceCoProcessor: public CoProcessor {
public:
	PressureDifferenceCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& path, 
                                   D3Q27IntegrateValuesHelperPtr h1, D3Q27IntegrateValuesHelperPtr h2, 
                                   LBMReal rhoReal, LBMReal uReal, LBMReal uLB,
                                   /*const LBMUnitConverterPtr conv,*/ CommunicatorPtr comm);
	virtual ~PressureDifferenceCoProcessor();
	void process(double step);
protected:
	D3Q27IntegrateValuesHelperPtr h1, h2;
   std::string path;
   LBMUnitConverterPtr conv;
	void collectData(double step);
   CommunicatorPtr comm;
   LBMReal factor1; //= (1/3)*rhoReal*(uReal/uLB)^2 for calculation pReal = rhoLB * (1/3)*rhoReal*(uReal/uLB)^2, rhoReal and uReal in SI
   LBMReal factor2; //= rhoReal*(uReal/uLB)^2       for calculation pReal = press * rhoReal*(uReal/uLB)^2,       rhoReal and uReal in SI
};


#endif /* D3Q27RHODIFFERENCEPOSTPROCESSOR_H_ */

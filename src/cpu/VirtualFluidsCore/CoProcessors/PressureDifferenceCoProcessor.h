/*
 *  D3Q27PressureDifferenceCoProcessor.h
 *
 *  Created on: 28.12.2010
 *  Author: kucher
 */

#ifndef D3Q27PRESSUREDIFFERENCECoProcessor_H
#define D3Q27PRESSUREDIFFERENCECoProcessor_H

#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"
#include "LBMSystem.h"

class Communicator;
class Grid3D;
class UbScheduler;
class LBMUnitConverter;
class IntegrateValuesHelper;

class PressureDifferenceCoProcessor: public CoProcessor {
public:
	PressureDifferenceCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path,
        SPtr<IntegrateValuesHelper> h1, SPtr<IntegrateValuesHelper> h2,
                                   LBMReal rhoReal, LBMReal uReal, LBMReal uLB,
                                   /*const SPtr<LBMUnitConverter> conv,*/ SPtr<Communicator> comm);
	~PressureDifferenceCoProcessor() override;

	void process(double step) override;

protected:
    SPtr<IntegrateValuesHelper> h1, h2;
   std::string path;
   SPtr<LBMUnitConverter> conv;
	void collectData(double step);
    SPtr<Communicator> comm;
   LBMReal factor1; //= (1/3)*rhoReal*(uReal/uLB)^2 for calculation pReal = rhoLB * (1/3)*rhoReal*(uReal/uLB)^2, rhoReal and uReal in SI
   LBMReal factor2; //= rhoReal*(uReal/uLB)^2       for calculation pReal = press * rhoReal*(uReal/uLB)^2,       rhoReal and uReal in SI
};


#endif /* D3Q27RHODIFFERENCECoProcessor_H_ */

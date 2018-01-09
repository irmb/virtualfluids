/*
 *  D3Q27PressureDifferenceCoProcessor.h
 *
 *  Created on: 28.12.2010
 *  Author: kucher
 */

#ifndef D3Q27PRESSUREDIFFERENCECoProcessor_H
#define D3Q27PRESSUREDIFFERENCECoProcessor_H

#include <memory>
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
	PressureDifferenceCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, const std::string& path,
        std::shared_ptr<IntegrateValuesHelper> h1, std::shared_ptr<IntegrateValuesHelper> h2,
                                   LBMReal rhoReal, LBMReal uReal, LBMReal uLB,
                                   /*const LBMUnitConverterPtr conv,*/ std::shared_ptr<Communicator> comm);
	virtual ~PressureDifferenceCoProcessor();

	void process(double step) override;

protected:
    std::shared_ptr<IntegrateValuesHelper> h1, h2;
   std::string path;
   std::shared_ptr<LBMUnitConverter> conv;
	void collectData(double step);
    std::shared_ptr<Communicator> comm;
   LBMReal factor1; //= (1/3)*rhoReal*(uReal/uLB)^2 for calculation pReal = rhoLB * (1/3)*rhoReal*(uReal/uLB)^2, rhoReal and uReal in SI
   LBMReal factor2; //= rhoReal*(uReal/uLB)^2       for calculation pReal = press * rhoReal*(uReal/uLB)^2,       rhoReal and uReal in SI
};


#endif /* D3Q27RHODIFFERENCECoProcessor_H_ */

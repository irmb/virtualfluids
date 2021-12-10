/*
 *  Author: K. Kutscher
 *  mail: kutscher@irmb.tu-bs.de
 */
#ifndef RestartDemObjectsCoProcessor_H
#define RestartDemObjectsCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class DemCoProcessor;
class CreateDemObjectsCoProcessor;

class RestartDemObjectsCoProcessor : public CoProcessor
{
public:
    RestartDemObjectsCoProcessor();
    RestartDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                 SPtr<DemCoProcessor> demCoProcessor,
                                 SPtr<CreateDemObjectsCoProcessor> createDemObjectsCoProcessor, double radius,
                                 std::shared_ptr<vf::mpi::Communicator> comm);
    ~RestartDemObjectsCoProcessor() {}
    void process(double step) override;
    void restart(double step);
    void write(int step);
    void read(int step);

private:
    std::string path;
    double radius;
    std::shared_ptr<vf::mpi::Communicator> comm;
    SPtr<DemCoProcessor> demCoProcessor;
    SPtr<CreateDemObjectsCoProcessor> createDemObjectsCoProcessor;
};
#endif

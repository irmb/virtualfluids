#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "mpi.h"
#include "time.h"

#include <Logger/Logger.h>
#include <Input/Input.h>

#include <GridGenerator/grid/GridWrapper/GridWrapperCPU/GridWrapperCPU.h>
#include <GridGenerator/grid/GridWrapper/GridWrapperGPU/GridWrapperGPU.h>
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>
#include <GridGenerator/grid/GridBuilder/GridBuilderImp.h>

#include <GridGenerator/utilities/transformator/Transformator.h>
#include <GridGenerator/utilities/transformator/ArrowTransformator.h>

#include <GridGenerator/io/GridVTKWriter/GridVTKWriter.h>
#include <GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h>

#include <VirtualFluids_GPU/LBM/Simulation.h>
#include <VirtualFluids_GPU/Parameter/Parameter.h>
#include <VirtualFluids_GPU/DataStructureInitializer/GridProvider.h>

#include <utilities/StringUtil.h>

std::string getGridPath(std::shared_ptr<Parameter> para, std::string Gridpath)
{
    if (para->getNumprocs() == 1)
        return Gridpath + "/";
    else
        return Gridpath + "/" + StringUtil::toString(para->getMyID()) + "/";
}

void setParameters(std::shared_ptr<Parameter> para, std::unique_ptr<input::Input> &input)
{
    para->setMaxDev(StringUtil::toInt(input->getValue("NumberOfDevices")));
    para->setDevices(StringUtil::toVector/*<int>*/(input->getValue("Devices")));

    std::string _path = input->getValue("Path");
    std::string _prefix = input->getValue("Prefix");
    std::string _gridpath = input->getValue("GridPath");
    std::string gridPath = getGridPath(para, _gridpath);
    para->setFName(_path + "/" + _prefix);
    para->setPrintFiles(false);
    para->setPrintFiles(StringUtil::toBool(input->getValue("WriteGrid")));
    para->setD3Qxx(StringUtil::toInt(input->getValue("D3Qxx")));
    para->setMaxLevel(StringUtil::toInt(input->getValue("NOGL")));
    para->setTEnd(StringUtil::toInt(input->getValue("TimeEnd")));
    para->setTOut(StringUtil::toInt(input->getValue("TimeOut")));

    para->setViscosity(StringUtil::toFloat(input->getValue("Viscosity_LB")));
    para->setVelocity(StringUtil::toFloat(input->getValue("Velocity_LB")));
    para->setViscosityRatio(StringUtil::toFloat(input->getValue("Viscosity_Ratio_World_to_LB")));
    para->setVelocityRatio(StringUtil::toFloat(input->getValue("Velocity_Ratio_World_to_LB")));
    para->setDensityRatio(StringUtil::toFloat(input->getValue("Density_Ratio_World_to_LB")));
    para->setFactorPressBC(StringUtil::toFloat(input->getValue("dfpbc")));

    para->setgeoVec(gridPath + input->getValue("geoVec"));
    para->setcoordX(gridPath + input->getValue("coordX"));
    para->setcoordY(gridPath + input->getValue("coordY"));
    para->setcoordZ(gridPath + input->getValue("coordZ"));
    para->setneighborX(gridPath + input->getValue("neighborX"));
    para->setneighborY(gridPath + input->getValue("neighborY"));
    para->setneighborZ(gridPath + input->getValue("neighborZ"));
    para->setgeomBoundaryBcQs(gridPath + input->getValue("geomBoundaryBcQs"));
    para->setgeomBoundaryBcValues(gridPath + input->getValue("geomBoundaryBcValues"));
    para->setinletBcQs(gridPath + input->getValue("inletBcQs"));
    para->setinletBcValues(gridPath + input->getValue("inletBcValues"));
    para->setoutletBcQs(gridPath + input->getValue("outletBcQs"));
    para->setoutletBcValues(gridPath + input->getValue("outletBcValues"));
    para->settopBcQs(gridPath + input->getValue("topBcQs"));
    para->settopBcValues(gridPath + input->getValue("topBcValues"));
    para->setbottomBcQs(gridPath + input->getValue("bottomBcQs"));
    para->setbottomBcValues(gridPath + input->getValue("bottomBcValues"));
    para->setfrontBcQs(gridPath + input->getValue("frontBcQs"));
    para->setfrontBcValues(gridPath + input->getValue("frontBcValues"));
    para->setbackBcQs(gridPath + input->getValue("backBcQs"));
    para->setbackBcValues(gridPath + input->getValue("backBcValues"));
    para->setnumberNodes(gridPath + input->getValue("numberNodes"));
    para->setLBMvsSI(gridPath + input->getValue("LBMvsSI"));


    para->setForcing(StringUtil::toFloat(input->getValue("ForcingX")), StringUtil::toFloat(input->getValue("ForcingY")), StringUtil::toFloat(input->getValue("ForcingZ")));
    
    para->setDistX(StringUtil::toVector/*<int>*/(input->getValue("DistX")));
    para->setDistY(StringUtil::toVector/*<int>*/(input->getValue("DistY")));
    para->setDistZ(StringUtil::toVector/*<int>*/(input->getValue("DistZ")));

    //////////////////////////////////////////////////////////////////////////
    //for Multi GPU
    if (para->getNumprocs() > 1)
    {
        //////////////////////////////////////////////////////////////////////////
        //3D domain decomposition
        std::vector<std::string> sendProcNeighborsX, sendProcNeighborsY, sendProcNeighborsZ;
        std::vector<std::string> recvProcNeighborsX, recvProcNeighborsY, recvProcNeighborsZ;
        for (int i = 0; i < para->getNumprocs(); i++)
        {
            sendProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xs.dat");
            sendProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Ys.dat");
            sendProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zs.dat");
            recvProcNeighborsX.push_back(gridPath + StringUtil::toString(i) + "Xr.dat");
            recvProcNeighborsY.push_back(gridPath + StringUtil::toString(i) + "Yr.dat");
            recvProcNeighborsZ.push_back(gridPath + StringUtil::toString(i) + "Zr.dat");
        }
        para->setPossNeighborFilesX(sendProcNeighborsX, "send");
        para->setPossNeighborFilesY(sendProcNeighborsY, "send");
        para->setPossNeighborFilesZ(sendProcNeighborsZ, "send");
        para->setPossNeighborFilesX(recvProcNeighborsX, "recv");
        para->setPossNeighborFilesY(recvProcNeighborsY, "recv");
        para->setPossNeighborFilesZ(recvProcNeighborsZ, "recv");
    }
}


static void buildGrid(std::shared_ptr<GridBuilder> builder, double length, double width, double high, double delta, std::string distribution, std::shared_ptr<Transformator> trans, std::string stlGeometryPath)
{
    MPI_Barrier(MPI_COMM_WORLD);
    clock_t begin = clock();

    
    builder->addGrid((doubflo)length, (doubflo)width, (doubflo)high, (doubflo)delta, distribution, trans);

    if(stlGeometryPath != "")
        builder->meshGeometry(stlGeometryPath, 0);

    MPI_Barrier(MPI_COMM_WORLD);
    clock_t end = clock();

    double time = double(end - begin) / CLOCKS_PER_SEC;
    *logging::out << logging::Logger::INTERMEDIATE << "Complete time calculation grid: " << std::to_string(time) << "s\n";
}

static void generateGrid(const std::string &configFile)
{
    *logging::out << logging::Logger::LOW << configFile << "\n";

    std::ifstream stream;
    stream.open(configFile.c_str(), std::ios::in);
    if (stream.fail())
        throw "can not open config file!\n";

    std::unique_ptr<input::Input> input = input::Input::makeInput(stream, "config");

    double length = std::stod(input->getValue("channelLength"));
    double width = std::stod(input->getValue("channelWidth"));
    double high = std::stod(input->getValue("channelHeight"));
    double delta = std::stod(input->getValue("channelDelta"));

    double transX = std::stod(input->getValue("translateX"));
    double transY = std::stod(input->getValue("translateY"));
    double transZ = std::stod(input->getValue("translateZ"));

    bool computeOnGPU    = StringUtil::toBool(input->getValue("computeOnGPU"));

    std::string distribution = input->getValue("D3Qxx");
    std::string inputGeometry = input->getValue("STL_INPUT");
    std::string arrowOutput = input->getValue("ARROW_OUTPUT");

    std::shared_ptr<GridBuilder> builder;
    if (computeOnGPU)
        builder = GridBuilderImp::make("gpu");
    else
        builder = GridBuilderImp::make("cpu");

    std::shared_ptr<Transformator> trans = Transformator::makeTransformator((doubflo)delta, (doubflo)transX, (doubflo)transY, (doubflo)transZ);

    buildGrid(builder, length, width, high, delta, distribution, trans, inputGeometry);

    builder->getKernel(0, 0)->copyDataFromGPU();
    if (StringUtil::toBool(input->getValue("write_VTK")))
        builder->writeGridToVTK(input->getValue("VTK_OUTPUT"), 0);

    if (StringUtil::toBool(input->getValue("startSimulation")))
    {
        builder->deleteSolidNodes();

        GridVTKWriter::writeSparseGridToVTK(builder->getKernel(0, 0)->grid, input->getValue("VTK_OUTPUT") ,trans);
        SimulationFileWriter::writeSimulationFiles("E:\\GRIDGENERATION\\files\\", builder, false, trans);

        auto arrowTransformator = ArrowTransformator::makeTransformator((doubflo)delta, (doubflo)transX, (doubflo)transY, (doubflo)transZ);
        if (StringUtil::toBool(input->getValue("write_arrows")))
            builder->writeArrows(arrowOutput, arrowTransformator);

        std::shared_ptr<Parameter> para = Parameter::make();
        std::shared_ptr<GridProvider> gridProvider = GridProvider::makeGridGenerator(builder, para);
        //std::shared_ptr<GridProvider> reader = GridProvider::makeGridReader(true, para);

        setParameters(para, input);
        Simulation sim;
        sim.init(para, gridProvider);
        sim.run();
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::ofstream os("log.txt");

    if (rank == 0)
        logging::Logger::setStream(std::cout);
    else
        logging::Logger::setStream(os);

    logging::Logger::setDebugLevel(logging::Logger::INTERMEDIATE);
    logging::Logger::enablePrintedRankNumbers(true);

    MPI_Barrier(MPI_COMM_WORLD);
    

    if (argc > 1)
        generateGrid(argv[1]);

    MPI_Finalize();
	return 0;
}

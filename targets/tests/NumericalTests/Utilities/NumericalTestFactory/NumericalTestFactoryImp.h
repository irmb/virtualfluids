#ifndef NUMERICAL_TEST_FACTORY_IMP_H
#define NUMERICAL_TEST_FACTORY_IMP_H

#include "NumericalTestFactory.h"

struct ConfigDataStruct;
struct GridInformationStruct;
struct L2NormTestStruct;
struct L2NormTestParameterStruct;
struct L2NormTestBetweenKernelsParameterStruct;
struct L2NormTestBetweenKernelsStruct;
struct LogFileParameterStruct;
struct NumericalTestStruct;
struct PhiAndNyTestStruct;
struct PhiAndNyTestParameterStruct;
struct ShearWaveParameterStruct;
struct SimulationDataStruct;
struct TaylorGreenVortexUxParameterStruct;
struct TaylorGreenVortexUzParameterStruct;
struct TestSimulationDataStruct;
struct VectorWriterInformationStruct;

class AnalyticalResults2DToVTKWriter;
class BasicTestLogFileInformation;
class ColorConsoleOutput;
class L2NormTest;
class L2NormPostProcessingStrategy;
class L2NormBetweenKernelPostProcessingStrategy;
class L2NormTestBetweenKernels;
class LogFileTimeInformation;
class LogFileQueueImp;
class LogFileWriter;
class PhiAndNyTest;
class PhiAndNyTestPostProcessingStrategy;
class SimulationInfo;
class SimulationLogFileInformation;
class SimulationParameter;
class SimulationResults;
class TestQueueImp;
class TestSimulation;
class TestSimulationImp;
class TestLogFileInformation;

class NumericalTestFactoryImp : public NumericalTestFactory
{
public:
	static std::shared_ptr<NumericalTestFactoryImp> getNewInstance(std::shared_ptr<ConfigDataStruct> configFileData);

	std::vector<std::shared_ptr<TestSimulation> > getTestSimulations();
	std::shared_ptr<TestQueue> getTestQueue();
	std::shared_ptr<LogFileQueue> getLogFileQueue();

private:
	NumericalTestFactoryImp() {};
	NumericalTestFactoryImp(std::shared_ptr<ConfigDataStruct> configFileData);

	void init(std::shared_ptr<ConfigDataStruct> configFileData);

	std::shared_ptr<NumericalTestStruct> makeNumericalTestStruct(std::shared_ptr<ConfigDataStruct> configFileData, std::shared_ptr<SimulationDataStruct> simDataStruct, std::string kernel, double viscosity, int basicTimeStepLength);
	void addNumericalTestStruct(std::shared_ptr<NumericalTestStruct> numericalTestStruct);

	std::shared_ptr<SimulationDataStruct> makeTaylorGreenUxSimulationData(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);
	std::shared_ptr<SimulationDataStruct> makeTaylorGreenUzSimulationData(std::string kernelName, double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);
	std::shared_ptr<SimulationDataStruct> makeShearWaveSimulationData(std::string kernelName, double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);

	std::vector<std::shared_ptr<TestSimulationImp> > makeTestSimulations(std::vector<std::shared_ptr<TestSimulationDataStruct> > testSimDataStruct, std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo, unsigned int ySliceForCalculation);

	std::shared_ptr<PhiAndNyTestStruct> makePhiAndNyTestsStructs(std::shared_ptr<PhiAndNyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations, double viscosity);
	std::shared_ptr<L2NormTestStruct> makeL2NormTestsStructs(std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSimumlations);
	std::shared_ptr<L2NormTestBetweenKernelsStruct> makeL2NormTestsBetweenKernelsStructs(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::string kernelName);
	std::vector<std::shared_ptr<PhiAndNyTest> > makePhiAndNyTests(std::shared_ptr<PhiAndNyTestParameterStruct> testParameter, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<PhiAndNyTestPostProcessingStrategy> > phiAndNuPostProStrategy, double viscosity, std::string dataToCalculate);
	std::vector<std::shared_ptr<L2NormTest> > makeL2NormTests(std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormPostProcessingStrategy> > postProStrategy, std::shared_ptr<L2NormTestParameterStruct> testParameter);
	std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > makeL2NormTestsBetweenKernels(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies);
	
	std::vector<std::shared_ptr<L2NormTestBetweenKernels> > linkL2NormTestsBetweenKernels(std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> testPara, std::vector<std::shared_ptr<TestSimulationImp> > testSim, std::vector<std::shared_ptr<L2NormBetweenKernelPostProcessingStrategy> > postProcessingStrategies);

	std::shared_ptr<LogFileWriter> makeLogFileWriter(std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::vector<std::shared_ptr<SimulationInfo> > simInfo, std::string kernelName, double viscosity, int basicTimeStepLength, std::shared_ptr<LogFileParameterStruct> logFilePara);

	std::vector<std::shared_ptr<TestSimulation> > myTestSimulations;
	std::shared_ptr<TestQueueImp> myTestQueue;
	std::shared_ptr<LogFileQueueImp> myLogFileWriterQueue;

	std::vector<std::vector<std::shared_ptr<L2NormTestBetweenKernels> > > l2NormTestsBetweenKernels;


	std::shared_ptr<ColorConsoleOutput> colorOutput;
	std::shared_ptr<AnalyticalResults2DToVTKWriter> anaResultWriter;
	std::shared_ptr<BasicTestLogFileInformation> basicTestLogFileInfo;

	int simID;
	int numberOfSimulations;
	int simPerKernel, numberOfTestGroupsBetweenKernels, numberOfTestsForOneSimulation, numberOfTestsBetweenKernels;
	int posBasicSimulationForL2Test, posDivergentSimulationForL2Test;
};
#endif
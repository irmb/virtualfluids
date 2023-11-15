#ifndef PHI_TEST_H
#define PHI_TEST_H

#include "Utilities/Test/TestImp.h"

#include <memory>
#include <vector>
#include <iostream>

class FFTCalculator;
class PhiTestPostProcessingStrategy;
struct PhiTestParameterStruct;

class PhiTest : public TestImp 
{
public:
    static std::shared_ptr<PhiTest> getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiTestParameterStruct> testPara, std::string dataToCalculate);
    
    void update();
    void addSimulation(std::shared_ptr<NumericalTestSimulation> sim, std::shared_ptr<SimulationInfo> simInfo, std::shared_ptr<PhiTestPostProcessingStrategy> postProStrategy);
    void evaluate();

    std::string getDataToCalculate();
    std::vector<int> getLx();
    std::vector<double> getPhiDiff();
    double getOrderOfAccuracy();



private:
    PhiTest(std::shared_ptr<ColorConsoleOutput> colorOutput, double viscosity, std::shared_ptr<PhiTestParameterStruct> testPara, std::string dataToCalculate);
    double calcOrderOfAccuracy(std::vector<double> data);
    TestStatus checkTestPassed(double orderOfAccuracy);
    std::vector<std::string> buildTestOutput();
    std::vector<std::string> buildBasicTestOutput();
    std::vector<std::string> buildErrorTestOutput();

    unsigned int startStepCalculation, endStepCalculation;
    std::vector<double> lx;
    std::vector<double> phiDiff;
    double orderOfAccuracy;
    double minOrderOfAccuracy;
    std::string dataToCalculate;

    std::vector<std::shared_ptr<PhiTestPostProcessingStrategy> > postProStrategies;

};
#endif

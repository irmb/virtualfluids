#ifndef NUMERICAL_TEST_STRUCT_H
#define NUMERICAL_TEST_STRUCT_H

#include <memory>
#include <vector>

class LogFileWriter;
class Test;
class TestSimulationImp;

struct NumericalTestStruct
{
	std::vector<std::shared_ptr<TestSimulationImp> > testSimulations;
	std::vector<std::shared_ptr<Test> > tests;
	std::shared_ptr<LogFileWriter> logFileWriter;
};
#endif
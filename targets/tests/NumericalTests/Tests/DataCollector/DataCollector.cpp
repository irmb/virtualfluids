#include "DataCollector.h"
#include "Tests\DataQueue\DataQueue.h"
#include "Utilities\EvaluationParameter\EvaluationParameter.h"


DataCollector::DataCollector(DataQueue * dataQueueNu, DataQueue * dataQueuePhi, const int arraySize, std::string testName)
{
	this->dataQueueNu = dataQueueNu;
	this->dataQueuePhi = dataQueuePhi;
	this->arraySize = arraySize;
	this->testName = testName;
	for (int i = 0; i < arraySize; i++) {
		dataQueueNu[i].valueName = "NuDiff";
		dataQueuePhi[i].valueName = "Phi";
	}
}

void DataCollector::addNuDiffAndPhi(double nudiff, double phi, std::shared_ptr<EvaluationParameter> evaPara)
{
	if (evaPara->getTestName() == testName) {
		addDataToQueue(nudiff, dataQueueNu, evaPara);
		addDataToQueue(phi, dataQueuePhi, evaPara);
	}
}

void DataCollector::addDataToQueue(double data, DataQueue * dataQueue, std::shared_ptr<EvaluationParameter> evaPara)
{
	int n = dataQueue[0].numberOfTests;
	if (n == 0) {
		dataQueue[n].a = data;
		dataQueue[n].la = evaPara->getLx();
	}
	else if (n == evaPara->getNumberOfEqualTests() - 1) {
		dataQueue[n - 1].b = data;
		dataQueue[n - 1].lb = evaPara->getLx();
		dataQueue[n - 1].orderOfAccuracy = log(dataQueue[n - 1].a / dataQueue[n - 1].b) / log(dataQueue[n - 1].lb / dataQueue[n - 1].la);
		dataQueue[n - 1].minOrderOfAccuracy = evaPara->getMinOrderOfAccuracy();
		dataQueue[n - 1].testName = evaPara->getTestName();
		dataQueue[n - 1].expected = true;
	}
	else {
		dataQueue[n - 1].b = data;
		dataQueue[n - 1].lb = evaPara->getLx();
		dataQueue[n - 1].orderOfAccuracy = log(dataQueue[n - 1].a / dataQueue[n - 1].b) / log(dataQueue[n - 1].lb / dataQueue[n - 1].la);
		dataQueue[n - 1].minOrderOfAccuracy = evaPara->getMinOrderOfAccuracy();
		dataQueue[n - 1].testName = evaPara->getTestName();
		dataQueue[n - 1].expected = true;
		dataQueue[n].a = data;
		dataQueue[n].la = evaPara->getLx();
	}
	increaseNumberOfTestByOne(dataQueue);
}

void DataCollector::increaseNumberOfTestByOne(DataQueue * dataQueue)
{
	for (int i = 0; i < arraySize; i++)
		dataQueue[i].numberOfTests++;
}



#ifndef DATACONTAINERINITIALIZER_H
#define DATACONTAINERINITIALIZER_H

#include <memory>
#include <string>

struct DataQueue;
class EvaluationParameter;

class DataCollector
{
public:
	DataCollector(DataQueue *DataQueueNu, DataQueue *DataQueuePhi, const int arraySize, std::string testName);
	void addNuDiffAndPhi(double nudiff, double phi, std::shared_ptr<EvaluationParameter> evaPara);
	

private:
	void addDataToQueue(double data, DataQueue *dataQueue, std::shared_ptr<EvaluationParameter> evaPara);
	void increaseNumberOfTestByOne(DataQueue *dataQueue);

	int arraySize;
	DataQueue *dataQueuePhi;
	DataQueue *dataQueueNu;
	std::string testName;
};
#endif
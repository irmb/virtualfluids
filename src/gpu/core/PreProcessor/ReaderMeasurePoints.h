#ifndef POSITION_READER_H
#define POSITION_READER_H

class Parameter;
class CudaMemoryManager;

class ReaderMeasurePoints
{
private:
    static void readMeasurePoints(Parameter* para);

public:
    static void readMeasurePoints(Parameter* para, CudaMemoryManager* cudaMemoryManager);
};

#endif
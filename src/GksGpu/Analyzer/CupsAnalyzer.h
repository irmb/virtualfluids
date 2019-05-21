#ifndef  CupsAnalyzer_H
#define  CupsAnalyzer_H

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/Timer/Timer.h"

struct DataBase;

class VF_PUBLIC CupsAnalyzer
{
private:
    Timer timer;
    Timer timerRestart;

    bool outputPerTime;

    bool outputPerIter;

    real outputTime;
    uint outputPerTimeCounter;

    uint outputIter;

    unsigned long long numberOfCellUpdatesPerTimeStep;

    uint counter;

public:

    CupsAnalyzer( SPtr<DataBase> dataBase, 
                  bool outputPerTime = true, real outputTime = 600.0,
                  bool outputPerIter = true, uint outputIter = 10000 );

    void start();

    void restart();

    void run( uint iter );

private:

    bool checkOutputPerTime( real currentRuntime );
    bool checkOutputPerIter( uint iter );

    void printCups(uint iter, real currentRunTime, real cups);
};

#endif

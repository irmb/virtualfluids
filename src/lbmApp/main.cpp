#include <iostream>
#include <stdio.h>

#include "lbm/CalcMac.h"
//#include <cuda_runtime.h>

/* __global__ */ void test()
{
    printf("Hello World from GPU!\n");
    real f[27] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    double density = LBM::getDensity(f);

    printf("Hello density: %f \n", density);
}

int main()
{
    std::cout << "hello world \n";
    test();
    //test<<<1,1>>>();
    //cudaDeviceSynchronize();

    return 0;
}
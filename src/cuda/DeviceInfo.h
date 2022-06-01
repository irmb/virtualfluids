#ifndef CUDA_DEVICEINFO_H
#define CUDA_DEVICEINFO_H

namespace vf::cuda
{

void verifyAndSetDevice(int deviceId);

void printCudaInformation(int deviceId);

}

#endif

#ifndef AccumulatorDataType_H
#define AccumulatorDataType_H

namespace GksGpu {

// This file is used to control the data type of accumulator variables.
// Accumulator variables are variables, where cell values are written 
// during the flux computation, which is per face. Since the face evaluation 
// order on GPUs is arbitrary, the cutoff errors for these accumulators are non 
// deterministic. This deficiency can be solved for single precision calculations
// by setting the accumulator data type to double. The deviations are then 
// so small, that they are cut off during the downcast to single.
// using double precision accumulators has some performance implications, 
// especially on consumer hardware.

//typedef float realAccumulator;
typedef double realAccumulator;

} // namespace GksGpu

#endif

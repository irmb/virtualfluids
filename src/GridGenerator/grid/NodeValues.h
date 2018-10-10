#ifndef NodeValues_H
#define NodeValues_H

#define FLUID 0

#define FLUID_CFC 1
#define FLUID_CFF 2

#define FLUID_FCC 3
#define FLUID_FCF 4

#define MULTI_GPU_SEND 10
#define MULTI_GPU_RECIEVE 11

#define BC_PRESSURE 20
#define BC_VELOCITY 21
#define BC_SOLID 22

#define BC_SLIP 23
#define BC_NOSLIP 24
#define BC_OUTFLOW 25

#define STOPPER_OUT_OF_GRID 30
#define STOPPER_COARSE_UNDER_FINE 31
#define STOPPER_SOLID 32
#define STOPPER_OUT_OF_GRID_BOUNDARY 33

#define INVALID_OUT_OF_GRID 40
#define INVALID_COARSE_UNDER_FINE 41
#define INVALID_SOLID 42

//????WTF?????
#define INSIDE 50
#define NEGATIVE_DIRECTION_BORDER 51
#define Q_DEPRECATED 52

#define OVERLAP_TMP 60

#endif

#ifndef NodeValues_H
#define NodeValues_H

namespace vf
{
namespace gpu
{

static constexpr char FLUID = 0;

static constexpr char FLUID_CFC = 1;
static constexpr char FLUID_CFF = 2;

static constexpr char FLUID_FCC = 3;
static constexpr char FLUID_FCF = 4;

static constexpr char MULTI_GPU_SEND = 10;
static constexpr char MULTI_GPU_RECIEVE = 11;

static constexpr char BC_PRESSURE = 20;
static constexpr char BC_VELOCITY = 21;
static constexpr char BC_SOLID = 22;

static constexpr char BC_SLIP = 23;
static constexpr char BC_NOSLIP = 24;
static constexpr char BC_OUTFLOW = 25;

static constexpr char STOPPER_OUT_OF_GRID = 30;
static constexpr char STOPPER_COARSE_UNDER_FINE = 31;
static constexpr char STOPPER_SOLID = 32;
static constexpr char STOPPER_OUT_OF_GRID_BOUNDARY = 33;

static constexpr char INVALID_OUT_OF_GRID = 40;
static constexpr char INVALID_COARSE_UNDER_FINE = 41;
static constexpr char INVALID_SOLID = 42;

//????WTF?????
static constexpr char INSIDE = 50;
static constexpr char NEGATIVE_DIRECTION_BORDER = 51;
static constexpr char Q_DEPRECATED = 52;

static constexpr char OVERLAP_TMP = 60;

}
}

#endif

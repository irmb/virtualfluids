#ifndef REAL_CONSTANT_H
#define REAL_CONSTANT_H


namespace vf
{
namespace lbm 
{
namespace constant
{

#ifdef VF_DOUBLE_ACCURACY
static constexpr double c1o2 = 0.5;
static constexpr double c3o2 = 1.5;
static constexpr double c1o3 = 0.333333333333333;
static constexpr double c2o3 = 0.666666666666667;
static constexpr double c1o4 = 0.25;
static constexpr double c3o4 = 0.75;
static constexpr double c1o6 = 0.166666666666667;
static constexpr double c1o7 = 0.142857142857143;
static constexpr double c1o8 = 0.125;
static constexpr double c1o9 = 0.111111111111111;
static constexpr double c2o9 = 0.222222222222222;
static constexpr double c4o9 = 0.444444444444444;
static constexpr double c1o10 = 0.1;
static constexpr double c1o12 = 0.083333333333333;
static constexpr double c1o16 = 0.0625;
static constexpr double c3o16 = 0.1875;
static constexpr double c9o16 = 0.5625;
static constexpr double c1o18 = 0.055555555555556;
static constexpr double c1o20 = 0.05;
static constexpr double c19o20 = 0.95;
static constexpr double c21o20 = 1.05;
static constexpr double c1o24 = 0.041666666666667;
static constexpr double c1o27 = 0.037037037037037;
static constexpr double c3o32 = 0.09375;
static constexpr double c4o32 = 0.125;
static constexpr double c1o36 = 0.027777777777778;
static constexpr double c1o48 = 0.020833333333333;
static constexpr double c1o64 = 0.015625;
static constexpr double c3o64 = 0.046875;
static constexpr double c9o64 = 0.140625;
static constexpr double c27o64 = 0.421875;
static constexpr double c1o66 = 0.015151515151515;
static constexpr double c1o72 = 0.013888888888889;
static constexpr double c1o264 = 0.003787878787879;
static constexpr double c8o27 = 0.296296296296296;
static constexpr double c2o27 = 0.074074074074074;
static constexpr double c1o54 = 0.018518518518519;
static constexpr double c1o100 = 0.01;
static constexpr double c99o100 = 0.99;
static constexpr double c1o126 = 0.007936507936508;
static constexpr double c1o216 = 0.004629629629630;
static constexpr double c5o4 = 1.25;
static constexpr double c9o4 = 2.25;
static constexpr double c5o2 = 2.5;
static constexpr double c9o2 = 4.5;

static constexpr double c0o1 = 0.;
static constexpr double c1o1 = 1.;
static constexpr double c2o1 = 2.;
static constexpr double c3o1 = 3.;
static constexpr double c4o1 = 4.;
static constexpr double c5o1 = 5.;
static constexpr double c6o1 = 6.;
static constexpr double c7o1 = 7.;
static constexpr double c8o1 = 8.;
static constexpr double c9o1 = 9.;
static constexpr double c10o1 = 10.;
static constexpr double c11o1 = 11.;
static constexpr double c12o1 = 12.;
static constexpr double c13o1 = 13.;
static constexpr double c14o1 = 14.;
static constexpr double c15o1 = 15.;
static constexpr double c16o1 = 16.;
static constexpr double c17o1 = 17.;
static constexpr double c18o1 = 18.;
static constexpr double c21o1 = 21.;
static constexpr double c24o1 = 24.;
static constexpr double c25o1 = 25.;
static constexpr double c26o1 = 26.;
static constexpr double c27o1 = 27.;
static constexpr double c28o1 = 28.;
static constexpr double c29o1 = 29.;
static constexpr double c30o1 = 30.;
static constexpr double c32o1 = 32.;
static constexpr double c33o1 = 33.;
static constexpr double c34o1 = 34.;
static constexpr double c36o1 = 36.;
static constexpr double c40o1 = 40.;
static constexpr double c42o1 = 42.;
static constexpr double c46o1 = 46.;
static constexpr double c48o1 = 48.;
static constexpr double c50o1 = 50.;
static constexpr double c52o1 = 52.;
static constexpr double c54o1 = 54.;
static constexpr double c56o1 = 56.;
static constexpr double c64o1 = 64.;
static constexpr double c66o1 = 66.;
static constexpr double c68o1 = 68.;
static constexpr double c69o1 = 69.;
static constexpr double c72o1 = 72.;
static constexpr double c84o1 = 84.;
static constexpr double c88o1 = 88.;
static constexpr double c96o1 = 96.;
static constexpr double c100o1 = 10.;
static constexpr double c130o1 = 13.;
static constexpr double c152o1 = 15.;
static constexpr double c166o1 = 16.;
static constexpr double c195o1 = 19.;
static constexpr double c216o1 = 21.;
static constexpr double c264o1 = 26.;
static constexpr double c290o1 = 29.;
static constexpr double c367o1 = 36.;

static constexpr double Op0000002 = 0.0000002;
static constexpr double c10eM30 = 1e-30;
static constexpr double c10eM10 = 1e-10;
static constexpr double smallSingle = 0.0000000002;

static constexpr double cPi = 3.1415926535;
static constexpr double cPio180 = 1.74532925199e-2;
static constexpr double c180oPi = 57.2957795131;

#else
static constexpr float c1o2 = 0.5f;
static constexpr float c3o2 = 1.5f;
static constexpr float c1o3 = (1.0f / 3.0f);
static constexpr float c2o3 = (2.0f / 3.0f);
static constexpr float c1o4 = 0.25f;
static constexpr float c3o4 = 0.75f;
static constexpr float c1o6 = (1.0f / 6.0f);
static constexpr float c1o7 = (1.0f / 7.0f);
static constexpr float c1o8 = 0.125f;
static constexpr float c1o9 = (1.0f / 9.0f);
static constexpr float c2o9 = (2.0f / 9.0f);
static constexpr float c4o9 = (4.0f / 9.0f);
static constexpr float c1o10 = 0.1f;
static constexpr float c1o12 = (1.0f / 12.0f);
static constexpr float c1o16 = 0.0625f;
static constexpr float c3o16 = 0.1875f;
static constexpr float c9o16 = 0.5625f;
static constexpr float c1o18 = (1.0f / 18.0f);
static constexpr float c1o20 = 0.05f;
static constexpr float c19o20 = 0.95f;
static constexpr float c21o20 = 1.05f;
static constexpr float c1o24 = (1.0f / 24.0f);
static constexpr float c1o27 = (1.0f / 27.0f);
static constexpr float c3o32 = 0.09375f;
static constexpr float c4o32 = 0.125f;
static constexpr float c1o36 = (1.0f / 36.0f);
static constexpr float c1o48 = (1.0f / 48.0f);
static constexpr float c1o64 = 0.015625f;
static constexpr float c3o64 = 0.046875f;
static constexpr float c9o64 = 0.140625f;
static constexpr float c27o64 = 0.421875f;
static constexpr float c1o66 = (1.0f / 66.0f);
static constexpr float c1o72 = (1.0f / 72.0f);
static constexpr float c1o264 = (1.0f / 264.0f);
static constexpr float c8o27 = (8.0f / 27.0f);
static constexpr float c2o27 = (2.0f / 27.0f);
static constexpr float c1o54 = (1.0f / 54.0f);
static constexpr float c1o100 = 0.01f;
static constexpr float c99o100 = 0.99f;
static constexpr float c1o126 = (1.0f / 126.0f);
static constexpr float c1o216 = (1.0f / 216.0f);
static constexpr float c5o4 = 1.25f;
static constexpr float c9o4 = 2.25f;
static constexpr float c5o2 = 2.5f;
static constexpr float c9o2 = 4.5f;

static constexpr float c0o1 = 0.f;
static constexpr float c1o1 = 1.f;
static constexpr float c2o1 = 2.f;
static constexpr float c3o1 = 3.f;
static constexpr float c4o1 = 4.f;
static constexpr float c5o1 = 5.f;
static constexpr float c6o1 = 6.f;
static constexpr float c7o1 = 7.f;
static constexpr float c8o1 = 8.f;
static constexpr float c9o1 = 9.f;
static constexpr float c10o1 = 10.f;
static constexpr float c11o1 = 11.f;
static constexpr float c12o1 = 12.f;
static constexpr float c13o1 = 13.f;
static constexpr float c14o1 = 14.f;
static constexpr float c15o1 = 15.f;
static constexpr float c16o1 = 16.f;
static constexpr float c17o1 = 17.f;
static constexpr float c18o1 = 18.f;
static constexpr float c21o1 = 21.f;
static constexpr float c24o1 = 24.f;
static constexpr float c25o1 = 25.f;
static constexpr float c26o1 = 26.f;
static constexpr float c27o1 = 27.f;
static constexpr float c28o1 = 28.f;
static constexpr float c29o1 = 29.f;
static constexpr float c30o1 = 30.f;
static constexpr float c32o1 = 32.f;
static constexpr float c33o1 = 33.f;
static constexpr float c34o1 = 34.f;
static constexpr float c36o1 = 36.f;
static constexpr float c40o1 = 40.f;
static constexpr float c42o1 = 42.f;
static constexpr float c46o1 = 46.f;
static constexpr float c48o1 = 48.f;
static constexpr float c50o1 = 50.f;
static constexpr float c52o1 = 52.f;
static constexpr float c54o1 = 54.f;
static constexpr float c56o1 = 56.f;
static constexpr float c64o1 = 64.f;
static constexpr float c66o1 = 66.f;
static constexpr float c68o1 = 68.f;
static constexpr float c69o1 = 69.f;
static constexpr float c72o1 = 72.f;
static constexpr float c84o1 = 84.f;
static constexpr float c88o1 = 88.f;
static constexpr float c96o1 = 96.f;
static constexpr float c100o1 = 100.0f;
static constexpr float c130o1 = 130.0f;
static constexpr float c152o1 = 152.0f;
static constexpr float c166o1 = 166.0f;
static constexpr float c195o1 = 195.0f;
static constexpr float c216o1 = 216.0f;
static constexpr float c264o1 = 264.0f;
static constexpr float c290o1 = 290.0f;
static constexpr float c367o1 = 367.0f;

static constexpr float Op0000002 = 0.0000002f;
static constexpr float c10eM30 = 1e-30f;
static constexpr float c10eM10 = 1e-10f;
static constexpr float smallSingle = 0.0000000002f;

static constexpr float cPi = 3.1415926535f;
static constexpr float cPio180 = 0.0174532925199f;
static constexpr float c180oPi = 57.2957795131f;

#endif

}
}
}

#endif

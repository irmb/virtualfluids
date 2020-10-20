#ifndef REAL_CONSTANT_H
#define REAL_CONSTANT_H

#ifdef VF_DOUBLE_ACCURACY
#define c1o2 0.5
#define c3o2 1.5
#define c1o3 0.333333333333333
#define c2o3 0.666666666666667
#define c1o4 0.25
#define c3o4 0.75
#define c1o6 0.166666666666667
#define c1o7 0.142857142857143
#define c1o8 0.125
#define c1o9 0.111111111111111
#define c2o9 0.222222222222222
#define c4o9 0.444444444444444
#define c1o10 0.1
#define c1o12 0.083333333333333
#define c1o16 0.0625
#define c3o16 0.1875
#define c9o16 0.5625
#define c1o18 0.055555555555556
#define c1o20 0.05
#define c19o20 0.95
#define c21o20 1.05
#define c1o24 0.041666666666667
#define c1o27 0.037037037037037
#define c3o32 0.09375
#define c4o32 0.125
#define c1o36 0.027777777777778
#define c1o48 0.020833333333333
#define c1o64 0.015625
#define c3o64 0.046875
#define c9o64 0.140625
#define c27o64 0.421875
#define c1o66 0.015151515151515
#define c1o72 0.013888888888889
#define c1o264 0.003787878787879
#define c8o27 0.296296296296296
#define c2o27 0.074074074074074
#define c1o54 0.018518518518519
#define c1o100 0.01
#define c99o100 0.99
#define c1o126 0.007936507936508
#define c1o216 0.004629629629630
#define c5o4 1.25
#define c9o4 2.25
#define c5o2 2.5
#define c9o2 4.5

#define c0o1 0.
#define c1o1 1.
#define c2o1 2.
#define c3o1 3.
#define c4o1 4.
#define c5o1 5.
#define c6o1 6.
#define c7o1 7.
#define c8o1 8.
#define c9o1 9.
#define c10o1 10.
#define c11o1 11.
#define c12o1 12.
#define c13o1 13.
#define c14o1 14.
#define c15o1 15.
#define c16o1 16.
#define c17o1 17.
#define c18o1 18.
#define c21o1 21.
#define c24o1 24.
#define c25o1 25.
#define c26o1 26.
#define c27o1 27.
#define c28o1 28.
#define c29o1 29.
#define c30o1 30.
#define c32o1 32.
#define c33o1 33.
#define c34o1 34.
#define c36o1 36.
#define c40o1 40.
#define c42o1 42.
#define c46o1 46.
#define c48o1 48.
#define c50o1 50.
#define c52o1 52.
#define c54o1 54.
#define c56o1 56.
#define c64o1 64.
#define c66o1 66.
#define c68o1 68.
#define c69o1 69.
#define c72o1 72.
#define c84o1 84.
#define c88o1 88.
#define c96o1 96.
#define c100o1 100.0
#define c130o1 130.0
#define c152o1 152.0
#define c166o1 166.0
#define c195o1 195.0
#define c216o1 216.0
#define c264o1 264.0
#define c290o1 290.0
#define c367o1 367.0

#define Op0000002 0.0000002
#define c10eM30 1e-30
#define c10eM10 1e-10
#define smallSingle 0.0000000002

#else
#define c1o2 0.5f
#define c3o2 1.5f
#define c1o3 (1.0f / 3.0f)
#define c2o3 (2.0f / 3.0f)
#define c1o4 0.25f
#define c3o4 0.75f
#define c1o6 (1.0f / 6.0f)
#define c1o7 (1.0f / 7.0f)
#define c1o8 0.125f
#define c1o9 (1.0f / 9.0f)
#define c2o9 (2.0f / 9.0f)
#define c4o9 (4.0f / 9.0f)
#define c1o10 0.1f
#define c1o12 (1.0f / 12.0f)
#define c1o16 0.0625f
#define c3o16 0.1875f
#define c9o16 0.5625f
#define c1o18 (1.0f / 18.0f)
#define c1o20 0.05f
#define c19o20 0.95f
#define c21o20 1.05f
#define c1o24 (1.0f / 24.0f)
#define c1o27 (1.0f / 27.0f)
#define c3o32 0.09375f
#define c4o32 0.125f
#define c1o36 (1.0f / 36.0f)
#define c1o48 (1.0f / 48.0f)
#define c1o64 0.015625f
#define c3o64 0.046875f
#define c9o64 0.140625f
#define c27o64 0.421875f
#define c1o66 (1.0f / 66.0f)
#define c1o72 (1.0f / 72.0f)
#define c1o264 (1.0f / 264.0f)
#define c8o27 (8.0f / 27.0f)
#define c2o27 (2.0f / 27.0f)
#define c1o54 (1.0f / 54.0f)
#define c1o100 0.01f
#define c99o100 0.99f
#define c1o126 (1.0f / 126.0f)
#define c1o216 (1.0f / 216.0f)
#define c5o4 1.25f
#define c9o4 2.25f
#define c5o2 2.5f
#define c9o2 4.5f

#define c0o1 0.f
#define c1o1 1.f
#define c2o1 2.f
#define c3o1 3.f
#define c4o1 4.f
#define c5o1 5.f
#define c6o1 6.f
#define c7o1 7.f
#define c8o1 8.f
#define c9o1 9.f
#define c10o1 10.f
#define c11o1 11.f
#define c12o1 12.f
#define c13o1 13.f
#define c14o1 14.f
#define c15o1 15.f
#define c16o1 16.f
#define c17o1 17.f
#define c18o1 18.f
#define c21o1 21.f
#define c24o1 24.f
#define c25o1 25.f
#define c26o1 26.f
#define c27o1 27.f
#define c28o1 28.f
#define c29o1 29.f
#define c30o1 30.f
#define c32o1 32.f
#define c33o1 33.f
#define c34o1 34.f
#define c36o1 36.f
#define c40o1 40.f
#define c42o1 42.f
#define c46o1 46.f
#define c48o1 48.f
#define c50o1 50.f
#define c52o1 52.f
#define c54o1 54.f
#define c56o1 56.f
#define c64o1 64.f
#define c66o1 66.f
#define c68o1 68.f
#define c69o1 69.f
#define c72o1 72.f
#define c84o1 84.f
#define c88o1 88.f
#define c96o1 96.f
#define c100o1 100.0f
#define c130o1 130.0f
#define c152o1 152.0f
#define c166o1 166.0f
#define c195o1 195.0f
#define c216o1 216.0f
#define c264o1 264.0f
#define c290o1 290.0f
#define c367o1 367.0f

#define Op0000002 0.0000002f
#define c10eM30 1e-30
#define c10eM10 1e-10
#define smallSingle 0.0000000002f
#endif

#endif
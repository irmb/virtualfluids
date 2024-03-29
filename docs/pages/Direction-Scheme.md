<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder -->

# Direction Scheme

| **Direction** | **Index** |
|---------------|-----------|
| 000           | 0         |
| P00           | 1         |
| M00           | 2         |
| 0P0           | 3         |
| 0M0           | 4         |
| 00P           | 5         |
| 00M           | 6         |
| PP0           | 7         |
| MM0           | 8         |
| PM0           | 9         |
| MP0           | 10        |
| P0P           | 11        |
| M0M           | 12        |
| P0M           | 13        |
| M0P           | 14        |
| 0PP           | 15        |
| 0MM           | 16        |
| 0PM           | 17        |
| 0MP           | 18        |
| PPP           | 19        |
| MPP           | 20        |
| PMP           | 21        |
| MMP           | 22        |
| PPM           | 23        |
| MPM           | 24        |
| PMM           | 25        |
| MMM           | 26        |


## GPU PointerChasing

| Direction          | Index | P_MMM | P_P00 | P_0P0 | P_00P |
|--------------------|-------|-------|-------|-------|-------|
| 000                | 0     | 0     | 0     | 0     | 0     |
| P00                | 1     | 0     | 1     | 0     | 0     |
| M00                | 2     | 1     | 0     | 1     | 1     |
| 0P0                | 3     | 0     | 0     | 1     | 0     |
| 0M0                | 4     | 1     | 1     | 0     | 1     |
| 00P                | 5     | 0     | 0     | 0     | 1     |
| 00M                | 6     | 1     | 1     | 1     | 0     |
| PP0                | 7     | 0     | 1     | 1     | 0     |
| MM0                | 8     | 1     | 0     | 0     | 1     |
| PM0                | 9     | 1     | 2     | 0     | 1     |
| MP0                | 10    | 1     | 0     | 2     | 1     |
| P0P                | 11    | 0     | 1     | 0     | 1     |
| M0M                | 12    | 1     | 0     | 1     | 0     |
| P0M                | 13    | 1     | 2     | 1     | 0     |
| M0P                | 14    | 1     | 0     | 1     | 2     |
| 0PP                | 15    | 0     | 0     | 1     | 1     |
| 0MM                | 16    | 1     | 1     | 0     | 0     |
| 0PM                | 17    | 1     | 1     | 2     | 0     |
| 0MP                | 18    | 1     | 1     | 0     | 2     |
| PPP                | 19    | 0     | 1     | 1     | 1     |
| MPP                | 20    | 1     | 0     | 2     | 2     |
| PMP                | 21    | 1     | 2     | 0     | 2     |
| MMP                | 22    | 1     | 0     | 0     | 2     |
| PPM                | 23    | 1     | 2     | 2     | 0     |
| MPM                | 24    | 1     | 0     | 2     | 0     |
| PMM                | 25    | 1     | 2     | 0     | 0     |
| MMM                | 26    | 1     | 0     | 0     | 0     |

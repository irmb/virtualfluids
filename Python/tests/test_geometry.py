r"""
=======================================================================================
 ____          ____    __    ______     __________   __      __       __        __
 \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
      \    \  |    |   ________________________________________________________________
       \    \ |    |  |  ______________________________________________________________|
        \    \|    |  |  |         __          __     __     __     ______      _______
         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/

  This file is part of VirtualFluids. VirtualFluids is free software: you can
  redistribute it and/or modify it under the terms of the GNU General Public
  License as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  SPDX-License-Identifier: GPL-3.0-or-later
  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

! \author Sven Marcus, Henry Korb
=======================================================================================
"""
import unittest

from pyfluids import cpu


class TestGeometry(unittest.TestCase):

    def test_when_setting_point_coordinates_in_constructor__point_should_have_coordinates(self):
        """
        WHEN setting point coordinates in constructor THEN point should have coordinates
        """
        sut = cpu.geometry.GbPoint3D(4, 8, 3)

        self.assertEqual(sut.x1, 4)
        self.assertEqual(sut.x2, 8)
        self.assertEqual(sut.x3, 3)

    def test_when_setting_point_coordinates__point_should_have_coordinates(self):
        """
        WHEN setting point coordinates THEN point should have coordinates
        """
        sut = cpu.geometry.GbPoint3D()

        sut.x1 = 4
        sut.x2 = 8
        sut.x3 = 3

        self.assertEqual(sut.x1, 4)
        self.assertEqual(sut.x2, 8)
        self.assertEqual(sut.x3, 3)

    def test_when_setting_line_points__line_should_have_points(self):
        """
        WHEN setting line points THEN line should have points
        """
        sut = cpu.geometry.GbLine3D()

        point1 = cpu.geometry.GbPoint3D()
        point2 = cpu.geometry.GbPoint3D()
        sut.point1 = point1
        sut.point2 = point2

        self.assertEqual(sut.point1, point1)
        self.assertEqual(sut.point2, point2)

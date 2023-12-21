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


class BoundaryConditionsTest(unittest.TestCase):

    def test__can_create_no_slip_bc(self):
        """
        Should be able to create NoSlipBoundaryCondition
        """
        sut = cpu.boundaryconditions.NoSlipBoundaryCondition()

    def test__can_create_velocity_bc(self):
        """
        Should be able to create VelocityBoundaryCondition
        """
        sut = cpu.boundaryconditions.VelocityBoundaryCondition()

    def test__can_create_velocity_bc_with_directions_function_and_time(self):
        """
        Should be able to create VelocityBoundaryCondition with directions, function and start/end time
        """
        from pymuparser import Parser

        parser = Parser()
        parser.expression = "1"
        sut = cpu.boundaryconditions.VelocityBoundaryCondition(True, True, True, parser, 0, 1)

    def test__can_create_velocity_bc_with_directions__function_per_direction__and__time(self):
        """
        Should be able to create VelocityBoundaryCondition with directions, function per direction and start/end time
        """
        from pymuparser import Parser

        f1 = Parser()
        f1.expression = "1"

        f2 = Parser()
        f2.expression = "1"

        f3 = Parser()
        f3.expression = "1"
        sut = cpu.boundaryconditions.VelocityBoundaryCondition(True, True, True, f1, f2, f3, 0, 1)

    def test__can_create_velocity_bc_with_speeds_and_times_per_direction(self):
        """
        Should be able to create VelocityBoundaryCondition with speeds and start/end times per direction
        """
        vx1, vx2, vx3 = 1, 2, 3
        start1, end1 = 0, 1
        start2, end2 = 1, 2
        start3, end3 = 2, 3

        sut = cpu.boundaryconditions.VelocityBoundaryCondition(vx1, start1, end1, vx2, start2, end2, vx3, start3, end3)

    def test__can_create_non_reflecting_outflow(self):
        """
        Should be able to create NonReflectingOutflow
        """

        sut = cpu.boundaryconditions.NonReflectingOutflow()

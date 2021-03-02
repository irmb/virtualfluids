import unittest
from pyfluids.boundaryconditions import *


class BoundaryConditionsTest(unittest.TestCase):

    def test__can_create_no_slip_bc(self):
        """
        Should be able to create NoSlipBoundaryCondition
        """
        sut = NoSlipBoundaryCondition()

    def test__can_create_velocity_bc(self):
        """
        Should be able to create VelocityBoundaryCondition
        """
        sut = VelocityBoundaryCondition()

    def test__can_create_velocity_bc_with_directions_function_and_time(self):
        """
        Should be able to create VelocityBoundaryCondition with directions, function and start/end time
        """
        from pymuparser import Parser

        parser = Parser()
        parser.expression = "1"
        sut = VelocityBoundaryCondition(True, True, True, parser, 0, 1)

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
        sut = VelocityBoundaryCondition(True, True, True, f1, f2, f3, 0, 1)

    def test__can_create_velocity_bc_with_speeds_and_times_per_direction(self):
        """
        Should be able to create VelocityBoundaryCondition with speeds and start/end times per direction
        """
        vx1, vx2, vx3 = 1, 2, 3
        start1, end1 = 0, 1
        start2, end2 = 1, 2
        start3, end3 = 2, 3

        sut = VelocityBoundaryCondition(vx1, start1, end1, vx2, start2, end2, vx3, start3, end3)

    def test__can_create_non_reflecting_outflow(self):
        """
        Should be able to create NonReflectingOutflow
        """

        sut = NonReflectingOutflow()
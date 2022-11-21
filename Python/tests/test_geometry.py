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

import unittest

from pyfluids.cpu.kernel import LBMKernel, KernelType


class TestLBMKernel(unittest.TestCase):

    def setUp(self) -> None:
        self.sut = LBMKernel(KernelType.BGK)

    def test_lbm_kernel__when_use_forcing_set_to_true__use_forcing_should_be_true(self) -> None:
        """
        WHEN use_forcing is set to true THEN use_forcing should be true
        """
        self.sut.use_forcing = True

        self.assertTrue(self.sut.use_forcing)

    def test_lbm_kernel__when_forcing_in_x1_set_to_five__forcing_in_x1_should_be_five(self) -> None:
        """
        WHEN forcing_in_x1 is set to 5 THEN forcing_in_x1 should be 5
        """
        self.sut.forcing_in_x1 = 5

        self.assertEqual(self.sut.forcing_in_x1, 5)

    def test_lbm_kernel__when_forcing_in_x2_set_to_five__forcing_in_x2_should_be_five(self) -> None:
        """
        WHEN forcing_in_x2 is set to 5 THEN forcing_in_x2 should be 5
        """
        self.sut.forcing_in_x2 = 5

        self.assertEqual(self.sut.forcing_in_x2, 5)

    def test_lbm_kernel__when_forcing_in_x3_set_to_five__forcing_in_x3_should_be_five(self) -> None:
        """
        WHEN forcing_in_x3 is set to 5 THEN forcing_in_x3 should be 5
        """
        self.sut.forcing_in_x3 = 5

        self.assertEqual(self.sut.forcing_in_x3, 5)

    def test_lbm_kernel__when_setting_forcing_in_all_directions__forcing_should_equal_set_values(self) -> None:
        """
        WHEN setting forcing in all directions THEN forcing should equal set values
        """

        self.sut.set_forcing(3, 8, 5)

        self.assertEqual(self.sut.forcing_in_x1, 3)
        self.assertEqual(self.sut.forcing_in_x2, 8)
        self.assertEqual(self.sut.forcing_in_x3, 5)

    def test_lbm_kernel__when_getting_type__should_equal_kernel_type_enum_value(self) -> None:
        """
        WHEN getting the kernel type IT should equal the corresponding KernelType enum value
        """

        actual = self.sut.type
        self.assertEqual(KernelType.BGK, actual)

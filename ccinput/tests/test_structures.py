from unittest import TestCase

from ccinput.utilities import standardize_xyz
from ccinput.exceptions import *

class XyzTests(TestCase):
    def test_standard_xyz(self):
        xyz = "1\n\nCl 0.0 0.0 0.0\n"
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")

    def test_headerless_xyz(self):
        xyz = "Cl 0.0 0.0 0.0\n"
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")

    def test_extra_line(self):
        xyz = "1\n\nCl 0.0 0.0 0.0\n\n"
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")

    def test_extra_lines(self):
        xyz = "1\n\nCl 0.0 0.0 0.0\n\n\n\n"
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")

    def test_invalid_header(self):
        xyz = "2\n\nCl 0.0 0.0 0.0\n"
        with self.assertRaises(InvalidXYZ):
            standardize_xyz(xyz)

    def test_invalid_header_space(self):
        xyz = "2\n\nCl 0.0 0.0 0.0\n\n"
        with self.assertRaises(InvalidXYZ):
            standardize_xyz(xyz)

    def test_extra_spacing(self):
        xyz = """1
        header
        Cl 0.0 0.0 0.0


        """
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")

    def test_multiple_atoms(self):
        xyz = """3
        header
        Cl 0.0 0.0 0.0
        Cl 0.0 0.0 0.0
        Cl 0.0 0.0 0.0
        """
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\nCl   0.00000000   0.00000000   0.00000000\nCl   0.00000000   0.00000000   0.00000000\n")

    def test_multiple_atoms_invalid_header(self):
        xyz = """2
        header
        Cl 0.0 0.0 0.0
        Cl 0.0 0.0 0.0
        Cl 0.0 0.0 0.0
        """
        with self.assertRaises(InvalidXYZ):
            standardize_xyz(xyz)

    def test_multiple_atoms_invalid_header2(self):
        xyz = """4
        header
        Cl 0.0 0.0 0.0
        Cl 0.0 0.0 0.0
        Cl 0.0 0.0 0.0
        """
        with self.assertRaises(InvalidXYZ):
            standardize_xyz(xyz)

    def test_invalid_element(self):
        xyz = "1\n\nBl 0.0 0.0 0.0\n"
        with self.assertRaises(InvalidXYZ):
            standardize_xyz(xyz)

    def test_lowercase(self):
        xyz = "1\n\ncl 0.0 0.0 0.0\n"
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")

    def test_uppercase(self):
        xyz = "1\n\nCL 0.0 0.0 0.0\n"
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")


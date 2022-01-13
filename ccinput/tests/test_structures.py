from unittest import TestCase

from ccinput.utilities import standardize_xyz
from ccinput.exceptions import *

class XyzTests(TestCase):
    def test_standard_xyz(self):
        xyz = "1\n\nCl 0.0 0.0 0.0\n"
        self.assertEqual(standardize_xyz(xyz), "Cl   0.00000000   0.00000000   0.00000000\n")

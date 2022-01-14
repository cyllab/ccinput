from unittest import TestCase

from ccinput.calculation import Calculation, Parameters
from ccinput.constants import CalcType
from ccinput.exceptions import *

class CalculationTests(TestCase):
    def setUp(self):
        self.xyz = "Cl 0 0 0\n"
        self.params = Parameters("gaussian", method="am1")

    def create_calc(self, type=CalcType.SP, nproc=1, mem="1000MB", charge=-1, multiplicity=1):
        return Calculation(self.xyz, self.params, type, nproc=nproc, mem=mem, \
                charge=charge, multiplicity=multiplicity, software="gaussian")

    def test_base(self):
        calc = self.create_calc()

    def test_correct_str_charge(self):
        calc = self.create_calc(charge="-1")

    def test_correct_str_charge2(self):
        calc = self.create_calc(charge="+1")

    def test_incorrect_str_charge(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(charge="a")

    def test_incorrect_str_charge2(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(charge="-")

    def test_incorrect_float_charge(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(charge=1.5)

    def test_incorrect_float_charge2(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(charge=0.001)

    def test_empty_str_charge(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(charge="")

    def test_correct_str_multiplicity(self):
        calc = self.create_calc(multiplicity="1")

    def test_correct_str_multiplicity2(self):
        calc = self.create_calc(multiplicity="+1")

    def test_empty_str_multiplicity(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(multiplicity="")

    def test_incorrect_str_multiplicity(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(multiplicity="a")

    def test_incorrect_str_multiplicity2(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(multiplicity="-")

    def test_incorrect_float_multiplicity(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(multiplicity=1.5)

    def test_incorrect_float_multiplicity2(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(multiplicity=0.001)

    def test_correct_str_nproc(self):
        calc = self.create_calc(nproc="1")

    def test_correct_str_nproc2(self):
        calc = self.create_calc(nproc="+1")

    def test_incorrect_str_nproc(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(nproc="a")

    def test_incorrect_str_nproc2(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(nproc="-1")

    def test_incorrect_str_nproc2(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(nproc="-")

    def test_incorrect_int_nproc(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(nproc=-1)

    def test_incorrect_float_nproc(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(nproc=1.5)

    def test_incorrect_float_nproc2(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(nproc=0.001)

    def test_memory_parsing_int(self):
        calc = self.create_calc(mem=1000)
        self.assertEqual(calc.mem, 1000)

    def test_memory_parsing_str_unitless(self):
        calc = self.create_calc(mem="1000")
        self.assertEqual(calc.mem, 1000)

    def test_memory_parsing_str_m(self):
        calc = self.create_calc(mem="1000m")
        self.assertEqual(calc.mem, 1000)

    def test_memory_parsing_str_space_m(self):
        calc = self.create_calc(mem="1000 m")
        self.assertEqual(calc.mem, 1000)

    def test_memory_parsing_str_mib(self):
        calc = self.create_calc(mem="1000 mib")
        self.assertEqual(calc.mem, 1049)

    def test_memory_parsing_str_gib(self):
        calc = self.create_calc(mem="1 gib")
        self.assertEqual(calc.mem, 1074)

    def test_memory_parsing_str_g(self):
        calc = self.create_calc(mem="1000g")
        self.assertEqual(calc.mem, 1000000)

    def test_memory_parsing_str_space_g(self):
        calc = self.create_calc(mem="1000 g")
        self.assertEqual(calc.mem, 1000000)

    def test_memory_parsing_str_t(self):
        calc = self.create_calc(mem="1t")
        self.assertEqual(calc.mem, 1000000)

    def test_memory_parsing_str_space_t(self):
        calc = self.create_calc(mem="1 t")
        self.assertEqual(calc.mem, 1000000)

    def test_memory_parsing_negative_str(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(mem="-1 g")

    def test_memory_parsing_negative_int(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(mem=-1000)

    def test_memory_parsing_too_much(self):
        with self.assertRaises(InvalidParameter):
            calc = self.create_calc(mem="1000 t")

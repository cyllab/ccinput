from unittest import TestCase

from ccinput.calculation import Calculation, Parameters, Constraint, parse_scan_constraints
from ccinput.constants import CalcType
from ccinput.exceptions import InvalidParameter

class CalculationTests(TestCase):
    def setUp(self):
        self.xyz = "Cl 0 0 0\n"
        self.params = Parameters("gaussian", method="am1")

    def create_calc(self, type=CalcType.SP, nproc=1, mem="1000MB", charge=-1, multiplicity=1):
        return Calculation(self.xyz, self.params, type, nproc=nproc, mem=mem, \
                charge=charge, multiplicity=multiplicity, software="gaussian")

    def test_base(self):
        self.create_calc()

    def test_correct_str_charge(self):
        self.create_calc(charge="-1")

    def test_correct_str_charge2(self):
        self.create_calc(charge="+1")

    def test_incorrect_str_charge(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(charge="a")

    def test_incorrect_str_charge2(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(charge="-")

    def test_incorrect_float_charge(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(charge=1.5)

    def test_incorrect_float_charge2(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(charge=0.001)

    def test_empty_str_charge(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(charge="")

    def test_correct_str_multiplicity(self):
        self.create_calc(multiplicity="1")

    def test_correct_str_multiplicity2(self):
        self.create_calc(multiplicity="+1")

    def test_empty_str_multiplicity(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(multiplicity="")

    def test_incorrect_str_multiplicity(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(multiplicity="a")

    def test_incorrect_str_multiplicity2(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(multiplicity="-")

    def test_incorrect_float_multiplicity(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(multiplicity=1.5)

    def test_incorrect_float_multiplicity2(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(multiplicity=0.001)

    def test_correct_str_nproc(self):
        self.create_calc(nproc="1")

    def test_correct_str_nproc2(self):
        self.create_calc(nproc="+1")

    def test_incorrect_str_nproc(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(nproc="a")

    def test_incorrect_str_nproc2(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(nproc="-1")

    def test_incorrect_str_nproc3(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(nproc="-")

    def test_incorrect_int_nproc(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(nproc=-1)

    def test_incorrect_float_nproc(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(nproc=1.5)

    def test_incorrect_float_nproc2(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(nproc=0.001)

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
            self.create_calc(mem="-1 g")

    def test_memory_parsing_negative_int(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(mem=-1000)

    def test_memory_parsing_too_much(self):
        with self.assertRaises(InvalidParameter):
            self.create_calc(mem="1000 t")

class ScanParsingTests(TestCase):
    def setUp(self):
        self.xyz = "Cl 0 0 0\nCl 1 0 0\nCl 2 0 0\nCl 3 0 0\nCl 4 0 0\n"

    def test_blank(self):
        scan = parse_scan_constraints(arr=[], sfrom=[], sto=[], snsteps=[],
                                      sstep=[], xyz_str=self.xyz, software="orca")
        self.assertEqual(scan, [])

    def test_correct_distance_nsteps(self):
        scan = parse_scan_constraints(arr=[['1', '2']], sfrom=[1], sto=[1.5], snsteps=[5],
                                      sstep=[], xyz_str=self.xyz, software="orca")
        self.assertTrue(len(scan) == 1)
        constr = scan[0]
        self.assertTrue(isinstance(constr, Constraint))
        self.assertEqual(constr.ids, [1, 2])
        self.assertEqual(constr.num_steps, 5)
        self.assertAlmostEqual(constr.start_d, 1.0)
        self.assertAlmostEqual(constr.end_d, 1.5)
        self.assertAlmostEqual(constr.step_size, 0.1)

    def test_correct_distance_step(self):
        scan = parse_scan_constraints(arr=[['1', '2']], sfrom=[1], sto=[1.5], snsteps=[],
                                      sstep=[0.1], xyz_str=self.xyz, software="orca")
        self.assertTrue(len(scan) == 1)
        constr = scan[0]
        self.assertTrue(isinstance(constr, Constraint))
        self.assertEqual(constr.ids, [1, 2])
        self.assertEqual(constr.num_steps, 5)
        self.assertAlmostEqual(constr.start_d, 1.0)
        self.assertAlmostEqual(constr.end_d, 1.5)
        self.assertAlmostEqual(constr.step_size, 0.1)

    def test_correct_distance_decrease(self):
        scan = parse_scan_constraints(arr=[['1', '2']], sfrom=[1.5], sto=[1], snsteps=[],
                                      sstep=[0.1], xyz_str=self.xyz, software="orca")
        self.assertTrue(len(scan) == 1)
        constr = scan[0]
        self.assertTrue(isinstance(constr, Constraint))
        self.assertEqual(constr.ids, [1, 2])
        self.assertEqual(constr.num_steps, 5)
        self.assertAlmostEqual(constr.start_d, 1.5)
        self.assertAlmostEqual(constr.end_d, 1.0)
        self.assertAlmostEqual(constr.step_size, -0.1)

    def test_correct_distance_start_gaussian(self):
        scan = parse_scan_constraints(arr=[['1', '2']], sfrom=[1.5], sto=[0.5], snsteps=[],
                                      sstep=[0.1], xyz_str=self.xyz, software="gaussian")
        self.assertTrue(len(scan) == 1)
        constr = scan[0]
        self.assertTrue(isinstance(constr, Constraint))
        self.assertEqual(constr.ids, [1, 2])
        self.assertEqual(constr.num_steps, 5)
        self.assertAlmostEqual(constr.start_d, 1.0)
        self.assertAlmostEqual(constr.end_d, 0.5)
        self.assertAlmostEqual(constr.step_size, -0.1)

    def test_correct_distance_multiple_same(self):
        scan = parse_scan_constraints(arr=[['1', '2'], ['3', '4']], sfrom=[1, 1], sto=[1.5, 1.5],
                                      snsteps=[], sstep=[0.1, 0.1], xyz_str=self.xyz,
                                      software="orca")
        self.assertTrue(len(scan) == 2)
        for constr in scan:
            self.assertTrue(isinstance(constr, Constraint))
            self.assertEqual(constr.num_steps, 5)
            self.assertAlmostEqual(constr.start_d, 1.0)
            self.assertAlmostEqual(constr.end_d, 1.5)
            self.assertAlmostEqual(constr.step_size, 0.1)

    def test_correct_distance_multiple_different(self):
        scan = parse_scan_constraints(arr=[['1', '2'], ['3', '4']], sfrom=[1, 1], sto=[1.5, 0.5],
                                      snsteps=[], sstep=[0.1, -0.1], xyz_str=self.xyz,
                                      software="orca")
        self.assertTrue(len(scan) == 2)
        for constr in scan:
            self.assertTrue(isinstance(constr, Constraint))
            self.assertEqual(constr.num_steps, 5)

        constr1, constr2 = scan
        self.assertAlmostEqual(constr1.start_d, 1.0)
        self.assertAlmostEqual(constr1.end_d, 1.5)
        self.assertAlmostEqual(constr1.step_size, 0.1)

        self.assertAlmostEqual(constr2.start_d, 1.0)
        self.assertAlmostEqual(constr2.end_d, 0.5)
        self.assertAlmostEqual(constr2.step_size, -0.1)

    def test_missing_ids(self):
        with self.assertRaises(InvalidParameter):
            parse_scan_constraints(arr=[['1', '2']], sfrom=[1, 1], sto=[1.5, 0.5],
                                      snsteps=[], sstep=[0.1, -0.1], xyz_str=self.xyz,
                                      software="orca")
    def test_missing_nsteps(self):
        with self.assertRaises(InvalidParameter):
            parse_scan_constraints(arr=[['1', '2'], ['3', '4']], sfrom=[1, 1], sto=[1.5, 0.5],
                                      snsteps=[5], sstep=[], xyz_str=self.xyz,
                                      software="orca")
    def test_missing_to(self):
        with self.assertRaises(InvalidParameter):
            parse_scan_constraints(arr=[['1', '2'], ['3', '4']], sfrom=[1, 1], sto=[1.5],
                                      snsteps=[], sstep=[0.1, -0.1], xyz_str=self.xyz,
                                      software="orca")
    def test_missing_step(self):
        with self.assertRaises(InvalidParameter):
            parse_scan_constraints(arr=[['1', '2'], ['3', '4']], sfrom=[1, 1], sto=[1.5, 0.5],
                                      snsteps=[], sstep=[0.1], xyz_str=self.xyz,
                                      software="orca")

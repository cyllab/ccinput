from ccinput.tests.testing_utilities import InputTests
from ccinput.packages.xtb import XtbCalculation
from ccinput.exceptions import InvalidParameter


class XtbTests(InputTests):
    def test_sp_basic(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "xtb",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_sp_charge(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "xtb",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --chrg -1"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_sp_multiplicity(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "xtb",
            "multiplicity": "2",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --uhf 2"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_sp_charge_multiplicity(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "xtb",
            "multiplicity": "3",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --chrg -1 --uhf 3"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_opt_charge(self):
        params = {
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "xtb",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --opt tight --chrg -1"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_freq_charge(self):
        params = {
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "xtb",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --hess --chrg -1"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_opt_freq(self):
        params = {
            "type": "Opt+Freq",
            "file": "Cl.xyz",
            "software": "xtb",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --ohess --chrg -1"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_solvent(self):
        params = {
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "xtb",
            "solvation_model": "GBSA",
            "solvent": "chcl3",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --hess -g chcl3 --chrg -1"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_solvent_ALPB(self):
        params = {
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "xtb",
            "solvent": "chcl3",
            "solvation_model": "ALPB",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --hess --alpb chcl3 --chrg -1"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_solvent_invalid_PCM(self):
        params = {
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "xtb",
            "solvent": "chcl3",
            "solvation_model": "PCM",
            "charge": "-1",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_solvent_synonym(self):
        params = {
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "xtb",
            "solvent": "chloroform",
            "solvation_model": "GBSA",
            "charge": "-1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb Cl.xyz --hess -g chcl3 --chrg -1"

        self.assertTrue(self.is_equivalent(REF, xtb.command))
        self.assertTrue(self.is_equivalent("", xtb.input_file))

    def test_solvent_invalid(self):
        params = {
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "xtb",
            "solvent": "octanol",
            "solvation_model": "GBSA",
            "charge": "-1",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_scan(self):
        params = {
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Scan_9_1.4_10/1_2;",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --opt tight --input input"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=1.0
        distance: 1, 2, auto
        $scan
        1: 9.0, 1.4, 10
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_constrained_opt_no_constraint(self):
        params = {
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_freeze(self):
        params = {
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/1_2;",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --opt tight --input input"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=1.0
        distance: 1, 2, auto
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_constraint_overlap(self):
        params = {
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/1_2;Freeze/2_3;",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --opt tight --input input"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=1.0
        distance: 1, 2, auto
        distance: 2, 3, auto
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_freeze_soft(self):
        params = {
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/1_2;",
            "specifications": "--forceconstant 0.1",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --opt tight --input input"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=0.1
        distance: 1, 2, auto
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_duplicate_specifications(self):
        params = {
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/1_2;",
            "specifications": "--forceconstant 0.1 --forceconstant 0.2",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --opt tight --input input"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=0.2
        distance: 1, 2, auto
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_conformational_search(self):
        params = {
            "type": "Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -rthr 0.6 -ewin 6"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        self.assertEqual("", xtb.input_file)

    def test_conformational_search_specs(self):
        params = {
            "type": "Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--rthr 0.8 --ewin 8",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -rthr 0.8 -ewin 8"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        self.assertEqual("", xtb.input_file)

    def test_conformational_search_nci(self):
        params = {
            "type": "Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--rthr 0.8 --nci",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -nci -rthr 0.8 -ewin 6"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

        self.assertEqual("", xtb.input_file)

    def test_constrained_conformational_search1(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/1_2;",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -cinp input -rthr 0.6 -ewin 6"
        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=1.0
        reference=ethanol.xyz
        distance: 1, 2, auto
        atoms: 1-2
        $metadyn
        atoms: 3-9
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_constrained_conformational_search2(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/1_4;Freeze/6_8;",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -cinp input -rthr 0.6 -ewin 6"
        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=1.0
        reference=ethanol.xyz
        distance: 1, 4, auto
        distance: 6, 8, auto
        atoms: 1,4,6,8
        $metadyn
        atoms: 2-3,5,7,9
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_constrained_conformational_search3(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/2_3;",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -cinp input -rthr 0.6 -ewin 6"
        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=1.0
        reference=ethanol.xyz
        distance: 2, 3, auto
        atoms: 2-3
        $metadyn
        atoms: 1,4-9
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_constrained_conformational_search4(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/2_3;",
            "specifications": "--force_constant 2.0",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -cinp input -rthr 0.6 -ewin 6"
        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=2.0
        reference=ethanol.xyz
        distance: 2, 3, auto
        atoms: 2-3
        $metadyn
        atoms: 1,4-9
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_constrained_conformational_search5(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/2_3;Freeze/3_4;",
            "specifications": "--force_constant 2.0",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -cinp input -rthr 0.6 -ewin 6"
        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=2.0
        reference=ethanol.xyz
        distance: 2, 3, auto
        distance: 3, 4, auto
        atoms: 2-4
        $metadyn
        atoms: 1,5-9
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_constrained_conformational_search_no_constraint(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--force_constant 2.0",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_constrained_conformational_search_equals(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/2_3;",
            "specifications": "--force_constant=2.0",
        }

        xtb = self.generate_calculation(**params)

        REF = "crest ethanol.xyz -cinp input -rthr 0.6 -ewin 6"
        self.assertTrue(self.is_equivalent(REF, xtb.command))

        INPUT = """$constrain
        force constant=2.0
        reference=ethanol.xyz
        distance: 2, 3, auto
        atoms: 2-3
        $metadyn
        atoms: 1,4-9
        """
        self.assertTrue(self.is_equivalent(INPUT, xtb.input_file))

    def test_invalid_specification(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/2_3;",
            "specifications": "--force 2.0",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_invalid_specification2(self):
        params = {
            "type": "Constrained Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "constraints": "Freeze/2_3;",
            "specifications": "-force_constant 2.0",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_invalid_specification3(self):
        params = {
            "type": "Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "-rthr 0.8 --ewin 8",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_invalid_specification4(self):
        params = {
            "type": "Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--rthr abc",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_invalid_specification_for_type1(self):
        params = {
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--rthr 0.8",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_invalid_specification_for_type2(self):
        params = {
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--ewin 8",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_unknown_specification(self):
        params = {
            "type": "Conformational Search",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--abc",
        }

        with self.assertRaises(InvalidParameter):
            xtb = self.generate_calculation(**params)

    def test_gfn0(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--gfn 0",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --gfn 0"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

    def test_gfn0_2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--gfn0",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --gfn 0"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

    def test_gfn1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--gfn 0",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --gfn 0"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

    def test_gfn2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--gfn 2",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

    def test_gfnff(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "xtb",
            "specifications": "--gfnff",
        }

        xtb = self.generate_calculation(**params)

        REF = "xtb ethanol.xyz --gfnff"

        self.assertTrue(self.is_equivalent(REF, xtb.command))

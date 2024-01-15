from ccinput.tests.testing_utilities import InputTests
from ccinput.exceptions import (
    InvalidParameter,
    ImpossibleCalculation,
    UnimplementedError,
    MissingParameter,
)


class NwchemTests(InputTests):
    def test_sp_HF(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
        }

        with self.assertRaises(UnimplementedError):
            self.generate_calculation(**params)

    def test_sp_DFT_SMD(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "B3LYP",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc b3lyp
        mult 1
        end

        cosmo
        minbem 3
        ificos 1
        solvent chcl3
        do_cosmo_smd
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_SMD18(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
            "solvation_radii": "SMD18",
        }

        with self.assertRaises(UnimplementedError):
            self.generate_calculation(**params)

    def test_sp_DFT_COSMO(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "camb3lyp",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "c6h6",
            "solvation_model": "COSMO",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.00
        cam 0.33 cam_alpha 0.19 cam_beta 0.46
        mult 1
        end

        cosmo
        minbem 3
        ificos 1
        solvent benzene
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_hf_COSMO(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "c6h6",
            "solvation_model": "COSMO",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        cosmo
        minbem 3
        ificos 1
        solvent benzene
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_hf_CPCM(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "def2tzvp",
            "charge": "-1",
            "solvent": "c6h6",
            "solvation_model": "CPCM",
        }

        with self.assertRaises(UnimplementedError):
            self.generate_calculation(**params)

    def test_solvent_synonyms1(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "chloroform",
            "solvation_model": "COSMO",
            "solvation_radii": "default",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        cosmo
        minbem 3
        ificos 1
        solvent chcl3
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvent_synonyms2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "meoh",
            "solvation_model": "Cosmo",
            "solvation_radii": "default",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        cosmo
        minbem 3
        ificos 1
        solvent methanol
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_solvation(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "ABC",
        }

        with self.assertRaises(Exception):
            self.generate_calculation(**params)

    def test_sp_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "scf(grid coarse)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        grid coarse
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_specifications2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "dft(grid coarse)",
        }

        inp = self.generate_calculation(**params)

        REF = """
            TITLE "File created by ccinput"
            start Cl
            memory total 1250 mb
            charge -1

            geometry units angstroms noautosym
            Cl   0.00000000   0.00000000   0.00000000
            end

            basis
            * library Def2-SVP
            end

            dft
            xc m06-2x
            mult 1
            grid coarse
            end

            task dft energy
            """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_superfluous_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "opt(loose)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.calculation_block.strip(), "")

    def test_superfluous_specifications2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "freq(noraman)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.calculation_block.strip(), "")

    def test_opt_HF(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        task scf optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_HF(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        task scf freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_bond_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Freeze/1_2;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        geometry adjust
        zcoord
        bond 1 2 constant
        end
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_angle_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Freeze/2_1_3;",
        }
        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        geometry adjust
        zcoord
        angle 2 1 3 constant
        end
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_opt_mod(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "",
        }

        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_no_method(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "",
            "basis_set": "6-31+G(d,p)",
            "constraints": "",
        }

        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_freeze_dihedral_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Freeze/4_1_5_8;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        geometry adjust
        zcoord
        torsion 4 1 5 8 constant
        end
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_dihedral_DFT2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Freeze/4_1_5_8;Freeze/1_2_3_4;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        geometry adjust
        zcoord
        torsion 4 1 5 8 constant
        torsion 1 2 3 4 constant
        end
        end

        task dft optimize


        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nmr_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "NMR Prediction",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        property
        shielding
        end

        task dft property
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "TS Optimisation",
            "file": "mini_ts.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start mini_ts
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        N    1.08764072  -0.33994563  -0.00972525
        H    1.99826837   0.05502843   0.00651241
        H    0.59453997  -0.48560162   0.83949232
        H    0.66998094  -0.58930117  -0.87511947
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft saddle
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_gen_bs(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "O=Def2-TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis spherical
        #BASIS SET: (12s,7p,3d,1f) -> [6s,4p,3d,1f]
        O    S
        27032.3826310              0.21726302465E-03
        4052.3871392              0.16838662199E-02
        922.32722710             0.87395616265E-02
        261.24070989             0.35239968808E-01
        85.354641351            0.11153519115
        31.035035245            0.25588953961
        O    S
        12.260860728            0.39768730901
        4.9987076005           0.24627849430
        O    S
        1.1703108158           1.0000000
        O    S
        0.46474740994          1.0000000
        O    S
        0.18504536357          1.0000000
        O    S
        0.70288026270E-01      1.0000000
        O    P
        63.274954801            0.60685103418E-02
        14.627049379            0.41912575824E-01
        4.4501223456           0.16153841088
        1.5275799647           0.35706951311
        O    P
        0.52935117943           .44794207502
        O    P
        0.17478421270           .24446069663
        O    P
        0.51112745706E-01      1.0000000
        O    D
        2.31400000             1.0000000
        O    D
        0.64500000             1.0000000
        O    D
        0.14696477366          1.0000000
        O    F
        1.42800000             1.0000000

        * library 6-31+G* except O
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_gen_bs_space(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": " O=Def2-TZVPD; ",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis spherical
        #BASIS SET: (12s,7p,3d,1f) -> [6s,4p,3d,1f]
        O    S
        27032.3826310              0.21726302465E-03
        4052.3871392              0.16838662199E-02
        922.32722710             0.87395616265E-02
        261.24070989             0.35239968808E-01
        85.354641351            0.11153519115
        31.035035245            0.25588953961
        O    S
        12.260860728            0.39768730901
        4.9987076005           0.24627849430
        O    S
        1.1703108158           1.0000000
        O    S
        0.46474740994          1.0000000
        O    S
        0.18504536357          1.0000000
        O    S
        0.70288026270E-01      1.0000000
        O    P
        63.274954801            0.60685103418E-02
        14.627049379            0.41912575824E-01
        4.4501223456           0.16153841088
        1.5275799647           0.35706951311
        O    P
        0.52935117943           .44794207502
        O    P
        0.17478421270           .24446069663
        O    P
        0.51112745706E-01      1.0000000
        O    D
        2.31400000             1.0000000
        O    D
        0.64500000             1.0000000
        O    D
        0.14696477366          1.0000000
        O    F
        1.42800000             1.0000000

        * library 6-31+G* except O
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_gen_bs_space_middle(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "O = Def2 -T ZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis spherical
        #BASIS SET: (12s,7p,3d,1f) -> [6s,4p,3d,1f]
        O    S
        27032.3826310              0.21726302465E-03
        4052.3871392              0.16838662199E-02
        922.32722710             0.87395616265E-02
        261.24070989             0.35239968808E-01
        85.354641351            0.11153519115
        31.035035245            0.25588953961
        O    S
        12.260860728            0.39768730901
        4.9987076005           0.24627849430
        O    S
        1.1703108158           1.0000000
        O    S
        0.46474740994          1.0000000
        O    S
        0.18504536357          1.0000000
        O    S
        0.70288026270E-01      1.0000000
        O    P
        63.274954801            0.60685103418E-02
        14.627049379            0.41912575824E-01
        4.4501223456           0.16153841088
        1.5275799647           0.35706951311
        O    P
        0.52935117943           .44794207502
        O    P
        0.17478421270           .24446069663
        O    P
        0.51112745706E-01      1.0000000
        O    D
        2.31400000             1.0000000
        O    D
        0.64500000             1.0000000
        O    D
        0.14696477366          1.0000000
        O    F
        1.42800000             1.0000000

        * library 6-31+G* except O
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_irrelevant_gen_bs(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "Cl=Def2-TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library 6-31+G*
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_gen_bs_synonym(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "O=Def2TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis spherical
        #BASIS SET: (12s,7p,3d,1f) -> [6s,4p,3d,1f]
        O    S
        27032.3826310              0.21726302465E-03
        4052.3871392              0.16838662199E-02
        922.32722710             0.87395616265E-02
        261.24070989             0.35239968808E-01
        85.354641351            0.11153519115
        31.035035245            0.25588953961
        O    S
        12.260860728            0.39768730901
        4.9987076005           0.24627849430
        O    S
        1.1703108158           1.0000000
        O    S
        0.46474740994          1.0000000
        O    S
        0.18504536357          1.0000000
        O    S
        0.70288026270E-01      1.0000000
        O    P
        63.274954801            0.60685103418E-02
        14.627049379            0.41912575824E-01
        4.4501223456           0.16153841088
        1.5275799647           0.35706951311
        O    P
        0.52935117943           .44794207502
        O    P
        0.17478421270           .24446069663
        O    P
        0.51112745706E-01      1.0000000
        O    D
        2.31400000             1.0000000
        O    D
        0.64500000             1.0000000
        O    D
        0.14696477366          1.0000000
        O    F
        1.42800000             1.0000000

        * library 6-31+G* except O
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_genecp_bs(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Ph2I_cation.xyz",
            "software": "nwchem",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=Def2-TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Ph2I_cation
        memory total 1250 mb
        charge 1

        geometry units angstroms noautosym
        C   -3.06870000  -2.28540000   0.00000000
        C   -1.67350000  -2.28540000   0.00000000
        C   -0.97600000  -1.07770000   0.00000000
        C   -1.67360000   0.13090000  -0.00120000
        C   -3.06850000   0.13080000  -0.00170000
        C   -3.76610000  -1.07740000  -0.00070000
        H   -3.61840000  -3.23770000   0.00040000
        H   -1.12400000  -3.23790000   0.00130000
        H    0.12370000  -1.07760000   0.00060000
        H   -1.12340000   1.08300000  -0.00130000
        H   -4.86570000  -1.07720000  -0.00090000
        I   -4.11890000   1.94920000  -0.00350000
        C   -4.64360000   2.85690000  -1.82310000
        C   -3.77180000   3.76300000  -2.42740000
        C   -5.86360000   2.55380000  -2.42750000
        C   -4.12020000   4.36650000  -3.63560000
        H   -2.81040000   4.00240000  -1.95030000
        C   -6.21180000   3.15650000  -3.63650000
        H   -6.55070000   1.83950000  -1.95140000
        C   -5.34050000   4.06290000  -4.24060000
        H   -3.43340000   5.08120000  -4.11170000
        H   -7.17360000   2.91710000  -4.11310000
        H   -5.61500000   4.53870000  -5.19320000
        end

        basis spherical
        #BASIS SET: (12s,11p,9d,2f) -> [7s,6p,4d,2f]
        I    S
        5899.5791533              0.24188269271E-03
        898.54238765             0.15474041742E-02
        200.37237912             0.42836684457E-02
        31.418053840           -0.39417936275E-01
        15.645987838            0.96086691992
        I    S
        11.815741857            0.75961524091
        6.4614458287           0.42495501835
        I    S
        2.3838067579           1.0000000
        I    S
        1.1712089662           1.0000000
        I    S
        0.32115875757          1.0000000
        I    S
        0.12387919364          1.0000000
        I    S
        0.43491550641E-01      1.0000000
        I    P
        197.30030547             0.73951226905E-03
        20.061411349            0.66168450008E-01
        9.7631460485          -0.28554662348
        I    P
        12.984316904           -0.49096186164E-01
        3.6199503008           0.38914432482
        2.0232273090           0.65610817262
        1.0367490559           0.31803551647
        I    P
        0.45937816000          1.0000000
        I    P
        0.19116532928          1.0000000
        I    P
        0.74878813023E-01      1.0000000
        I    P
        0.21653491846E-01      1.0000000
        I    D
        119.12671745             0.82596039573E-03
        33.404240134            0.68377675770E-02
        17.805918203           -0.10308158997E-01
        4.8990510353           0.22670457658
        2.4516753106           0.44180113937
        1.1820693432           0.36775472225
        I    D
        0.52923110068          1.0000000
        I    D
        0.17000000000          1.0000000
        I    D
        0.61341708807E-01      1.0000000
        I    F
        2.1800000              1.0000000
        I    F
        0.44141808             1.0000000

        * library 6-31+G* except I
        end

        ecp
        I nelec 28
        I ul
        2      19.45860900           -21.84204000
        2      19.34926000           -28.46819100
        2       4.82376700            -0.24371300
        2       4.88431500            -0.32080400
        I S
        2      40.01583500            49.99429300
        2      17.42974700           281.02531700
        2       9.00548400            61.57332600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I P
        2      15.35546600            67.44284100
        2      14.97183300           134.88113700
        2       8.96016400            14.67505100
        2       8.25909600            29.37566600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I D
        2      15.06890800            35.43952900
        2      14.55532200            53.17605700
        2       6.71864700             9.06719500
        2       6.45639300            13.20693700
        2       1.19177900             0.08933500
        2       1.29115700             0.05238000
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_genecp_bs_multiple_atoms(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Ph2I_cation.xyz",
            "software": "nwchem",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=Def2-TZVPD;H=Def2-TZVPD;C=Def2-TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Ph2I_cation
        memory total 1250 mb
        charge 1

        geometry units angstroms noautosym
        C   -3.06870000  -2.28540000   0.00000000
        C   -1.67350000  -2.28540000   0.00000000
        C   -0.97600000  -1.07770000   0.00000000
        C   -1.67360000   0.13090000  -0.00120000
        C   -3.06850000   0.13080000  -0.00170000
        C   -3.76610000  -1.07740000  -0.00070000
        H   -3.61840000  -3.23770000   0.00040000
        H   -1.12400000  -3.23790000   0.00130000
        H    0.12370000  -1.07760000   0.00060000
        H   -1.12340000   1.08300000  -0.00130000
        H   -4.86570000  -1.07720000  -0.00090000
        I   -4.11890000   1.94920000  -0.00350000
        C   -4.64360000   2.85690000  -1.82310000
        C   -3.77180000   3.76300000  -2.42740000
        C   -5.86360000   2.55380000  -2.42750000
        C   -4.12020000   4.36650000  -3.63560000
        H   -2.81040000   4.00240000  -1.95030000
        C   -6.21180000   3.15650000  -3.63650000
        H   -6.55070000   1.83950000  -1.95140000
        C   -5.34050000   4.06290000  -4.24060000
        H   -3.43340000   5.08120000  -4.11170000
        H   -7.17360000   2.91710000  -4.11310000
        H   -5.61500000   4.53870000  -5.19320000
        end

        basis spherical
        #BASIS SET: (12s,11p,9d,2f) -> [7s,6p,4d,2f]
        I    S
        5899.5791533              0.24188269271E-03
        898.54238765             0.15474041742E-02
        200.37237912             0.42836684457E-02
        31.418053840           -0.39417936275E-01
        15.645987838            0.96086691992
        I    S
        11.815741857            0.75961524091
        6.4614458287           0.42495501835
        I    S
        2.3838067579           1.0000000
        I    S
        1.1712089662           1.0000000
        I    S
        0.32115875757          1.0000000
        I    S
        0.12387919364          1.0000000
        I    S
        0.43491550641E-01      1.0000000
        I    P
        197.30030547             0.73951226905E-03
        20.061411349            0.66168450008E-01
        9.7631460485          -0.28554662348
        I    P
        12.984316904           -0.49096186164E-01
        3.6199503008           0.38914432482
        2.0232273090           0.65610817262
        1.0367490559           0.31803551647
        I    P
        0.45937816000          1.0000000
        I    P
        0.19116532928          1.0000000
        I    P
        0.74878813023E-01      1.0000000
        I    P
        0.21653491846E-01      1.0000000
        I    D
        119.12671745             0.82596039573E-03
        33.404240134            0.68377675770E-02
        17.805918203           -0.10308158997E-01
        4.8990510353           0.22670457658
        2.4516753106           0.44180113937
        1.1820693432           0.36775472225
        I    D
        0.52923110068          1.0000000
        I    D
        0.17000000000          1.0000000
        I    D
        0.61341708807E-01      1.0000000
        I    F
        2.1800000              1.0000000
        I    F
        0.44141808             1.0000000

        #BASIS SET: (5s,2p) -> [3s,2p]
        H    S
        34.0613410              0.60251978E-02
        5.1235746              0.45021094E-01
        1.1646626              0.20189726
        H    S
        0.32723041             1.0000000
        H    S
        0.10307241             1.0000000
        H    P
        0.8000000              1.0000000
        H    P
        0.95774129632E-01      1.0000000

        #BASIS SET: (12s,6p,3d,1f) -> [6s,3p,3d,1f]
        C    S
        13575.3496820              0.22245814352E-03
        2035.2333680              0.17232738252E-02
        463.22562359             0.89255715314E-02
        131.20019598             0.35727984502E-01
        42.853015891            0.11076259931
        15.584185766            0.24295627626
        C    S
        6.2067138508           0.41440263448
        2.5764896527           0.23744968655
        C    S
        0.57696339419          1.0000000
        C    S
        0.22972831358          1.0000000
        C    S
        0.95164440028E-01      1.0000000
        C    S
        0.48475401370E-01      1.0000000
        C    P
        34.697232244            0.53333657805E-02
        7.9582622826           0.35864109092E-01
        2.3780826883           0.14215873329
        0.81433208183          0.34270471845
        C    P
        0.28887547253           .46445822433
        C    P
        0.10056823671           .24955789874
        C    D
        1.09700000             1.0000000
        C    D
        0.31800000             1.0000000
        C    D
        0.90985336424E-01      1.0000000
        C    F
        0.76100000             1.0000000

        end

        ecp
        I nelec 28
        I ul
        2      19.45860900           -21.84204000
        2      19.34926000           -28.46819100
        2       4.82376700            -0.24371300
        2       4.88431500            -0.32080400
        I S
        2      40.01583500            49.99429300
        2      17.42974700           281.02531700
        2       9.00548400            61.57332600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I P
        2      15.35546600            67.44284100
        2      14.97183300           134.88113700
        2       8.96016400            14.67505100
        2       8.25909600            29.37566600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I D
        2      15.06890800            35.43952900
        2      14.55532200            53.17605700
        2       6.71864700             9.06719500
        2       6.45639300            13.20693700
        2       1.19177900             0.08933500
        2       1.29115700             0.05238000
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_genecp_bs_multiple_atoms_synonyms(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Ph2I_cation.xyz",
            "software": "nwchem",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=def2-tzvpd;H=def2TZVPD;C=def2 TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Ph2I_cation
        memory total 1250 mb
        charge 1

        geometry units angstroms noautosym
        C   -3.06870000  -2.28540000   0.00000000
        C   -1.67350000  -2.28540000   0.00000000
        C   -0.97600000  -1.07770000   0.00000000
        C   -1.67360000   0.13090000  -0.00120000
        C   -3.06850000   0.13080000  -0.00170000
        C   -3.76610000  -1.07740000  -0.00070000
        H   -3.61840000  -3.23770000   0.00040000
        H   -1.12400000  -3.23790000   0.00130000
        H    0.12370000  -1.07760000   0.00060000
        H   -1.12340000   1.08300000  -0.00130000
        H   -4.86570000  -1.07720000  -0.00090000
        I   -4.11890000   1.94920000  -0.00350000
        C   -4.64360000   2.85690000  -1.82310000
        C   -3.77180000   3.76300000  -2.42740000
        C   -5.86360000   2.55380000  -2.42750000
        C   -4.12020000   4.36650000  -3.63560000
        H   -2.81040000   4.00240000  -1.95030000
        C   -6.21180000   3.15650000  -3.63650000
        H   -6.55070000   1.83950000  -1.95140000
        C   -5.34050000   4.06290000  -4.24060000
        H   -3.43340000   5.08120000  -4.11170000
        H   -7.17360000   2.91710000  -4.11310000
        H   -5.61500000   4.53870000  -5.19320000
        end

        basis spherical
        #BASIS SET: (12s,11p,9d,2f) -> [7s,6p,4d,2f]
        I    S
        5899.5791533              0.24188269271E-03
        898.54238765             0.15474041742E-02
        200.37237912             0.42836684457E-02
        31.418053840           -0.39417936275E-01
        15.645987838            0.96086691992
        I    S
        11.815741857            0.75961524091
        6.4614458287           0.42495501835
        I    S
        2.3838067579           1.0000000
        I    S
        1.1712089662           1.0000000
        I    S
        0.32115875757          1.0000000
        I    S
        0.12387919364          1.0000000
        I    S
        0.43491550641E-01      1.0000000
        I    P
        197.30030547             0.73951226905E-03
        20.061411349            0.66168450008E-01
        9.7631460485          -0.28554662348
        I    P
        12.984316904           -0.49096186164E-01
        3.6199503008           0.38914432482
        2.0232273090           0.65610817262
        1.0367490559           0.31803551647
        I    P
        0.45937816000          1.0000000
        I    P
        0.19116532928          1.0000000
        I    P
        0.74878813023E-01      1.0000000
        I    P
        0.21653491846E-01      1.0000000
        I    D
        119.12671745             0.82596039573E-03
        33.404240134            0.68377675770E-02
        17.805918203           -0.10308158997E-01
        4.8990510353           0.22670457658
        2.4516753106           0.44180113937
        1.1820693432           0.36775472225
        I    D
        0.52923110068          1.0000000
        I    D
        0.17000000000          1.0000000
        I    D
        0.61341708807E-01      1.0000000
        I    F
        2.1800000              1.0000000
        I    F
        0.44141808             1.0000000

        #BASIS SET: (5s,2p) -> [3s,2p]
        H    S
        34.0613410              0.60251978E-02
        5.1235746              0.45021094E-01
        1.1646626              0.20189726
        H    S
        0.32723041             1.0000000
        H    S
        0.10307241             1.0000000
        H    P
        0.8000000              1.0000000
        H    P
        0.95774129632E-01      1.0000000

        #BASIS SET: (12s,6p,3d,1f) -> [6s,3p,3d,1f]
        C    S
        13575.3496820              0.22245814352E-03
        2035.2333680              0.17232738252E-02
        463.22562359             0.89255715314E-02
        131.20019598             0.35727984502E-01
        42.853015891            0.11076259931
        15.584185766            0.24295627626
        C    S
        6.2067138508           0.41440263448
        2.5764896527           0.23744968655
        C    S
        0.57696339419          1.0000000
        C    S
        0.22972831358          1.0000000
        C    S
        0.95164440028E-01      1.0000000
        C    S
        0.48475401370E-01      1.0000000
        C    P
        34.697232244            0.53333657805E-02
        7.9582622826           0.35864109092E-01
        2.3780826883           0.14215873329
        0.81433208183          0.34270471845
        C    P
        0.28887547253           .46445822433
        C    P
        0.10056823671           .24955789874
        C    D
        1.09700000             1.0000000
        C    D
        0.31800000             1.0000000
        C    D
        0.90985336424E-01      1.0000000
        C    F
        0.76100000             1.0000000

        end

        ecp
        I nelec 28
        I ul
        2      19.45860900           -21.84204000
        2      19.34926000           -28.46819100
        2       4.82376700            -0.24371300
        2       4.88431500            -0.32080400
        I S
        2      40.01583500            49.99429300
        2      17.42974700           281.02531700
        2       9.00548400            61.57332600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I P
        2      15.35546600            67.44284100
        2      14.97183300           134.88113700
        2       8.96016400            14.67505100
        2       8.25909600            29.37566600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I D
        2      15.06890800            35.43952900
        2      14.55532200            53.17605700
        2       6.71864700             9.06719500
        2       6.45639300            13.20693700
        2       1.19177900             0.08933500
        2       1.29115700             0.05238000
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_ecp(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "AuI.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=Def2-TZVPD;Au=Def2-TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start AuI
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        Au  -9.27600000  -1.06330000   0.00000000
        I   -6.60600000  -1.06330000   0.00000000
        end

        basis spherical
        #BASIS SET: (12s,11p,9d,2f) -> [7s,6p,4d,2f]
        I    S
        5899.5791533              0.24188269271E-03
        898.54238765             0.15474041742E-02
        200.37237912             0.42836684457E-02
        31.418053840           -0.39417936275E-01
        15.645987838            0.96086691992
        I    S
        11.815741857            0.75961524091
        6.4614458287           0.42495501835
        I    S
        2.3838067579           1.0000000
        I    S
        1.1712089662           1.0000000
        I    S
        0.32115875757          1.0000000
        I    S
        0.12387919364          1.0000000
        I    S
        0.43491550641E-01      1.0000000
        I    P
        197.30030547             0.73951226905E-03
        20.061411349            0.66168450008E-01
        9.7631460485          -0.28554662348
        I    P
        12.984316904           -0.49096186164E-01
        3.6199503008           0.38914432482
        2.0232273090           0.65610817262
        1.0367490559           0.31803551647
        I    P
        0.45937816000          1.0000000
        I    P
        0.19116532928          1.0000000
        I    P
        0.74878813023E-01      1.0000000
        I    P
        0.21653491846E-01      1.0000000
        I    D
        119.12671745             0.82596039573E-03
        33.404240134            0.68377675770E-02
        17.805918203           -0.10308158997E-01
        4.8990510353           0.22670457658
        2.4516753106           0.44180113937
        1.1820693432           0.36775472225
        I    D
        0.52923110068          1.0000000
        I    D
        0.17000000000          1.0000000
        I    D
        0.61341708807E-01      1.0000000
        I    F
        2.1800000              1.0000000
        I    F
        0.44141808             1.0000000

        #BASIS SET: (8s,8p,6d,1f) -> [6s,5p,3d,1f]
        Au    S
        30.000000000            0.20749231108
        27.000000000           -0.33267893394
        14.746824331            0.38302817958
        Au    S
        5.6017248938           1.0000000
        Au    S
        1.3874162443           1.0000000
        Au    S
        0.62923031957          1.0000000
        Au    S
        0.14027517613          1.0000000
        Au    S
        0.49379413761E-01      1.0000000
        Au    P
        15.500000000            0.15001711880
        14.000000000           -0.23609813183
        6.4227368205           0.31458896948
        1.6595601681          -0.57279670446
        Au    P
        0.79402913993          1.0000000
        Au    P
        0.35125155397          1.0000000
        Au    P
        0.11801737494          1.0000000
        Au    P
        0.45000000000E-01      1.0000000
        Au    D
        9.5524098656           0.40145559502E-01
        7.2698886937          -0.93690906606E-01
        1.7746496789           0.31746282317
        0.79960541055          0.46795192483
        Au    D
        0.33252279372          1.0000000
        Au    D
        0.12445133105          1.0000000
        Au    F
        0.7248200              1.0000000

        end

        ecp
        I nelec 28
        I ul
        2      19.45860900           -21.84204000
        2      19.34926000           -28.46819100
        2       4.82376700            -0.24371300
        2       4.88431500            -0.32080400
        I S
        2      40.01583500            49.99429300
        2      17.42974700           281.02531700
        2       9.00548400            61.57332600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I P
        2      15.35546600            67.44284100
        2      14.97183300           134.88113700
        2       8.96016400            14.67505100
        2       8.25909600            29.37566600
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400
        I D
        2      15.06890800            35.43952900
        2      14.55532200            53.17605700
        2       6.71864700             9.06719500
        2       6.45639300            13.20693700
        2       1.19177900             0.08933500
        2       1.29115700             0.05238000
        2      19.45860900            21.84204000
        2      19.34926000            28.46819100
        2       4.82376700             0.24371300
        2       4.88431500             0.32080400

        Au nelec 60
        Au ul
        2       4.78982000            30.49008890
        2       2.39491000             5.17107381
        Au S
        2      13.20510000           426.84667920
        2       6.60255000            37.00708285
        2       4.78982000           -30.49008890
        2       2.39491000            -5.17107381
        Au P
        2      10.45202000           261.19958038
        2       5.22601000            26.96249604
        2       4.78982000           -30.49008890
        2       2.39491000            -5.17107381
        Au D
        2       7.85110000           124.79066561
        2       3.92555000            16.30072573
        2       4.78982000           -30.49008890
        2       2.39491000            -5.17107381
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_builtin_custom_bs(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "O=Z3Pol;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis spherical

        * library 6-31+G* except O
        O library Z3Pol
        end

        dft
        xc b3lyp
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_global_specification(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "SCF(Tight)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_confirmed_command_specification(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "opt(maxiter 5)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        driver
        maxiter 5
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_global_specification(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "SCF(Tight);SCF(Direct)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        direct
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_global_specification_alternative(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "SCF(Tight,Direct)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        direct
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_global_specification_alternative2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "SCF(Tight, Direct)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        direct
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_global_specification_alternative3(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "SCF(Tight)   scf(direct)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        direct
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_global_specification2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "0",
            "specifications": "SCF(Tight);opt(maxiter 5)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        end

        driver
        maxiter 5
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_alternative(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "0",
            "specifications": "SCF(Tight) opt(maxiter 5)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        end

        driver
        maxiter 5
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_alternative2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "0",
            "specifications": "SCF(Tight); opt(maxiter 5) scf(direct)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        direct
        end

        driver
        maxiter 5
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_mixed(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "optfreq",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "0",
            "specifications": "opt(maxiter 5); scf(Tight); freq(temp 298)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        tight
        end

        driver
        maxiter 5
        end

        freq
        temp 298
        end

        task dft optimize
        task dft freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_neb(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "mep",
            "name": "neb_calculation",
            "file": "elimination_substrate.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "neb(xyz_path xyz_path.xyz)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start neb_calculation
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        C   -0.74277000   0.14309000   0.12635000
        C    0.71308000  -0.12855000  -0.16358000
        Cl   0.90703000  -0.47793000  -1.61303000
        H   -0.84928000   0.38704000   1.20767000
        H   -1.36298000  -0.72675000  -0.06978000
        H   -1.11617000   0.99405000  -0.43583000
        H    1.06397000  -0.95639000   0.44985000
        H    1.30839000   0.75217000   0.07028000
        O   -0.91651000   0.74066000   3.00993000
        H   -1.82448000   0.94856000   3.28105000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        neb
        xyz_path xyz_path.xyz
        end

        task dft neb ignore
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_neb_synonym(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "mep",
            "name": "neb_calculation",
            "file": "elimination_substrate.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "mep(xyz_path xyz_path.xyz)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start neb_calculation
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        C   -0.74277000   0.14309000   0.12635000
        C    0.71308000  -0.12855000  -0.16358000
        Cl   0.90703000  -0.47793000  -1.61303000
        H   -0.84928000   0.38704000   1.20767000
        H   -1.36298000  -0.72675000  -0.06978000
        H   -1.11617000   0.99405000  -0.43583000
        H    1.06397000  -0.95639000   0.44985000
        H    1.30839000   0.75217000   0.07028000
        O   -0.91651000   0.74066000   3.00993000
        H   -1.82448000   0.94856000   3.28105000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        neb
        xyz_path xyz_path.xyz
        end

        task dft neb ignore
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_string_method(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "mep",
            "name": "neb_calculation",
            "file": "elimination_substrate.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "neb(xyz_path xyz_path.xyz);string;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start neb_calculation
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        C   -0.74277000   0.14309000   0.12635000
        C    0.71308000  -0.12855000  -0.16358000
        Cl   0.90703000  -0.47793000  -1.61303000
        H   -0.84928000   0.38704000   1.20767000
        H   -1.36298000  -0.72675000  -0.06978000
        H   -1.11617000   0.99405000  -0.43583000
        H    1.06397000  -0.95639000   0.44985000
        H    1.30839000   0.75217000   0.07028000
        O   -0.91651000   0.74066000   3.00993000
        H   -1.82448000   0.94856000   3.28105000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        string
        xyz_path xyz_path.xyz
        end

        task dft string ignore
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_special_char(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "!#",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library Def2-SVP
        end

        dft
        xc m06-2x
        mult 1
        end

        task dft optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_mem_m(self):
        params = {
            "nproc": 8,
            "mem": "10000 m",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_mem_GB(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_exchange_correlation_functional(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "xm05 cm05-2x",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc xm05 cm05-2x
        mult 1
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_exchange_correlation_functional2(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem ",
            "method": "becke88 perdew91",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc becke88 perdew91
        mult 1
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_restricted_HF(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "rHF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        rhf
        singlet
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_unrestricted_HF(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "uHF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        uhf
        singlet
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_unrestricted_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M062x",
            "basis_set": "3-21G",
            "charge": "-1",
            "specifications": "dft(odft)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc m06-2x
        mult 1
        odft
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_restricted_DFT(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M062x",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc m06-2x
        mult 1
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_freq(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "opt+freq",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "hf",
            "charge": "-1",
            "basis_set": "3-21G",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        task scf optimize
        task scf freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_freq2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "optfreq",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "hf",
            "charge": "-1",
            "basis_set": "3-21G",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        task scf optimize
        task scf freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_freq_spec(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "opt+freq",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "HF",
            "charge": "-1",
            "specifications": "freq(animate)",
            "basis_set": "3-21G",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        freq
        animate
        end

        task scf optimize
        task scf freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_d3(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "M062X",
            "basis_set": "3-21G",
            "charge": "-1",
            "d3": True,
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc m06-2x
        mult 1
        disp vdw 3
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_d3bj(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "PBE0",
            "basis_set": "3-21G",
            "charge": "-1",
            "d3bj": True,
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc pbe0
        mult 1
        disp vdw 4
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_d3_d3bj_crash(self):
        params = {
            "nproc": 8,
            "mem": "10GB",
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "PBE0",
            "basis_set": "3-21G",
            "charge": "-1",
            "d3": True,
            "d3bj": True,
        }

        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_smd_custom_radius(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "nwchem",
            "method": "B3LYP",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "custom_solvation_radii": "Cl=1.00",
            "solvation_model": "SMD",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start I
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        I    0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc b3lyp
        mult 1
        end

        cosmo
        minbem 3
        ificos 1
        solvent chcl3
        do_cosmo_smd
        parameters I_sol.parameters
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.radii_parameters.strip(), "Cl 1.0")

    def test_smd_custom_radii(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "nwchem",
            "method": "B3LYP",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "custom_solvation_radii": "Cl=1.00;I=1.55;",
            "solvation_model": "SMD",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start I
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        I    0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        dft
        xc b3lyp
        mult 1
        end

        cosmo
        minbem 3
        ificos 1
        solvent chcl3
        do_cosmo_smd
        parameters I_sol.parameters
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.radii_parameters.strip(), "Cl 1.0\nI 1.55")

    def test_pcm_custom_radius(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "nwchem",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_radii": "default",
            "custom_solvation_radii": "Cl=1.00",
            "solvation_model": "COSMO",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start I
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        I    0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        end

        cosmo
        minbem 3
        ificos 1
        solvent chcl3
        parameters I_sol.parameters
        end

        task scf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.radii_parameters.strip(), "Cl 1.0")

    def test_unavailable_calc_type(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "tda",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "b3lyp",
            "charge": "-1",
            "basis_set": "3-21G",
        }

        with self.assertRaises(ImpossibleCalculation):
            self.generate_calculation(**params)

    def test_sp_MP2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "mp2",
            "basis_set": "cc-pVDZ",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end

        
        task mp2 energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_MP2_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "mp2",
            "basis_set": "cc-pVDZ",
            "specifications": "mp2(freeze atomic)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end

        mp2
        freeze atomic
        end
        
        task mp2 energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_MP2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "mp2",
            "basis_set": "cc-pVDZ",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end


        task mp2 freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_MP2_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Opt",
            "file": "ethanol.xyz",
            "software": "nwchem",
            "charge": "0",
            "method": "mp2",
            "basis_set": "aug-cc-pVDZ",
            "specifications": "opt(maxiter 20);mp2(tight)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start ethanol
        memory total 1250 mb
        charge 0

        geometry units angstroms noautosym
        C   -1.31970000  -0.64380000   0.00000000
        H   -0.96310000  -1.65260000   0.00000000
        H   -0.96310000  -0.13940000  -0.87370000
        H   -2.38970000  -0.64380000   0.00000000
        C   -0.80640000   0.08220000   1.25740000
        H   -1.16150000   1.09160000   1.25640000
        H   -1.16470000  -0.42110000   2.13110000
        O    0.62360000   0.07990000   1.25870000
        H    0.94410000   0.53240000   2.04240000
        end

        basis
        * library aug-cc-pVDZ
        end

        mp2
        tight
        end

        driver
        maxiter 20
        end

        task mp2 optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_freq_multiple_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "opt+freq",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "hf",
            "charge": "-1",
            "basis_set": "3-21G",
            "specifications": "scf(maxiter 30);opt(maxiter 20);freq(temp 373)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library 3-21G
        end

        scf
        singlet
        maxiter 30
        end

        driver
        maxiter 20
        end

        freq
        temp 373
        end

        task scf optimize
        task scf freq
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_CCSD(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "ccsd",
            "basis_set": "cc-pVDZ",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end

        
        task ccsd energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_CCSDT(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "CCSD(T)",
            "basis_set": "cc-pVDZ",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end

        
        task ccsd(t) energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_CCSD_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "ccsd",
            "basis_set": "cc-pVDZ",
            "specifications": "cc(maxiter 14)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end

        ccsd
        maxiter 14
        end
        
        task ccsd energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_CCSDT_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Single Point Energy",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "CCSD(T)",
            "basis_set": "cc-pVDZ",
            "specifications": "cc(maxiter 14)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end

        ccsd
        maxiter 14
        end
        
        task ccsd(t) energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_CCSDT_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "CCSD(T)",
            "basis_set": "cc-pVDZ",
            "specifications": "cc(maxiter 14);opt(tight)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVDZ
        end

        ccsd
        maxiter 14
        end

        driver
        tight
        end

        task ccsd(t) optimize
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_RIMP2(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "sp",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "RI-MP2",
            "basis_set": "cc-pVTZ",
            "specifications": "mp2(freeze virtual 5)",
            "density_fitting": "cc-pVDZ",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVTZ
        end

        basis "cd basis"
        * library cc-pVDZ
        end

        mp2
        freeze virtual 5
        end

        task rimp2 energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_density_fitting(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "sp",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "PBE",
            "basis_set": "cc-pVTZ",
            "density_fitting": "cc-pVDZ",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVTZ
        end

        basis "cd basis"
        * library cc-pVDZ
        end

        dft
        xc PBE
        mult 1
        adft
        end

        task dft energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_CASSCF(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "sp",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "mcscf",
            "basis_set": "cc-pVTZ",
            "specifications": "mcscf(active 5);mcscf(actelec 2);",
        }

        inp = self.generate_calculation(**params)

        REF = """
        TITLE "File created by ccinput"
        start Cl
        memory total 1250 mb
        charge -1

        geometry units angstroms noautosym
        Cl   0.00000000   0.00000000   0.00000000
        end

        basis
        * library cc-pVTZ
        end

        mcscf
        multiplicity 1
        active 5
        actelec 2
        end

        task mcscf energy
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_CASSCF_no_specifications(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "sp",
            "file": "Cl.xyz",
            "software": "nwchem",
            "charge": "-1",
            "method": "mcscf",
            "basis_set": "cc-pVTZ",
        }

        with self.assertRaises(MissingParameter):
            self.generate_calculation(**params)

from ccinput.tests.testing_utilities import InputTests
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation, UnimplementedError


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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
            memory total 10000 mb
            charge -1

            geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge 0

        geometry units angstroms
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

        constraints
        spring bond 1 2 200.0 2.02195471
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
        with self.assertRaises(UnimplementedError):
            self.generate_calculation(**params)

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
        memory total 10000 mb
        charge 0

        geometry units angstroms
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

        constraints
        spring dihedral 4 1 5 8 200.0 -179.88772953
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
        memory total 10000 mb
        charge 0

        geometry units angstroms
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

        constraints
        spring dihedral 4 1 5 8 200.0 -179.88772953
        spring dihedral 1 2 3 4 200.0 35.26896033
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge 0

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge 0

        geometry units angstroms
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
        memory total 10000 mb
        charge 0

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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
        memory total 10000 mb
        charge -1

        geometry units angstroms
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


    def test_unavailable_calc_type(self):
        params = {
            "nproc": 8,
            "mem": "10000MB",
            "type": "Minimum Energy Path",
            "file": "Cl.xyz",
            "software": "nwchem",
            "method": "b3lyp",
            "charge": "-1",
            "basis_set": "3-21G",
        }

        with self.assertRaises(ImpossibleCalculation):
            self.generate_calculation(**params)


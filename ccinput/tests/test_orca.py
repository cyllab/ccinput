from ccinput.tests.testing_utilities import InputTests
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation


class OrcaTests(InputTests):
    def test_sp_SE(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "AM1",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP AM1
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 1000
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD_lowercase(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "smd",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD18(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
            "solvation_radii": "SMD18",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        I 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        radius[53] 2.74
        radius[35] 2.60
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD18_lowercase(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "smd",
            "solvation_radii": "smd18",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        I 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        radius[53] 2.74
        radius[35] 2.60
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvation_octanol_smd(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "octanol",
            "solvation_model": "SMD",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "1-octanol"
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_CPCM(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "CPCM",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(chloroform)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvent_synonym(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "CHCL3",
            "solvation_model": "SMD",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvation_octanol_cpcm1(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "octanol",
            "solvation_model": "CPCM",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(octanol)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvation_octanol_cpcm2(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "1-octanol",
            "solvation_model": "CPCM",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(octanol)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_solvation(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
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
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_specifications(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "TightSCF",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP tightscf
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_multiple_specifications(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "TightSCF GRID6",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP tightscf grid6
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_duplicate_specifications(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "tightscf TightSCF GRID6",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP tightscf grid6
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_MP2(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "RI-MP2",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "cc-pvtz/C",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP RI-MP2 cc-pVTZ cc-pVTZ/C
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_MP2_default_aux(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "RI-MP2",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP RI-MP2 cc-pVTZ cc-pVTZ/C
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_MP2_override_aux(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "RI-MP2",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "cc-pVDZ/C",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP RI-MP2 cc-pVTZ cc-pVDZ/C
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_MP2_without_RI(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "MP2",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP MP2 cc-pVTZ
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DLPNO_CCSDT(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "DLPNO-CCSD(T)",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "cc-pVTZ/C",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DLPNO_CCSDT_default_aux(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "DLPNO-CCSD(T)",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DLPNO_CCSDT_default_aux_specifications(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "DLPNO-CCSD(T)",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "tightscf",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C tightscf
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DLPNO_CCSDT_default_aux_solvation(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "DLPNO-CCSD(T)",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "solvent": "ch2cl2",
            "solvation_model": "CPCM",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C CPCM(ch2cl2)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_CCSDT(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "CCSD(T)",
            "basis_set": "cc-pVTZ",
            "charge": "-1",
            "specifications": "",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP CCSD(T) cc-pVTZ
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_SE(self):
        params = {
            "nproc": 8,
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "AM1",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT AM1
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 1000
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_HF(self):
        params = {
            "nproc": 8,
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT(self):
        params = {
            "nproc": 8,
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "ORCA",
            "charge": "-1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_SE(self):
        params = {
            "nproc": 8,
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "AM1",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !FREQ AM1
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 1000
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_HF(self):
        params = {
            "nproc": 8,
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !FREQ HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_DFT(self):
        params = {
            "nproc": 8,
            "type": "Frequency Calculation",
            "file": "Cl.xyz",
            "software": "ORCA",
            "charge": "-1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !FREQ B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    # opt mod SE and HF

    def test_scan_bond_DFT(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Scan_9_1.4_10/1_2;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240
        *
        %geom Scan
        B 0 1 = 9.00, 1.40, 10
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_bond_DFT_auto(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Scan_auto_1.4_10/1_2;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240
        *
        %geom Scan
        B 0 1 = 1.07, 1.40, 10
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_opt_mod(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "",
        }

        with self.assertRaises(Exception):
            self.generate_calculation(**params)

    def test_no_method(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "",
            "basis_set": "6-31+G(d,p)",
            "constraints": "",
        }

        with self.assertRaises(Exception):
            self.generate_calculation(**params)

    def test_scan_angle_DFT(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Scan_9_90_10/2_1_3;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240
        *
        %geom Scan
        A 1 0 2 = 9.00, 90.00, 10
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_dihedral_DFT(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Scan_9_1_10/4_1_5_8;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240
        *
        %geom Scan
        D 3 0 4 7 = 9.00, 1.00, 10
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_no_constraint(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_freeze_bond_DFT(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Freeze/1_2;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240
        *
        %geom Constraints
        { B 0 1 C }
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_angle_DFT(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Freeze/2_1_3;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240
        *
        %geom Constraints
        { A 1 0 2 C }
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_dihedral_DFT(self):
        params = {
            "nproc": 8,
            "type": "Constrained Optimisation",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "constraints": "Freeze/4_1_5_8;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240
        *
        %geom Constraints
        { D 3 0 4 7 C }
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nmr_DFT(self):
        params = {
            "nproc": 8,
            "type": "NMR Prediction",
            "file": "Cl.xyz",
            "software": "ORCA",
            "charge": "-1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !NMR B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_irrelevant_gen_bs(self):
        params = {
            "nproc": 8,
            "type": "NMR Prediction",
            "file": "Cl.xyz",
            "software": "ORCA",
            "charge": "-1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "N=Def2-SVP;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !NMR B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts_DFT(self):
        params = {
            "nproc": 8,
            "type": "TS Optimisation",
            "file": "mini_ts.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPTTS B3LYP 6-31+G(d,p)
        *xyz 0 1
        N   1.08764072053386     -0.33994563112543     -0.00972525479568
        H   1.99826836912112      0.05502842705407      0.00651240826058
        H   0.59453997172323     -0.48560162159600      0.83949232123172
        H   0.66998093862168     -0.58930117433261     -0.87511947469677
        *
        %geom
        Calc_Hess true
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts_xtb(self):
        params = {
            "nproc": 8,
            "type": "TS Optimisation",
            "file": "mini_ts.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "xtb",
        }

        inp = self.generate_calculation(**params)

        # One cannot calculate the Hessian for use in a TS optimization
        # when using xtb as QM engine.
        REF = """
        !OPTTS xtb
        *xyz 0 1
        N   1.08764072053386     -0.33994563112543     -0.00972525479568
        H   1.99826836912112      0.05502842705407      0.00651240826058
        H   0.59453997172323     -0.48560162159600      0.83949232123172
        H   0.66998093862168     -0.58930117433261     -0.87511947469677
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

        # combination tests

    def test_ts_DFT_custom_bs(self):
        params = {
            "nproc": 8,
            "type": "TS Optimisation",
            "file": "mini_ts.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "N=Def2-SVP;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPTTS B3LYP 6-31+G(d,p)
        *xyz 0 1
        N   1.08764072053386     -0.33994563112543     -0.00972525479568
        H   1.99826836912112      0.05502842705407      0.00651240826058
        H   0.59453997172323     -0.48560162159600      0.83949232123172
        H   0.66998093862168     -0.58930117433261     -0.87511947469677
        *
        %geom
        Calc_Hess true
        end
        %basis
        NewGTO N "Def2-SVP" end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts_DFT_custom_bs_synonym(self):
        params = {
            "nproc": 8,
            "type": "TS Optimisation",
            "file": "mini_ts.xyz",
            "software": "ORCA",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "N=def2Svp;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPTTS B3LYP 6-31+G(d,p)
        *xyz 0 1
        N   1.08764072053386     -0.33994563112543     -0.00972525479568
        H   1.99826836912112      0.05502842705407      0.00651240826058
        H   0.59453997172323     -0.48560162159600      0.83949232123172
        H   0.66998093862168     -0.58930117433261     -0.87511947469677
        *
        %geom
        Calc_Hess true
        end
        %basis
        NewGTO N "Def2-SVP" end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT_custom_bs_ecp(self):
        params = {
            "nproc": 8,
            "type": "Geometrical Optimisation",
            "file": "Ph2I_cation.xyz",
            "software": "ORCA",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=Def2-TZVPD;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 1 1
        C         -3.06870       -2.28540        0.00000
        C         -1.67350       -2.28540        0.00000
        C         -0.97600       -1.07770        0.00000
        C         -1.67360        0.13090       -0.00120
        C         -3.06850        0.13080       -0.00170
        C         -3.76610       -1.07740       -0.00070
        H         -3.61840       -3.23770        0.00040
        H         -1.12400       -3.23790        0.00130
        H          0.12370       -1.07760        0.00060
        H         -1.12340        1.08300       -0.00130
        H         -4.86570       -1.07720       -0.00090
        I         -4.11890        1.94920       -0.00350
        C         -4.64360        2.85690       -1.82310
        C         -3.77180        3.76300       -2.42740
        C         -5.86360        2.55380       -2.42750
        C         -4.12020        4.36650       -3.63560
        H         -2.81040        4.00240       -1.95030
        C         -6.21180        3.15650       -3.63650
        H         -6.55070        1.83950       -1.95140
        C         -5.34050        4.06290       -4.24060
        H         -3.43340        5.08120       -4.11170
        H         -7.17360        2.91710       -4.11310
        H         -5.61500        4.53870       -5.19320
        *
        %basis
        NewGTO I "Def2-TZVPD" end
        NewECP I "def2-ECP" end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT_two_custom_bs_ecp(self):
        params = {
            "nproc": 8,
            "type": "Geometrical Optimisation",
            "file": "Ph2I_cation.xyz",
            "software": "ORCA",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=Def2-TZVPD;H=Def2-TZVP",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 1 1
        C         -3.06870       -2.28540        0.00000
        C         -1.67350       -2.28540        0.00000
        C         -0.97600       -1.07770        0.00000
        C         -1.67360        0.13090       -0.00120
        C         -3.06850        0.13080       -0.00170
        C         -3.76610       -1.07740       -0.00070
        H         -3.61840       -3.23770        0.00040
        H         -1.12400       -3.23790        0.00130
        H          0.12370       -1.07760        0.00060
        H         -1.12340        1.08300       -0.00130
        H         -4.86570       -1.07720       -0.00090
        I         -4.11890        1.94920       -0.00350
        C         -4.64360        2.85690       -1.82310
        C         -3.77180        3.76300       -2.42740
        C         -5.86360        2.55380       -2.42750
        C         -4.12020        4.36650       -3.63560
        H         -2.81040        4.00240       -1.95030
        C         -6.21180        3.15650       -3.63650
        H         -6.55070        1.83950       -1.95140
        C         -5.34050        4.06290       -4.24060
        H         -3.43340        5.08120       -4.11170
        H         -7.17360        2.91710       -4.11310
        H         -5.61500        4.53870       -5.19320
        *
        %basis
        NewGTO I "Def2-TZVPD" end
        NewECP I "def2-ECP" end
        NewGTO H "Def2-TZVP" end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT_custom_bs_ecp_synonym(self):
        params = {
            "nproc": 8,
            "type": "Geometrical Optimisation",
            "file": "Ph2I_cation.xyz",
            "software": "ORCA",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=Def2 tzvpd;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 1 1
        C         -3.06870       -2.28540        0.00000
        C         -1.67350       -2.28540        0.00000
        C         -0.97600       -1.07770        0.00000
        C         -1.67360        0.13090       -0.00120
        C         -3.06850        0.13080       -0.00170
        C         -3.76610       -1.07740       -0.00070
        H         -3.61840       -3.23770        0.00040
        H         -1.12400       -3.23790        0.00130
        H          0.12370       -1.07760        0.00060
        H         -1.12340        1.08300       -0.00130
        H         -4.86570       -1.07720       -0.00090
        I         -4.11890        1.94920       -0.00350
        C         -4.64360        2.85690       -1.82310
        C         -3.77180        3.76300       -2.42740
        C         -5.86360        2.55380       -2.42750
        C         -4.12020        4.36650       -3.63560
        H         -2.81040        4.00240       -1.95030
        C         -6.21180        3.15650       -3.63650
        H         -6.55070        1.83950       -1.95140
        C         -5.34050        4.06290       -4.24060
        H         -3.43340        5.08120       -4.11170
        H         -7.17360        2.91710       -4.11310
        H         -5.61500        4.53870       -5.19320
        *
        %basis
        NewGTO I "Def2-TZVPD" end
        NewECP I "def2-ECP" end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT_custom_bs_explicit(self):
        params = {
            "nproc": 8,
            "type": "Geometrical Optimisation",
            "file": "Ph2I_cation.xyz",
            "software": "ORCA",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "I=psbkjc;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz 1 1
        C         -3.06870       -2.28540        0.00000
        C         -1.67350       -2.28540        0.00000
        C         -0.97600       -1.07770        0.00000
        C         -1.67360        0.13090       -0.00120
        C         -3.06850        0.13080       -0.00170
        C         -3.76610       -1.07740       -0.00070
        H         -3.61840       -3.23770        0.00040
        H         -1.12400       -3.23790        0.00130
        H          0.12370       -1.07760        0.00060
        H         -1.12340        1.08300       -0.00130
        H         -4.86570       -1.07720       -0.00090
        I         -4.11890        1.94920       -0.00350
        C         -4.64360        2.85690       -1.82310
        C         -3.77180        3.76300       -2.42740
        C         -5.86360        2.55380       -2.42750
        C         -4.12020        4.36650       -3.63560
        H         -2.81040        4.00240       -1.95030
        C         -6.21180        3.15650       -3.63650
        H         -6.55070        1.83950       -1.95140
        C         -5.34050        4.06290       -4.24060
        H         -3.43340        5.08120       -4.11170
        H         -7.17360        2.91710       -4.11310
        H         -5.61500        4.53870       -5.19320
        *
        %basis
        newgto I
        S   4
        1         2.6250000              0.07366000
        2         1.0140000             -0.83687000
        3         0.5009000              0.65624700
        4         0.2023000              0.90074400
        S   2
        1         0.2023000              1.00000000
        2         0.0780000              1.00000000
        S   2
        1         0.0780000              1.00000000
        2         0.0300741              1.00000000
        P   4
        1         2.6250000             -0.00888000
        2         1.0140000             -0.25735100
        3         0.5009000              0.45536800
        4         0.2023000              0.76010700
        P   2
        1         0.2023000              1.00000000
        2         0.0780000              1.00000000
        P   2
        1         0.0780000              1.00000000
        2         0.0300741              1.00000000
        D   2
        1         0.5009000              0.36111300
        2         0.2023000              0.93252210
        D   2
        1         0.0780000              0.93733150
        2         0.0300741              0.34843880
        end

        NewECP I
        N_core 46
        lmax f
        s 3
        1      4.1307100     12.1112300 0
        2      1.3337500    -41.0920600 2
        3      1.4912100     70.7376100 2
        p 3
        1      3.0469200     10.5927100 0
        2      1.0633900    -46.0227300 2
        3      1.1440500     65.0504700 2
        d 2
        1      3.9306300      9.7308900 0
        2      1.0692000     13.9888000 2
        f 2
        1      0.9762800     -3.6963900 1
        2      4.3334300    -14.0030500 1
        end
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_NEB(self):
        params = {
            "nproc": 8,
            "type": "Minimum Energy Path",
            "file": "elimination_substrate.xyz",
            "auxiliary_file": "elimination_product.xyz",
            "software": "orca",
            "charge": -1,
            "method": "gfn2-xtb",
        }

        inp = self.generate_calculation(**params)

        REF = """!NEB xtb2
        *xyz -1 1
        C         -0.74277        0.14309        0.12635
        C          0.71308       -0.12855       -0.16358
        Cl         0.90703       -0.47793       -1.61303
        H         -0.84928        0.38704        1.20767
        H         -1.36298       -0.72675       -0.06978
        H         -1.11617        0.99405       -0.43583
        H          1.06397       -0.95639        0.44985
        H          1.30839        0.75217        0.07028
        O         -0.91651        0.74066        3.00993
        H         -1.82448        0.94856        3.28105
        *
        %neb
        product "calc2.xyz"
        nimages 8
        end
        %MaxCore 125
        %pal
        nprocs 8
        end

        """
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_NEB2(self):
        params = {
            "nproc": 8,
            "type": "Minimum Energy Path",
            "file": "elimination_substrate.xyz",
            "auxiliary_file": "elimination_product.xyz",
            "software": "ORCA",
            "specifications": "--nimages 12",
            "charge": -1,
            "method": "gfn2-xtb",
        }

        inp = self.generate_calculation(**params)

        REF = """!NEB xtb2
        *xyz -1 1
        C         -0.74277        0.14309        0.12635
        C          0.71308       -0.12855       -0.16358
        Cl         0.90703       -0.47793       -1.61303
        H         -0.84928        0.38704        1.20767
        H         -1.36298       -0.72675       -0.06978
        H         -1.11617        0.99405       -0.43583
        H          1.06397       -0.95639        0.44985
        H          1.30839        0.75217        0.07028
        O         -0.91651        0.74066        3.00993
        H         -1.82448        0.94856        3.28105
        *
        %neb
        product "calc2.xyz"
        nimages 12
        end
        %MaxCore 125
        %pal
        nprocs 8
        end

        """
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_NEB_aux_name(self):
        params = {
            "nproc": 8,
            "type": "Minimum Energy Path",
            "file": "elimination_substrate.xyz",
            "auxiliary_file": "elimination_product.xyz",
            "software": "ORCA",
            "specifications": "--nimages 12",
            "charge": -1,
            "method": "gfn2-xtb",
            "aux_name": "product",
        }

        inp = self.generate_calculation(**params)

        REF = """!NEB xtb2
        *xyz -1 1
        C         -0.74277        0.14309        0.12635
        C          0.71308       -0.12855       -0.16358
        Cl         0.90703       -0.47793       -1.61303
        H         -0.84928        0.38704        1.20767
        H         -1.36298       -0.72675       -0.06978
        H         -1.11617        0.99405       -0.43583
        H          1.06397       -0.95639        0.44985
        H          1.30839        0.75217        0.07028
        O         -0.91651        0.74066        3.00993
        H         -1.82448        0.94856        3.28105
        *
        %neb
        product "product.xyz"
        nimages 12
        end
        %MaxCore 125
        %pal
        nprocs 8
        end

        """
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_hirshfeld_pop(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "--phirshfeld",
        }

        inp = self.generate_calculation(**params)
        REF = """
        !SP M062X Def2-SVP
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %output
        Print[ P_Hirshfeld] 1
        end
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_mo(self):
        params = {
            "nproc": 8,
            "type": "MO Calculation",
            "file": "Ph2I_cation.xyz",
            "software": "ORCA",
            "charge": "+1",
            "method": "B3LYP",
            "basis_set": "def2tzvp",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP B3LYP Def2-TZVP
        *xyz 1 1
        C         -3.06870       -2.28540        0.00000
        C         -1.67350       -2.28540        0.00000
        C         -0.97600       -1.07770        0.00000
        C         -1.67360        0.13090       -0.00120
        C         -3.06850        0.13080       -0.00170
        C         -3.76610       -1.07740       -0.00070
        H         -3.61840       -3.23770        0.00040
        H         -1.12400       -3.23790        0.00130
        H          0.12370       -1.07760        0.00060
        H         -1.12340        1.08300       -0.00130
        H         -4.86570       -1.07720       -0.00090
        I         -4.11890        1.94920       -0.00350
        C         -4.64360        2.85690       -1.82310
        C         -3.77180        3.76300       -2.42740
        C         -5.86360        2.55380       -2.42750
        C         -4.12020        4.36650       -3.63560
        H         -2.81040        4.00240       -1.95030
        C         -6.21180        3.15650       -3.63650
        H         -6.55070        1.83950       -1.95140
        C         -5.34050        4.06290       -4.24060
        H         -3.43340        5.08120       -4.11170
        H         -7.17360        2.91710       -4.11310
        H         -5.61500        4.53870       -5.19320
        *
        %plots
        dim1 45
        dim2 45
        dim3 45
        min1 0
        max1 0
        min2 0
        max2 0
        min3 0
        max3 0
        Format Gaussian_Cube
        MO("in-HOMO.cube",66,0);
        MO("in-LUMO.cube",67,0);
        MO("in-LUMOA.cube",68,0);
        MO("in-LUMOB.cube",69,0);
        end

        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nproc(self):
        params = {
            "nproc": 1,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 1000
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_freq(self):
        params = {
            "nproc": 8,
            "type": "opt+freq",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "AM1",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT FREQ AM1
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 1000
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_d3(self):
        params = {
            "nproc": 8,
            "type": "opt",
            "file": "Cl.xyz",
            "software": "ORCA",
            "basis_set": "Def2TZVP",
            "method": "M062X",
            "charge": "-1",
            "d3": True,
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT M062X Def2-TZVP d3zero
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_d3bj(self):
        params = {
            "nproc": 8,
            "type": "opt",
            "file": "Cl.xyz",
            "software": "ORCA",
            "basis_set": "Def2TZVP",
            "method": "PBE0",
            "charge": "-1",
            "d3bj": True,
        }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT PBE0 Def2-TZVP d3bj
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_d3_d3bj_crash(self):
        params = {
            "nproc": 8,
            "type": "opt",
            "file": "Cl.xyz",
            "software": "ORCA",
            "basis_set": "Def2TZVP",
            "method": "PBE0",
            "charge": "-1",
            "d3": True,
            "d3bj": True,
        }

        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_SMD_custom_radius(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
            "custom_solvation_radii": "Cl=1.00;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        radius[17] 1.00
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_SMD_custom_radii(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
            "custom_solvation_radii": "Cl=1.00;Br=2.00;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        radius[17] 1.00
        radius[35] 2.00
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_SMD_custom_radius_and_SMD18(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
            "solvation_radii": "SMD18",
            "custom_solvation_radii": "Cl=1.00;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        radius[17] 1.00
        radius[53] 2.74
        radius[35] 2.60
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_SMD_custom_radius_and_SMD18_clash(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
            "solvation_radii": "SMD18",
            "custom_solvation_radii": "I=3.00;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        smd true
        SMDsolvent "chloroform"
        radius[53] 3.00
        radius[35] 2.60
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_CPCM_custom_radius(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "CPCM",
            "custom_solvation_radii": "Cl=1.00;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(chloroform)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        radius[17] 1.00
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_CPCM_custom_radii(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "CPCM",
            "custom_solvation_radii": "Cl=1.00;Br=2.00;",
        }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(chloroform)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %MaxCore 125
        %pal
        nprocs 8
        end
        %cpcm
        radius[17] 1.00
        radius[35] 2.00
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_unavailable_calc_type(self):
        params = {
            "nproc": 8,
            "type": "Conformational Search",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        with self.assertRaises(ImpossibleCalculation):
            self.generate_calculation(**params)

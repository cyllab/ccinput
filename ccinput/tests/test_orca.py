from ccinput.tests.testing_utilities import InputTests
from ccinput.packages.orca import OrcaCalculation
from ccinput.wrapper import generate_calculation

class OrcaTests(InputTests):

    def test_sp_SE(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'AM1',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP AM1
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'SMD',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
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
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'smd',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
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
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'I.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'SMD',
                'solvation_radii': 'SMD18',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        I 0.0 0.0 0.0
        *
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
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'I.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'smd',
                'solvation_radii': 'smd18',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        I 0.0 0.0 0.0
        *
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
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'octanol',
                'solvation_model': 'SMD',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
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
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'CPCM',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(chloroform)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvent_synonym(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'CHCL3',
                'solvation_model': 'SMD',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
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
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'octanol',
                'solvation_model': 'CPCM',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(octanol)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvation_octanol_cpcm2(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': '1-octanol',
                'solvation_model': 'CPCM',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G CPCM(octanol)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_solvation(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'ABC',
                }

        with self.assertRaises(Exception):
            inp = self.generate_calculation(**params)

    def test_sp_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_specifications(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'TightSCF',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP tightscf
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_multiple_specifications(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'TightSCF GRID6',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP tightscf grid6
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_duplicate_specifications(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'tightscf TightSCF GRID6',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP M062X Def2-SVP tightscf grid6
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_MP2(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'RI-MP2',
                'basis_set': 'cc-pVTZ',
                'charge': '-1',
                'specifications': 'cc-pVTZ/C',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP RI-MP2 cc-pVTZ cc-pvtz/c
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_SE(self):
        params = {
                'nproc': 8,
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'AM1',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT AM1
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_HF(self):
        params = {
                'nproc': 8,
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'charge': '-1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !OPT B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_SE(self):
        params = {
                'nproc': 8,
                'type': 'Frequency Calculation',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'AM1',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !FREQ AM1
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_HF(self):
        params = {
                'nproc': 8,
                'type': 'Frequency Calculation',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !FREQ HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Frequency Calculation',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'charge': '-1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !FREQ B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    #opt mod SE and HF

    def test_scan_bond_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_9_1.4_10/1_2;',
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
        B 0 1 = 9, 1.4, 10
        end
        end
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_opt_mod(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': '',
                }

        with self.assertRaises(Exception):
            inp = self.generate_calculation(**params)

    def test_no_method(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': '',
                'basis_set': '6-31+G(d,p)',
                'constraints': '',
                }

        with self.assertRaises(Exception):
            inp = self.generate_calculation(**params)

    def test_scan_angle_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_9_90_10/2_1_3;',
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
        A 1 0 2 = 9, 90, 10
        end
        end
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_dihedral_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_9_0_10/4_1_5_8;',
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
        D 3 0 4 7 = 9, 0, 10
        end
        end
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_bond_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Freeze/1_2;',
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
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_angle_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Freeze/2_1_3;',
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
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))


    def test_freeze_dihedral_DFT(self):
        params = {
                'nproc': 8,
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Freeze/4_1_5_8;',
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
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nmr_DFT(self):
        params = {
                'nproc': 8,
                'type': 'NMR Prediction',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'charge': '-1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !NMR B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_irrelevant_gen_bs(self):
        params = {
                'nproc': 8,
                'type': 'NMR Prediction',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'charge': '-1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'N=Def2-SVP;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !NMR B3LYP 6-31+G(d,p)
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts_DFT(self):
        params = {
                'nproc': 8,
                'type': 'TS Optimisation',
                'in_file': 'mini_ts.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
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
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

        #combination tests

    def test_ts_DFT_custom_bs(self):
        params = {
                'nproc': 8,
                'type': 'TS Optimisation',
                'in_file': 'mini_ts.xyz',
                'software': 'ORCA',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'N=Def2-SVP;',
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
        %basis
        newgto N
        S   5
        1      1712.8415853             -0.53934125305E-02
        2       257.64812677            -0.40221581118E-01
        3        58.458245853           -0.17931144990
        4        16.198367905           -0.46376317823
        5         5.0052600809          -0.44171422662
        S   1
        1         0.58731856571          1.0000000
        S   1
        1         0.18764592253          1.0000000
        P   3
        1        13.571470233           -0.40072398852E-01
        2         2.9257372874          -0.21807045028
        3         0.79927750754         -0.51294466049
        P   1
        1         0.21954348034          1.0000000
        D   1
        1         1.0000000              1.0000000
        end
        end
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT_custom_bs_ecp(self):
        params = {
                'nproc': 8,
                'type': 'Geometrical Optimisation',
                'in_file': 'Ph2I_cation.xyz',
                'software': 'ORCA',
                'charge': '+1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'I=Def2-TZVPD;',
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
        S   5
        1      5899.5791533              0.24188269271E-03
        2       898.54238765             0.15474041742E-02
        3       200.37237912             0.42836684457E-02
        4        31.418053840           -0.39417936275E-01
        5        15.645987838            0.96086691992
        S   2
        1        11.815741857            0.75961524091
        2         6.4614458287           0.42495501835
        S   1
        1         2.3838067579           1.0000000
        S   1
        1         1.1712089662           1.0000000
        S   1
        1         0.32115875757          1.0000000
        S   1
        1         0.12387919364          1.0000000
        S   1
        1         0.43491550641E-01      1.0000000
        P   3
        1       197.30030547             0.73951226905E-03
        2        20.061411349            0.66168450008E-01
        3         9.7631460485          -0.28554662348
        P   4
        1        12.984316904           -0.49096186164E-01
        2         3.6199503008           0.38914432482
        3         2.0232273090           0.65610817262
        4         1.0367490559           0.31803551647
        P   1
        1         0.45937816000          1.0000000
        P   1
        1         0.19116532928          1.0000000
        P   1
        1         0.74878813023E-01      1.0000000
        P   1
        1         0.21653491846E-01      1.0000000
        D   6
        1       119.12671745             0.82596039573E-03
        2        33.404240134            0.68377675770E-02
        3        17.805918203           -0.10308158997E-01
        4         4.8990510353           0.22670457658
        5         2.4516753106           0.44180113937
        6         1.1820693432           0.36775472225
        D   1
        1         0.52923110068          1.0000000
        D   1
        1         0.17000000000          1.0000000
        D   1
        1         0.61341708807E-01      1.0000000
        F   1
        1         2.1800000              1.0000000
        F   1
        1         0.44141808             1.0000000
        end

        NewECP I
          N_core 28
          lmax f
          s 7
           1     40.01583500    49.99429300 2
           2     17.42974700   281.02531700 2
           3      9.00548400    61.57332600 2
           4     19.45860900    21.84204000 2
           5     19.34926000    28.46819100 2
           6      4.82376700     0.24371300 2
           7      4.88431500     0.32080400 2
          p 8
           1     15.35546600    67.44284100 2
           2     14.97183300   134.88113700 2
           3      8.96016400    14.67505100 2
           4      8.25909600    29.37566600 2
           5     19.45860900    21.84204000 2
           6     19.34926000    28.46819100 2
           7      4.82376700     0.24371300 2
           8      4.88431500     0.32080400 2
          d 10
           1     15.06890800    35.43952900 2
           2     14.55532200    53.17605700 2
           3      6.71864700     9.06719500 2
           4      6.45639300    13.20693700 2
           5      1.19177900     0.08933500 2
           6      1.29115700     0.05238000 2
           7     19.45860900    21.84204000 2
           8     19.34926000    28.46819100 2
           9      4.82376700     0.24371300 2
           10     4.88431500     0.32080400 2
          f 4
           1     19.45860900   -21.84204000 2
           2     19.34926000   -28.46819100 2
           3      4.82376700    -0.24371300 2
           4      4.88431500    -0.32080400 2
        end

        end
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_NEB(self):
        params = {
                'nproc': 8,
                'type': 'Minimum Energy Path',
                'in_file': 'elimination_substrate.xyz',
                'auxiliary_file': 'elimination_product.xyz',
                'software': 'orca',
                'charge': -1,
                'method': 'gfn2-xtb',
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
        product "struct2.xyz"
        nimages 8
        end
        %pal
        nprocs 8
        end

        """
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_NEB2(self):
        params = {
                'nproc': 8,
                'type': 'Minimum Energy Path',
                'in_file': 'elimination_substrate.xyz',
                'auxiliary_file': 'elimination_product.xyz',
                'software': 'ORCA',
                'specifications': '--nimages 12',
                'charge': -1,
                'method': 'gfn2-xtb',
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
        product "struct2.xyz"
        nimages 12
        end
        %pal
        nprocs 8
        end

        """
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_hirshfeld_pop(self):
        params = {
                'nproc': 8,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': '--phirshfeld',
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
        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_mo(self):
        params = {
                'nproc': 8,
                'type': 'MO Calculation',
                'in_file': 'Ph2I_cation.xyz',
                'software': 'ORCA',
                'charge': '+1',
                'method': 'B3LYP',
                'basis_set': 'def2tzvp',
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

        %pal
        nprocs 8
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nproc(self):
        params = {
                'nproc': 1,
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'ORCA',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        !SP HF 3-21G
        *xyz -1 1
        Cl 0.0 0.0 0.0
        *
        %pal
        nprocs 1
        end
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))


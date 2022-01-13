from ccinput.packages.gaussian import GaussianCalculation
from ccinput.tests.testing_utilities import InputTests

class GaussianTests(InputTests):

    def test_sp_SE(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'AM1',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp AM1

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'SMD',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(SMD, Solvent=chloroform)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_SMD18(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'I.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'SMD',
                'solvation_radii': 'SMD18',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(SMD, Solvent=chloroform, Read)

        File created by ccinput

        -1 1
        I 0.0 0.0 0.0

        modifysph

        Br 2.60
        I 2.74

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_PCM_Bondi(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'chloroform',
                'solvation_model': 'PCM',
                'solvation_radii': 'Bondi',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(PCM, Solvent=chloroform, Read)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        Radii=bondi

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_PCM_UFF(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'PCM',
                'solvation_radii': 'UFF',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(PCM, Solvent=chloroform)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_CPCM_Bondi(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'CPCM',
                'solvation_radii': 'Bondi',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(CPCM, Solvent=chloroform, Read)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        Radii=bondi

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF_CPCM_UFF(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'Chloroform',
                'solvation_model': 'CPCM',
                'solvation_radii': 'UFF',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(CPCM, Solvent=chloroform)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvent_synonyms1(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'chcl3',
                'solvation_model': 'CPCM',
                'solvation_radii': 'UFF',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(CPCM, Solvent=chloroform)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_solvent_synonyms2(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                'solvent': 'meoh',
                'solvation_model': 'CPCM',
                'solvation_radii': 'UFF',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G SCRF(CPCM, Solvent=methanol)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_solvation(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
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
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT_specifications(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'nosymm 5D',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp M062X/Def2SVP nosymm 5d

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_superfluous_specifications(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'opt(loose)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.confirmed_specifications.strip(), '')


    def test_superfluous_specifications2(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'freq(noraman)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.confirmed_specifications.strip(), '')

    def test_duplicate_specifications(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'NoSymm NoSymm',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt M062X/Def2SVP nosymm

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_duplicate_specifications2(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'NoSymm NOSYMM',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt M062X/Def2SVP nosymm

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_SE(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'AM1',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt AM1

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_HF(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt HF/3-21G

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'AM1',
                'charge': '-1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt B3LYP/6-31+G(d,p)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.confirmed_specifications.strip(), '')

    def test_freq_SE(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Frequency Calculation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'AM1',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p freq AM1

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_HF(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Frequency Calculation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p freq HF/3-21G

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Frequency Calculation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'charge': '-1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p freq B3LYP/6-31+G(d,p)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    #opt mod SE and HF

    def test_scan_bond_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_9_1.4_10/1_2;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        B 1 2 S 10 0.03

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_angle_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_9_90_10/2_1_3;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        A 2 1 3 S 10 -1.95

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_dihedral_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_9_0_10/4_1_5_8;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        D 4 1 5 8 S 10 -17.99

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_dihedral_DFT2(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_2_0_10/2_1_5_8;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        D 2 1 5 8 S 10 -5.99

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_scan_dihedral_DFT3(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Scan_2_0_10/3_1_5_8;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        D 3 1 5 8 S 10 6.01

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_bond_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Freeze/1_2;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        B 1 2 F

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_angle_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Freeze/2_1_3;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        A 2 1 3 F

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_opt_mod(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
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
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': '',
                'basis_set': '6-31+G(d,p)',
                'constraints': '',
                }

        with self.assertRaises(Exception):
            inp = self.generate_calculation(**params)

    def test_freeze_dihedral_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Freeze/4_1_5_8;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        D 4 1 5 8 F

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freeze_dihedral_DFT2(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'constraints': 'Freeze/4_1_5_8;Freeze/1_2_3_4;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(modredundant) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        D 4 1 5 8 F
        D 1 2 3 4 F

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nmr_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'NMR Prediction',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'charge': '-1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p nmr B3LYP/6-31+G(d,p)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'TS Optimisation',
                'in_file': 'mini_ts.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(ts, NoEigenTest, CalcFC) B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        N   1.08764072053386     -0.33994563112543     -0.00972525479568
        H   1.99826836912112      0.05502842705407      0.00651240826058
        H   0.59453997172323     -0.48560162159600      0.83949232123172
        H   0.66998093862168     -0.58930117433261     -0.87511947469677

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts_DFT_df(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'TS Optimisation',
                'in_file': 'mini_ts.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': 'Def2SVP',
                'density_fitting': 'W06',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(ts, NoEigenTest, CalcFC) B3LYP/Def2SVP/W06

        File created by ccinput

        0 1
        N   1.08764072053386     -0.33994563112543     -0.00972525479568
        H   1.99826836912112      0.05502842705407      0.00651240826058
        H   0.59453997172323     -0.48560162159600      0.83949232123172
        H   0.66998093862168     -0.58930117433261     -0.87511947469677

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

        #combination tests

    def test_td_DFT(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'UV-Vis Calculation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p td M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_gen_bs(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'O=Def2-TZVPD;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt B3LYP/Gen

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        C H 0
        6-31+G(d,p)
        ****
        O     0
        S    6   1.00
          27032.3826310              0.21726302465D-03
           4052.3871392              0.16838662199D-02
            922.32722710             0.87395616265D-02
            261.24070989             0.35239968808D-01
             85.354641351            0.11153519115
             31.035035245            0.25588953961
        S    2   1.00
             12.260860728            0.39768730901
              4.9987076005           0.24627849430
        S    1   1.00
              1.1703108158           1.0000000
        S    1   1.00
              0.46474740994          1.0000000
        S    1   1.00
              0.18504536357          1.0000000
        S    1   1.00
              0.70288026270D-01      1.0000000
        P    4   1.00
             63.274954801            0.60685103418D-02
             14.627049379            0.41912575824D-01
              4.4501223456           0.16153841088
              1.5275799647           0.35706951311
        P    1   1.00
              0.52935117943           .44794207502
        P    1   1.00
              0.17478421270           .24446069663
        P    1   1.00
              0.51112745706D-01      1.0000000
        D    1   1.00
              2.31400000             1.0000000
        D    1   1.00
              0.64500000             1.0000000
        D    1   1.00
              0.14696477366          1.0000000
        F    1   1.00
              1.42800000             1.0000000
        ****

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_irrelevant_gen_bs(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'ethanol.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'Cl=Def2-TZVPD;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt B3LYP/6-31+G(d,p)

        File created by ccinput

        0 1
        C         -1.31970       -0.64380        0.00000
        H         -0.96310       -1.65260        0.00000
        H         -0.96310       -0.13940       -0.87370
        H         -2.38970       -0.64380        0.00000
        C         -0.80640        0.08220        1.25740
        H         -1.16150        1.09160        1.25640
        H         -1.16470       -0.42110        2.13110
        O          0.62360        0.07990        1.25870
        H          0.94410        0.53240        2.04240

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_genecp_bs(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Ph2I_cation.xyz',
                'software': 'Gaussian',
                'charge': '+1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'I=Def2-TZVPD;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt B3LYP/GenECP

        File created by ccinput

        1 1
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

        C H 0
        6-31+G(d,p)
        ****
        I     0
        S    5   1.00
           5899.5791533              0.24188269271D-03
            898.54238765             0.15474041742D-02
            200.37237912             0.42836684457D-02
             31.418053840           -0.39417936275D-01
             15.645987838            0.96086691992
        S    2   1.00
             11.815741857            0.75961524091
              6.4614458287           0.42495501835
        S    1   1.00
              2.3838067579           1.0000000
        S    1   1.00
              1.1712089662           1.0000000
        S    1   1.00
              0.32115875757          1.0000000
        S    1   1.00
              0.12387919364          1.0000000
        S    1   1.00
              0.43491550641D-01      1.0000000
        P    3   1.00
            197.30030547             0.73951226905D-03
             20.061411349            0.66168450008D-01
              9.7631460485          -0.28554662348
        P    4   1.00
             12.984316904           -0.49096186164D-01
              3.6199503008           0.38914432482
              2.0232273090           0.65610817262
              1.0367490559           0.31803551647
        P    1   1.00
              0.45937816000          1.0000000
        P    1   1.00
              0.19116532928          1.0000000
        P    1   1.00
              0.74878813023D-01      1.0000000
        P    1   1.00
              0.21653491846D-01      1.0000000
        D    6   1.00
            119.12671745             0.82596039573D-03
             33.404240134            0.68377675770D-02
             17.805918203           -0.10308158997D-01
              4.8990510353           0.22670457658
              2.4516753106           0.44180113937
              1.1820693432           0.36775472225
        D    1   1.00
              0.52923110068          1.0000000
        D    1   1.00
              0.17000000000          1.0000000
        D    1   1.00
              0.61341708807D-01      1.0000000
        F    1   1.00
              2.1800000              1.0000000
        F    1   1.00
              0.44141808             1.0000000
        ****

        I     0
        I-ECP     3     28
        f potential
          4
        2     19.45860900           -21.84204000
        2     19.34926000           -28.46819100
        2      4.82376700            -0.24371300
        2      4.88431500            -0.32080400
        s-f potential
          7
        2     40.01583500            49.99429300
        2     17.42974700           281.02531700
        2      9.00548400            61.57332600
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400
        p-f potential
          8
        2     15.35546600            67.44284100
        2     14.97183300           134.88113700
        2      8.96016400            14.67505100
        2      8.25909600            29.37566600
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400
        d-f potential
          10
        2     15.06890800            35.43952900
        2     14.55532200            53.17605700
        2      6.71864700             9.06719500
        2      6.45639300            13.20693700
        2      1.19177900             0.08933500
        2      1.29115700             0.05238000
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_genecp_bs_multiple_atoms(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Ph2I_cation.xyz',
                'software': 'Gaussian',
                'charge': '+1',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'I=Def2-TZVPD;H=Def2-TZVPD;C=Def2-TZVPD;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt B3LYP/GenECP

        File created by ccinput

        1 1
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

        I     0
        S    5   1.00
           5899.5791533              0.24188269271D-03
            898.54238765             0.15474041742D-02
            200.37237912             0.42836684457D-02
             31.418053840           -0.39417936275D-01
             15.645987838            0.96086691992
        S    2   1.00
             11.815741857            0.75961524091
              6.4614458287           0.42495501835
        S    1   1.00
              2.3838067579           1.0000000
        S    1   1.00
              1.1712089662           1.0000000
        S    1   1.00
              0.32115875757          1.0000000
        S    1   1.00
              0.12387919364          1.0000000
        S    1   1.00
              0.43491550641D-01      1.0000000
        P    3   1.00
            197.30030547             0.73951226905D-03
             20.061411349            0.66168450008D-01
              9.7631460485          -0.28554662348
        P    4   1.00
             12.984316904           -0.49096186164D-01
              3.6199503008           0.38914432482
              2.0232273090           0.65610817262
              1.0367490559           0.31803551647
        P    1   1.00
              0.45937816000          1.0000000
        P    1   1.00
              0.19116532928          1.0000000
        P    1   1.00
              0.74878813023D-01      1.0000000
        P    1   1.00
              0.21653491846D-01      1.0000000
        D    6   1.00
            119.12671745             0.82596039573D-03
             33.404240134            0.68377675770D-02
             17.805918203           -0.10308158997D-01
              4.8990510353           0.22670457658
              2.4516753106           0.44180113937
              1.1820693432           0.36775472225
        D    1   1.00
              0.52923110068          1.0000000
        D    1   1.00
              0.17000000000          1.0000000
        D    1   1.00
              0.61341708807D-01      1.0000000
        F    1   1.00
              2.1800000              1.0000000
        F    1   1.00
              0.44141808             1.0000000
        ****
        H     0
        S    3   1.00
             34.0613410              0.60251978D-02
              5.1235746              0.45021094D-01
              1.1646626              0.20189726
        S    1   1.00
              0.32723041             1.0000000
        S    1   1.00
              0.10307241             1.0000000
        P    1   1.00
              0.8000000              1.0000000
        P    1   1.00
              0.95774129632D-01      1.0000000
        ****
        C     0
        S    6   1.00
          13575.3496820              0.22245814352D-03
           2035.2333680              0.17232738252D-02
            463.22562359             0.89255715314D-02
            131.20019598             0.35727984502D-01
             42.853015891            0.11076259931
             15.584185766            0.24295627626
        S    2   1.00
              6.2067138508           0.41440263448
              2.5764896527           0.23744968655
        S    1   1.00
              0.57696339419          1.0000000
        S    1   1.00
              0.22972831358          1.0000000
        S    1   1.00
              0.95164440028D-01      1.0000000
        S    1   1.00
              0.48475401370D-01      1.0000000
        P    4   1.00
             34.697232244            0.53333657805D-02
              7.9582622826           0.35864109092D-01
              2.3780826883           0.14215873329
              0.81433208183          0.34270471845
        P    1   1.00
              0.28887547253           .46445822433
        P    1   1.00
              0.10056823671           .24955789874
        D    1   1.00
              1.09700000             1.0000000
        D    1   1.00
              0.31800000             1.0000000
        D    1   1.00
              0.90985336424D-01      1.0000000
        F    1   1.00
              0.76100000             1.0000000
        ****

        I     0
        I-ECP     3     28
        f potential
          4
        2     19.45860900           -21.84204000
        2     19.34926000           -28.46819100
        2      4.82376700            -0.24371300
        2      4.88431500            -0.32080400
        s-f potential
          7
        2     40.01583500            49.99429300
        2     17.42974700           281.02531700
        2      9.00548400            61.57332600
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400
        p-f potential
          8
        2     15.35546600            67.44284100
        2     14.97183300           134.88113700
        2      8.96016400            14.67505100
        2      8.25909600            29.37566600
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400
        d-f potential
          10
        2     15.06890800            35.43952900
        2     14.55532200            53.17605700
        2      6.71864700             9.06719500
        2      6.45639300            13.20693700
        2      1.19177900             0.08933500
        2      1.29115700             0.05238000
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_ecp(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'AuI.xyz',
                'software': 'Gaussian',
                'charge': '0',
                'method': 'B3LYP',
                'basis_set': '6-31+G(d,p)',
                'custom_basis_sets': 'I=Def2-TZVPD;Au=Def2-TZVPD;',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt B3LYP/GenECP

        File created by ccinput

        0 1
        Au        -9.27600       -1.06330        0.00000
        I         -6.60600       -1.06330        0.00000

        I     0
        S    5   1.00
        5899.5791533              0.24188269271D-03
        898.54238765             0.15474041742D-02
        200.37237912             0.42836684457D-02
        31.418053840           -0.39417936275D-01
        15.645987838            0.96086691992
        S    2   1.00
        11.815741857            0.75961524091
        6.4614458287           0.42495501835
        S    1   1.00
        2.3838067579           1.0000000
        S    1   1.00
        1.1712089662           1.0000000
        S    1   1.00
        0.32115875757          1.0000000
        S    1   1.00
        0.12387919364          1.0000000
        S    1   1.00
        0.43491550641D-01      1.0000000
        P    3   1.00
        197.30030547             0.73951226905D-03
        20.061411349            0.66168450008D-01
        9.7631460485          -0.28554662348
        P    4   1.00
        12.984316904           -0.49096186164D-01
        3.6199503008           0.38914432482
        2.0232273090           0.65610817262
        1.0367490559           0.31803551647
        P    1   1.00
        0.45937816000          1.0000000
        P    1   1.00
        0.19116532928          1.0000000
        P    1   1.00
        0.74878813023D-01      1.0000000
        P    1   1.00
        0.21653491846D-01      1.0000000
        D    6   1.00
        119.12671745             0.82596039573D-03
        33.404240134            0.68377675770D-02
        17.805918203           -0.10308158997D-01
        4.8990510353           0.22670457658
        2.4516753106           0.44180113937
        1.1820693432           0.36775472225
        D    1   1.00
        0.52923110068          1.0000000
        D    1   1.00
        0.17000000000          1.0000000
        D    1   1.00
        0.61341708807D-01      1.0000000
        F    1   1.00
        2.1800000              1.0000000
        F    1   1.00
        0.44141808             1.0000000
        ****
        Au     0
        S    3   1.00
        30.000000000            0.20749231108
        27.000000000           -0.33267893394
        14.746824331            0.38302817958
        S    1   1.00
        5.6017248938           1.0000000
        S    1   1.00
        1.3874162443           1.0000000
        S    1   1.00
        0.62923031957          1.0000000
        S    1   1.00
        0.14027517613          1.0000000
        S    1   1.00
        0.49379413761D-01      1.0000000
        P    4   1.00
        15.500000000            0.15001711880
        14.000000000           -0.23609813183
        6.4227368205           0.31458896948
        1.6595601681          -0.57279670446
        P    1   1.00
        0.79402913993          1.0000000
        P    1   1.00
        0.35125155397          1.0000000
        P    1   1.00
        0.11801737494          1.0000000
        P    1   1.00
        0.45000000000D-01      1.0000000
        D    4   1.00
        9.5524098656           0.40145559502D-01
        7.2698886937          -0.93690906606D-01
        1.7746496789           0.31746282317
        0.79960541055          0.46795192483
        D    1   1.00
        0.33252279372          1.0000000
        D    1   1.00
        0.12445133105          1.0000000
        F    1   1.00
        0.7248200              1.0000000
        ****

        I     0
        I-ECP     3     28
        f potential
        4
        2     19.45860900           -21.84204000
        2     19.34926000           -28.46819100
        2      4.82376700            -0.24371300
        2      4.88431500            -0.32080400
        s-f potential
        7
        2     40.01583500            49.99429300
        2     17.42974700           281.02531700
        2      9.00548400            61.57332600
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400
        p-f potential
        8
        2     15.35546600            67.44284100
        2     14.97183300           134.88113700
        2      8.96016400            14.67505100
        2      8.25909600            29.37566600
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400
        d-f potential
        10
        2     15.06890800            35.43952900
        2     14.55532200            53.17605700
        2      6.71864700             9.06719500
        2      6.45639300            13.20693700
        2      1.19177900             0.08933500
        2      1.29115700             0.05238000
        2     19.45860900            21.84204000
        2     19.34926000            28.46819100
        2      4.82376700             0.24371300
        2      4.88431500             0.32080400
        AU     0
        AU-ECP     3     60
        f potential
        2
        2      4.78982000            30.49008890
        2      2.39491000             5.17107381
        s-f potential
        4
        2     13.20510000           426.84667920
        2      6.60255000            37.00708285
        2      4.78982000           -30.49008890
        2      2.39491000            -5.17107381
        p-f potential
        4
        2     10.45202000           261.19958038
        2      5.22601000            26.96249604
        2      4.78982000           -30.49008890
        2      2.39491000            -5.17107381
        d-f potential
        4
        2      7.85110000           124.79066561
        2      3.92555000            16.30072573
        2      4.78982000           -30.49008890
        2      2.39491000            -5.17107381

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_global_specification(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'SCF(Tight)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp M062X/Def2SVP scf(tight)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.confirmed_specifications.strip(), 'scf(tight)')

    def test_multiple_global_specification(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'SCF(Tight) SCF(XQC)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp M062X/Def2SVP scf(tight,xqc)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.confirmed_specifications.strip(), 'scf(tight,xqc)')

    def test_multiple_global_specification2(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'SCF(Tight,XQC)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp M062X/Def2SVP scf(tight,xqc)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.confirmed_specifications.strip(), 'scf(tight,xqc)')

    def test_multiple_global_specification3(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'SCF(Tight, XQC)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp M062X/Def2SVP scf(tight,xqc)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))
        self.assertEqual(inp.confirmed_specifications.strip(), 'scf(tight,xqc)')


    def test_cmd_specification(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'opt(MaxStep=5)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(maxstep=5) M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_multiple_cmd_specifications(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'opt(MaxStep=5) opt(Tight)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(maxstep=5, tight) M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_both_specifications(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'opt(MaxStep=5) SCF(Tight)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(maxstep=5) M062X/Def2SVP scf(tight)

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_mixed(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'opt(MaxStep=5) opt(Tight) nosymm 5D',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(maxstep=5, tight) M062X/Def2SVP nosymm 5d

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_specifications_mixed2(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'opt(MaxStep=5,Tight) nosymm 5D',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt(maxstep=5, tight) M062X/Def2SVP nosymm 5d

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_invalid_specification(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'opt(MaxStep=5,Tight nosymm 5D',
                }

        with self.assertRaises(Exception):
            inp = self.generate_calculation(**params)

    def test_special_char(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Geometrical Optimisation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': '!#',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p opt M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_confirmed_specification_not_step(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'Constrained Optimisation',
                'constraints': 'Freeze/1_2_3;',
                'in_file': 'CH4.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '0',
                'specifications': 'pop(nbo)',
                }

        inp = self.generate_calculation(**params)

        self.assertEqual(inp.confirmed_specifications, 'pop(nbo)')

    def test_cmd_specification_td(self):
        params = {
                'nproc': 8,
                'mem': '10000MB',
                'type': 'UV-Vis Calculation',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'M06-2X',
                'basis_set': 'Def2-SVP',
                'charge': '-1',
                'specifications': 'td(NStates=5)',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p td(nstates=5) M062X/Def2SVP

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nproc(self):
        params = {
                'nproc': 1,
                'mem': '10000MB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=1
        %mem=10000MB
        #p sp HF/3-21G

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_mem_m(self):
        params = {
                'nproc': 8,
                'mem': '10000 m',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_mem_GB(self):
        params = {
                'nproc': 8,
                'mem': '10GB',
                'type': 'Single-Point Energy',
                'in_file': 'Cl.xyz',
                'software': 'Gaussian',
                'method': 'HF',
                'basis_set': '3-21G',
                'charge': '-1',
                }

        inp = self.generate_calculation(**params)

        REF = """
        %chk=calc.chk
        %nproc=8
        %mem=10000MB
        #p sp HF/3-21G

        File created by ccinput

        -1 1
        Cl 0.0 0.0 0.0

        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))


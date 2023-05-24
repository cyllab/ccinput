from ccinput.tests.testing_utilities import (
    GaussianCalculationTests,
    OrcaCalculationTests,
    XtbCalculationTests,
)


class GaussianSanityTests(GaussianCalculationTests):
    def test_sp_SE(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "Semi-empirical",
            "method": "AM1",
            "charge": "0",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

        E = self.run_calc(**params)
        self.assertTrue(self.known_energy(E, params, fail=True))

    def test_sp_SE2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "Semi-empirical",
            "method": "PM6",
            "charge": "0",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_SE3(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "Semi-empirical",
            "method": "PM7",
            "charge": "0",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "0",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_SMD(self):
        params = {
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_SMD2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
        }

        E = self.run_calc(**params)
        self.assertTrue(self.known_energy(E, params, fail=True))

    def test_sp_HF_SMD18(self):
        params = {
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
            "solvation_radii": "SMD18",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_PCM(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "0",
            "solvent": "Chloroform",
            "solvation_model": "PCM",
            "solvation_radii": "Bondi",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_CPCM(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "0",
            "solvent": "Chloroform",
            "solvation_model": "CPCM",
            "solvation_radii": "Bondi",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_PCM_UFF(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "0",
            "solvent": "Chloroform",
            "solvation_model": "PCM",
            "solvation_radii": "UFF",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_CPCM_UFF(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "Gaussian",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "CPCM",
            "solvation_radii": "UFF",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_DFT_specifications(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "Gaussian",
            "method": "DFT",
            "method": "M06-2X",
            "basis_set": "Def2SVP",
            "charge": "-1",
            "specifications": "nosymm 5D",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_DFT_specifications2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "Gaussian",
            "method": "DFT",
            "method": "M06-2X",
            "basis_set": "Def2SVP",
            "charge": "-1",
            "specifications": "nosymm 6D",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_gen_bs1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "DFT",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_gen_bs2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "CH4.xyz",
            "software": "Gaussian",
            "method": "DFT",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "C=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_genecp_bs1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Me2I_cation.xyz",
            "software": "Gaussian",
            "method": "HF",
            "charge": "+1",
            "basis_set": "STO-3G",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_genecp_bs2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Me2I_cation.xyz",
            "software": "Gaussian",
            "method": "HF",
            "charge": "+1",
            "basis_set": "STO-3G",
            "custom_basis_sets": "I=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_genecp_bs3(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Me2I_cation.xyz",
            "software": "Gaussian",
            "method": "HF",
            "charge": "+1",
            "basis_set": "STO-3G",
            "custom_basis_sets": "I=Def2-TZVPD;H=Def2-TZVP;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "Gaussian",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "Gaussian",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
            "custom_basis_sets": "I=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs3(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "Gaussian",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
            "custom_basis_sets": "Te=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs4(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "Gaussian",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
            "custom_basis_sets": "Te=Def2-TZVPD;I=Def2-TZVPP;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))


class OrcaSanityTests(OrcaCalculationTests):
    def test_sp_SE1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "method": "Semi-empirical",
            "method": "AM1",
            "charge": "0",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

        self.assertTrue(self.known_energy(E, params, fail=True))

    def test_sp_SE2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "method": "Semi-empirical",
            "method": "PM3",
            "charge": "0",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_SMD(self):
        params = {
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_SMD2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "I.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "SMD",
        }

        E = self.run_calc(**params)
        self.assertTrue(self.known_energy(E, params, fail=True))

    def test_sp_HF_SMD18(self):
        params = {
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

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_HF_CPCM(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "ORCA",
            "method": "HF",
            "basis_set": "3-21G",
            "charge": "-1",
            "solvent": "Chloroform",
            "solvation_model": "CPCM",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_gen_bs1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "method": "DFT",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_gen_bs2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "ethanol.xyz",
            "software": "ORCA",
            "method": "DFT",
            "charge": "0",
            "method": "B3LYP",
            "basis_set": "6-31+G(d,p)",
            "custom_basis_sets": "O=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_genecp_bs1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Me2I_cation.xyz",
            "software": "ORCA",
            "method": "HF",
            "charge": "+1",
            "basis_set": "STO-3G",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_genecp_bs2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Me2I_cation.xyz",
            "software": "ORCA",
            "method": "HF",
            "charge": "+1",
            "basis_set": "STO-3G",
            "custom_basis_sets": "I=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_genecp_bs3(self):
        params = {
            "type": "Single-Point Energy",
            "file": "Me2I_cation.xyz",
            "software": "ORCA",
            "method": "HF",
            "charge": "+1",
            "basis_set": "STO-3G",
            "custom_basis_sets": "I=Def2-TZVPD;H=Def2-TZVP;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs1(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "ORCA",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs2(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "ORCA",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
            "custom_basis_sets": "I=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs3(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "ORCA",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
            "custom_basis_sets": "Te=Def2-TZVPD;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_multiple_genecp_bs4(self):
        params = {
            "type": "Single-Point Energy",
            "file": "TeI2.xyz",
            "software": "ORCA",
            "method": "HF",
            "charge": "0",
            "basis_set": "STO-3G",
            "custom_basis_sets": "Te=Def2-TZVPD;I=Def2-TZVPP;",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))


class XtbSanityTests(XtbCalculationTests):
    def test_sp_gfn2(self):
        params = {
            "calc_name": "test",
            "type": "Single-Point Energy",
            "project": "New Project",
            "software": "xtb",
            "new_project_name": "SeleniumProject",
            "file": "ethanol.xyz",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_gfn2_explicit(self):
        params = {
            "calc_name": "test",
            "type": "Single-Point Energy",
            "project": "New Project",
            "software": "xtb",
            "new_project_name": "SeleniumProject",
            "file": "ethanol.xyz",
            "specifications": "--gfn2",
        }

        E = self.run_calc(**params)
        self.assertTrue(self.known_energy(E, params, fail=True))

    def test_sp_gfn1(self):
        params = {
            "calc_name": "test",
            "type": "Single-Point Energy",
            "project": "New Project",
            "software": "xtb",
            "new_project_name": "SeleniumProject",
            "file": "ethanol.xyz",
            "specifications": "--gfn1",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_gfn0(self):
        params = {
            "calc_name": "test",
            "type": "Single-Point Energy",
            "project": "New Project",
            "software": "xtb",
            "new_project_name": "SeleniumProject",
            "file": "ethanol.xyz",
            "specifications": "--gfn0",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_gfnff(self):
        params = {
            "calc_name": "test",
            "type": "Single-Point Energy",
            "project": "New Project",
            "software": "xtb",
            "new_project_name": "SeleniumProject",
            "file": "ethanol.xyz",
            "specifications": "--gfnff",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_GBSA(self):
        params = {
            "calc_name": "test",
            "type": "Single-Point Energy",
            "project": "New Project",
            "software": "xtb",
            "solvent": "Benzene",
            "solvation_model": "GBSA",
            "new_project_name": "SeleniumProject",
            "file": "ethanol.xyz",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

    def test_sp_ALPB(self):
        params = {
            "calc_name": "test",
            "type": "Single-Point Energy",
            "project": "New Project",
            "software": "xtb",
            "solvent": "Benzene",
            "solvation_model": "ALPB",
            "new_project_name": "SeleniumProject",
            "file": "ethanol.xyz",
        }

        E = self.run_calc(**params)
        self.assertFalse(self.known_energy(E, params))

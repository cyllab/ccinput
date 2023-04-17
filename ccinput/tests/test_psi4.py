from ccinput.tests.testing_utilities import InputTests
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation


class Psi4Tests(InputTests):
    # def test_sp_SE(self):
    #     params = {

    #     }

    #     inp = self.generate_calculation(**params)

    #     REF = """

    #     """

    #     self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_HF(self):
        params = {
            "type": "Energy",
            "file": "Cl.xyz",
            "software": "psi4",
            "method": "hf",
            "basis_set": "sto3g",
            "charge": "-1",
            "multiplicity": "1",
        }
        inp = self.generate_calculation(**params)

        REF = """
        molecule {
        -1 1  
        Cl   0.00000000   0.00000000   0.00000000

        }

        energy('scf/sto-3g')
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_DFT(self):
        params = {
            "type": "Energy",
            "file": "Cl.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "def2-tzvp",
            "charge": "-1",
            "multiplicity": "1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        molecule {
        -1 1  
        Cl   0.00000000   0.00000000   0.00000000

        }

        energy('b3lyp/def2-tzvp')
        """
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_HF(self):
        params = {
            "type": "Geometrical Optimisation",
            "file": "Cl.xyz",
            "software": "psi4",
            "method": "HF",
            "basis_set": "ccpVDZ",
            "charge": "-1",
            "multiplicity": "1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        molecule {
        -1 1  
        Cl   0.00000000   0.00000000   0.00000000

        }

        optimize('scf/cc-pVDZ')
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT(self):
        params = {
            "type": "Opt",
            "file": "Cl.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "631g",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        molecule {
        -1 1  
        Cl   0.00000000   0.00000000   0.00000000

        }

        optimize('b3lyp/6-31G')
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_unavailable_calc_type(self):
        params = {
            "type": "freq",
            "file": "Cl.xyz",
            "software": "psi4",
            "method": "HF",
            "basis_set": "cc-pvdz",
            "charge": "-1",
        }

        self.assertRaises(
            ImpossibleCalculation,
            msg=f"Calculation Type of {params['type']} not implemented yet for psi4",
        )

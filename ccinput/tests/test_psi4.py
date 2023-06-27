from ccinput.tests.testing_utilities import InputTests
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation
from ccinput.constants import CalcType


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
        #File created by ccinput

        memory 1000 mb

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
        #File created by ccinput

        memory 1000 mb

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
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "HF",
            "basis_set": "ccpVDZ",
            "charge": "0",
            "multiplicity": "1",
        }

        inp = self.generate_calculation(**params)

        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        optimize('scf/cc-pVDZ')
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_opt_DFT(self):
        params = {
            "type": "Opt",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "631g",
            "charge": "0",
        }

        inp = self.generate_calculation(**params)

        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        optimize('b3lyp/6-31G')
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_HF(self):
        params = {
            "type": "freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "HF",
            "basis_set": "ccpvdz",
        }

        inp = self.generate_calculation(**params)

        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        frequencies('scf/cc-pVDZ')
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_DFT(self):
        params = {
            "type": "freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "631g",
        }
        inp = self.generate_calculation(**params)
        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        frequencies('b3lyp/6-31G')
        """
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_optfreq_HF(self):
        params = {
            "type": "opt+freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "HF",
            "basis_set": "ccpvdz",
        }

        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        optimize('scf/cc-pVDZ')
        frequencies('scf/cc-pVDZ')
        """
        inp = self.generate_calculation(**params)
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_optfreq_DFT(self):
        params = {
            "type": "opt+freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "6-31G",
        }

        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        optimize('b3lyp/6-31G')
        frequencies('b3lyp/6-31G')
        """
        inp = self.generate_calculation(**params)
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_optfreq_single_specification(self):
        params = {
            "type": "opt+freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "6-31G",
            "specifications": "opt(return_wfn=True)",
        }

        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        optimize('b3lyp/6-31G', return_wfn=True)
        frequencies('b3lyp/6-31G')
        """
        inp = self.generate_calculation(**params)
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_optfreq_multispecs(self):
        params = {
            "type": "opt+freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "6-31G",
            "specifications": "opt(return_wfn=True) freq(dertype=1 return_wfn=True)",
        }

        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        optimize('b3lyp/6-31G', return_wfn=True)
        frequencies('b3lyp/6-31G', dertype=1, return_wfn=True)
        """
        inp = self.generate_calculation(**params)
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq_with_specs(self):
        params = {
            "type": "freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "6-31G",
            "specifications": "dertype=1 return_wfn=True",
        }
        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        frequencies('b3lyp/6-31G', dertype=1, return_wfn=True)
        """
        inp = self.generate_calculation(**params)
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_sp_memory_change(self):
        params = {
            "type": "energy",
            "mem": "500",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "HF",
            "basis_set": "ccpvdz",
        }

        REF = """
        #File created by ccinput

        memory 500 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        energy('scf/cc-pVDZ')
        """
        inp = self.generate_calculation(**params)
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_gen_calc_with_specs(self):
        params = {
            "type": "energy",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "6-31G",
            "specifications": "return_wfn=True",
        }
        REF = """
        #File created by ccinput

        memory 1000 mb

        molecule {
        0 1  
        O   -0.03589000  -0.05389000   0.00000000
        H   -0.33841000   0.84272000   0.00000000
        H    0.90805000   0.01251000   0.00000000
        }

        energy('b3lyp/6-31G', return_wfn=True)
        """
        inp = self.generate_calculation(**params)
        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_optfreq_specs_invalid(self):
        params = {
            "type": "optfreq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "6-31G",
            "specifications": "dertype=1 return_wfn=True opt(return_wfn=True)",
        }
        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_optfreq_specs_invalid_wihtout_calctype(self):
        params = {
            "type": "opt+freq",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "b3lyp",
            "basis_set": "6-31G",
            "specifications": "opt(return_wfn=True) (return_wfn=True)",
        }
        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_unavailable_calc_type(self):
        params = {
            "type": "TS",
            "file": "h2o.xyz",
            "software": "psi4",
            "method": "HF",
            "basis_set": "cc-pvdz",
        }

        self.assertRaises(
            ImpossibleCalculation,
            msg=f"Calculation Type of {params['type']} not implemented yet for psi4",
        )

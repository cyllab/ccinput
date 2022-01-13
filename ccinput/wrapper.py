import os

from ccinput.packages.gaussian import GaussianCalculation
from ccinput.packages.orca import OrcaCalculation

from ccinput.calculation import Calculation, Parameters
from ccinput.utilities import get_abs_type, get_abs_software

SOFTWARE_CLASSES = {
        'gaussian': GaussianCalculation,
        'orca': OrcaCalculation,
        }

def process_calculation(calc):
    cls = SOFTWARE_CLASSES[calc.parameters.software](calc)
    return cls

def generate_input(software, type, theory_level="", method="", basis_set="", \
            solvent="", solvation_model="", solvation_radii="",  specifications="", \
            density_fitting="", custom_basis_sets="", xyz="", in_file="", \
            constraints="", nproc=0, mem=0, charge=0, multiplicity=1, **kwargs):

    if xyz != "":
        xyz_structure = xyz # Verify that it is valid...
    elif in_file != "":
        if not os.path.isfile(in_file):
            raise Exception("Test input file not found: {}".format(in_file))###

        with open(in_file) as f:
            lines = f.readlines()

        xyz_structure = ''.join(lines)
    else:
        raise Exception("No input")

    abs_software = get_abs_software(software)
    if abs_software == -1:
        print("sfjlkkj")###
        return

    calc_type = get_abs_type(type)### error catch

    #verification/warnings (e.g. using default solvation model...)

    params = Parameters(abs_software, solvent, solvation_model, solvation_radii, \
            basis_set, theory_level, method, specifications, density_fitting, \
            custom_basis_sets, **kwargs)

    calc = Calculation(xyz_structure, params, calc_type, constraints=constraints, \
            nproc=nproc, mem=0, charge=charge, multiplicity=multiplicity)

    return process_calculation(calc)

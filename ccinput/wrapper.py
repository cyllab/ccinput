import os

from ccinput.packages.gaussian import GaussianCalculation
from ccinput.packages.orca import OrcaCalculation

from ccinput.calculation import Calculation, Parameters
from ccinput.utilities import get_abs_type, get_abs_software, standardize_xyz, parse_xyz_from_file
from ccinput.exceptions import *

SOFTWARE_CLASSES = {
        'gaussian': GaussianCalculation,
        'orca': OrcaCalculation,
        }

def process_calculation(calc):
    cls = SOFTWARE_CLASSES[calc.parameters.software](calc)
    return cls

def generate_calculation(software=None, type=None, method="", basis_set="", \
            solvent="", solvation_model="", solvation_radii="",  specifications="", \
            density_fitting="", custom_basis_sets="", xyz="", in_file="", \
            constraints="", nproc=0, mem="", charge=0, multiplicity=1, name="calc", \
            header="File created by ccinput", **kwargs):

    if software is None:
        raise InvalidParameter("Specify a software package to use (software=...)")

    if type is None:
        raise InvalidParameter("Specify a calculation type (type='...')")

    if xyz != "":
        xyz_structure = standardize_xyz(xyz)
    elif in_file != "":
        xyz_structure = parse_xyz_from_file(in_file)
    else:
        raise Exception("No input")

    abs_software = get_abs_software(software)

    calc_type = get_abs_type(type)

    params = Parameters(abs_software, solvent, solvation_model, solvation_radii, \
            basis_set, method, specifications, density_fitting, \
            custom_basis_sets, **kwargs)

    calc = Calculation(xyz_structure, params, calc_type, constraints=constraints, \
            nproc=nproc, mem=mem, charge=charge, multiplicity=multiplicity, \
            name=name, header=header, software=abs_software)

    return process_calculation(calc)

def gen_input(**args):
    return generate_calculation(**args).input_file

def write_input(filename, **args):
    inp = generate_calculation(**args).input_file
    with open(filename, 'w') as out:
        out.write(inp)



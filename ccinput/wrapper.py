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
        raise InvalidParameter("No input")

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

def cmd():
    import argparse
    parser = argparse.ArgumentParser(description='Generation an input for a computational chemistry package')
    parser.add_argument('software', help='Desired software package (Gaussian or ORCA)')

    parser.add_argument('type', help='Calculation type (opt, freq, sp, ...)')

    parser.add_argument('method', help='Computational method (HF, AM1, B3LYP, ...)')

    parser.add_argument('--basis_set', '-bs',  default="", type=str, help='Basis set')

    parser.add_argument('--solvent', '-s',  default="", type=str,
            help='Solvent for implicit solvation')

    parser.add_argument('--solvation_model', '-sm',  default="",
            type=str, help='Solvation model for implicit solvation')

    parser.add_argument('--solvation_radii', '-sr',  default="",
            type=str, help='Solvation radii for implicit solvation')

    parser.add_argument('--specifications', '-spec', default="",
            type=str, help='Additional commands')

    parser.add_argument('--density_fitting', '-df', default="",
            type=str, help='Basis set for density fitting (Gaussian only)')

    parser.add_argument('--custom_basis_sets', '-cbs', default="",
            type=str, help='String specification of custom basis sets on a per-atom basis')

    parser.add_argument('--xyz', '-x', default="",
            type=str, help='XYZ structure as string')

    parser.add_argument('--file', '-f', default="",
            type=str, help='XYZ structure as file')

    parser.add_argument('--output', '-o', default="",
            type=str, help='Write the result to the specified file')

    parser.add_argument('--constraints', '-co', default="",
            type=str, help='String specification of constraints for certain calculations')

    parser.add_argument('--nproc', '-n', default=1,
            type=int, help='Number of CPU cores to use')

    parser.add_argument('--mem', default="1000MB",
            help='Amount of memory to use (if no unit is provided, assumed to be in MB)')

    parser.add_argument('--charge', '-c', default=0, type=int, help='Charge of the system')

    parser.add_argument('--mult', '-m', default=1, type=int, help='Multiplicity of the system')

    parser.add_argument('--name', default="calc", type=str,
            help='Name of the produced file (unused by some packages)')

    parser.add_argument('--header', default="File created by ccinput", type=str,
            help='Header in produced file (unused by some packages)')

    args = parser.parse_args()

    inp = gen_input(software=args.software, type=args.type, method=args.method, \
            basis_set=args.basis_set, solvent=args.solvent, solvation_model=args.solvation_model, \
            solvation_radii=args.solvation_radii,  specifications=args.specifications, density_fitting=args.density_fitting, \
            custom_basis_sets=args.custom_basis_sets, xyz=args.xyz, in_file=args.file, constraints=args.constraints, \
            nproc=args.nproc, mem=args.mem, charge=args.charge, multiplicity=args.mult, \
            name=args.name, header=args.header)

    if args.output != "":
        with open(args.output, 'w') as out:
            out.write(inp)
        print("Input file written to {}".format(args.output))
    else:
        print(inp)

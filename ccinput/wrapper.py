import os
import sys
import shlex

from ccinput.__init__ import __version__
from ccinput.packages.gaussian import GaussianCalculation
from ccinput.packages.orca import OrcaCalculation
from ccinput.packages.xtb import XtbCalculation

from ccinput.calculation import (
    Calculation,
    Parameters,
    Constraint,
    parse_str_constraints,
    parse_freeze_constraints,
    parse_scan_constraints,
)
from ccinput.utilities import (
    get_abs_type,
    get_abs_software,
    standardize_xyz,
    parse_xyz_from_file,
    warn,
)
from ccinput.exceptions import *
from ccinput.presets import (
    save_preset,
    load_preset,
    list_presets,
    print_preset,
    is_preset,
)

SOFTWARE_CLASSES = {
    "gaussian": GaussianCalculation,
    "orca": OrcaCalculation,
    "xtb": XtbCalculation,
}


def process_calculation(calc):
    cls = SOFTWARE_CLASSES[calc.parameters.software](calc)
    return cls


def generate_calculation(
    software=None,
    type=None,
    method="",
    basis_set="",
    solvent="",
    solvation_model="",
    solvation_radii="",
    custom_solvation_radii="",
    specifications="",
    freeze=[],
    scan=[],
    sfrom=[],
    sto=[],
    snsteps=[],
    sstep=[],
    density_fitting="",
    custom_basis_sets="",
    xyz="",
    constraints="",
    nproc=1,
    mem=1000,
    charge=0,
    multiplicity=1,
    d3=False,
    d3bj=False,
    aux_name="calc2",
    name="calc",
    header="File created by ccinput",
    file=None,
    **kwargs,
):

    if software is None:
        raise InvalidParameter("Specify a software package to use")

    if type is None:
        raise InvalidParameter("Specify a calculation type")

    if method is None:
        raise InvalidParameter("Specify a calculation method")

    if xyz == "":
        raise InvalidParameter("No input structure")

    xyz_structure = standardize_xyz(xyz)

    abs_software = get_abs_software(software)

    calc_type = get_abs_type(type)

    params = Parameters(
        abs_software,
        solvent,
        solvation_model,
        solvation_radii,
        custom_solvation_radii,
        basis_set,
        method,
        specifications,
        density_fitting,
        custom_basis_sets,
        d3,
        d3bj,
        **kwargs,
    )

    _constraints = parse_str_constraints(
        constraints, xyz_structure, software=abs_software
    )
    _constraints += parse_freeze_constraints(
        freeze, xyz_structure, software=abs_software
    )
    _constraints += parse_scan_constraints(
        scan, sfrom, sto, snsteps, sstep, xyz_structure, software=abs_software
    )

    calc = Calculation(
        xyz_structure,
        params,
        calc_type,
        constraints=_constraints,
        nproc=nproc,
        mem=mem,
        charge=charge,
        multiplicity=multiplicity,
        aux_name=aux_name,
        name=name,
        header=header,
        software=abs_software,
        file=file,
    )

    return process_calculation(calc)


def gen_obj(**args):
    if "file" in args:
        if isinstance(args["file"], list):
            if len(args["file"]) > 1:
                print("file", args["file"])
                raise UnimplementedError(
                    "No support for multiple input files at once except from the command line"
                )
            elif len(args["file"]) == 0:
                pass
            else:
                xyz = parse_xyz_from_file(args["file"][0])
                args["xyz"] = xyz
        else:
            xyz = parse_xyz_from_file(args["file"])
            args["xyz"] = xyz

    return generate_calculation(**args)


def gen_input(**args):
    return gen_obj(**args).input_file


def write_input(filename, **args):
    inp = gen_input(**args)
    with open(filename, "w") as out:
        out.write(inp)


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generates an input for a computational chemistry package"
    )
    parser.add_argument(
        "software", nargs="?", help="Desired software package (Gaussian or ORCA)"
    )

    parser.add_argument("type", nargs="?", help="Calculation type (opt, freq, sp, ...)")

    parser.add_argument(
        "method", nargs="?", help="Computational method (HF, AM1, B3LYP, ...)"
    )

    parser.add_argument("--basis_set", "-bs", default="", type=str, help="Basis set")

    parser.add_argument(
        "--solvent", "-s", default="", type=str, help="Solvent for implicit solvation"
    )

    parser.add_argument(
        "--solvation_model",
        "-sm",
        default="",
        type=str,
        help="Solvation model for implicit solvation",
    )

    parser.add_argument(
        "--solvation_radii",
        "-sr",
        default="",
        type=str,
        help="Set of solvation radii for implicit solvation",
    )

    parser.add_argument(
        "--custom_solvation_radii",
        "-csr",
        default="",
        type=str,
        help="Specific solvation radii to modify",
    )

    parser.add_argument(
        "--specifications", "-spec", default="", type=str, help="Additional commands"
    )

    parser.add_argument(
        "--density_fitting",
        "-df",
        default="",
        type=str,
        help="Basis set for density fitting (Gaussian only)",
    )

    parser.add_argument(
        "--custom_basis_sets",
        "-cbs",
        default="",
        type=str,
        help="String specification of custom basis sets on a per-atom basis",
    )

    parser.add_argument(
        "--xyz", "-x", default="", type=str, help="XYZ structure as string"
    )

    parser.add_argument(
        "--file",
        "-f",
        default=[],
        nargs="+",
        type=str,
        help="XYZ structure(s) as file(s)",
    )

    parser.add_argument(
        "--output",
        "-o",
        default="",
        type=str,
        help="Write the result to the specified file (single file) or using the specified pattern (multiple files)",
    )

    parser.add_argument(
        "--constraints",
        "-co",
        default="",
        type=str,
        help="String specification of constraints for certain calculations",
    )

    parser.add_argument(
        "--freeze",
        action="append",
        nargs="+",
        default=[],
        metavar="ATOM",
        help="Freeze specified distance/angle/dihedral between the specified atoms",
    )

    parser.add_argument(
        "--scan",
        action="append",
        nargs="+",
        default=[],
        metavar="ATOM",
        help="Scan specified distance/angle/dihedral between the specified atoms (additional parameters required)",
    )

    parser.add_argument(
        "--from",
        action="append",
        dest="sfrom",
        default=[],
        metavar="FROM",
        help="Initial distance (for --scan)",
    )

    parser.add_argument(
        "--to",
        action="append",
        dest="sto",
        default=[],
        metavar="TO",
        help="Final distance (for --scan)",
    )

    parser.add_argument(
        "--nsteps",
        action="append",
        dest="snsteps",
        default=[],
        metavar="NSTEPS",
        help="Number of steps (for --scan)",
    )

    parser.add_argument(
        "--step",
        action="append",
        dest="sstep",
        default=[],
        metavar="STEP",
        help="Step size in Å (for --scan)",
    )

    parser.add_argument(
        "--nproc", "-n", default=1, type=int, help="Number of CPU cores to use"
    )

    parser.add_argument(
        "--mem",
        default="1000MB",
        help="Amount of memory to use (if no unit is provided, assumed to be in MB)",
    )

    parser.add_argument(
        "--charge", "-c", default=0, type=int, help="Charge of the system"
    )

    parser.add_argument(
        "--mult", "-m", default=1, type=int, help="Multiplicity of the system"
    )

    dispersion_group = parser.add_mutually_exclusive_group()

    dispersion_group.add_argument(
        "--d3", action="store_true", help="Use the D3(0) dispersion " "correction"
    )

    dispersion_group.add_argument(
        "--d3bj",
        action="store_true",
        help="Use the D3 dispersion " "correction with Becke-Johnson damping",
    )

    parser.add_argument(
        "--name",
        default="calc",
        type=str,
        help="Name of the produced file (unused by some packages)",
    )

    parser.add_argument(
        "--aux_name",
        default="calc2",
        type=str,
        help="Name of the auxiliary file (some calculation types only)",
    )

    parser.add_argument(
        "--header",
        default="File created by ccinput",
        type=str,
        help="Header in produced file (unused by some packages)",
    )

    parser.add_argument(
        "--save",
        type=str,
        help="Save the current parameters as a preset of the given name",
    )

    parser.add_argument(
        "--preset",
        type=str,
        nargs="?",
        const="",
        help="Load parameters from the chosen preset",
    )

    parser.add_argument(
        "--version", "-v", action="version", version=f"%(prog)s {__version__}"
    )
    return parser


def cmd(cmd_line=None):
    parser = get_parser()

    if cmd_line:
        args = parser.parse_args(shlex.split(cmd_line))
    else:
        args = parser.parse_args()

    if args.save:
        if args.preset:
            warn("Cannot both save and load a preset")
            exit(0)
        preset_name = save_preset(args, parser.parse_args(["", "", ""]))
        print_preset(preset_name)
        return

    if (args.preset and not is_preset(args.preset)) or args.preset == "":
        list_presets()
        exit(0)

    if args.preset:
        calcs, outputs = get_input_from_args(
            args, default_params=vars(parser.parse_args([]))
        )
    else:
        calcs, outputs = get_input_from_args(args)

    if args.output != "":
        for calc, outp in zip(calcs, outputs):
            with open(outp, "w") as out:
                out.write(calc.input_file)
            print(f"Input file written to {outp}")
    else:
        if len(calcs) == 1:
            print(calcs[0].output)
        else:
            for calc in calcs:
                n = max(int((39 - len(calc.calc.name)) / 2), 4)
                header = "-" * n + f" {calc.calc.name} " + "-" * n
                print(header)
                print(calc.output)
                print("\n\n")


def get_input_from_args(args, default_params=None):
    xyzs = []
    names = []
    outputs = []
    files = []

    if args.file:
        xyzs = [parse_xyz_from_file(f) for f in args.file]
        files = args.file
        if len(args.file) > 1 or args.name == "calc":
            names = [os.path.basename(f).split(".")[0] for f in args.file]
        else:
            names = [args.name]

        if args.output != "":
            head, tail = os.path.split(args.output)
            prefix, ext = tail.split(".")
            if prefix != "":
                prefix += "_"
            if len(args.file) > 1:
                outputs = [
                    os.path.join(head, prefix + name + "." + ext) for name in names
                ]
            else:
                outputs = [args.output]
    else:
        xyzs = [args.xyz]
        names = [args.name]
        outputs = [args.output]
        files = [args.file]

    params = {
        "software": args.software,
        "type": args.type,
        "method": args.method,
        "basis_set": args.basis_set,
        "solvent": args.solvent,
        "solvation_model": args.solvation_model,
        "solvation_radii": args.solvation_radii,
        "custom_solvation_radii": args.custom_solvation_radii,
        "specifications": args.specifications,
        "freeze": args.freeze,
        "scan": args.scan,
        "sfrom": args.sfrom,
        "sto": args.sto,
        "snsteps": args.snsteps,
        "sstep": args.sstep,
        "density_fitting": args.density_fitting,
        "custom_basis_sets": args.custom_basis_sets,
        "constraints": args.constraints,
        "nproc": args.nproc,
        "mem": args.mem,
        "charge": args.charge,
        "multiplicity": args.mult,
        "d3": args.d3,
        "d3bj": args.d3bj,
        "aux_name": args.aux_name,
        "header": args.header,
    }

    if args.preset:
        try:
            preset_params = load_preset(args.preset)
        except InvalidParameter as e:
            warn(str(e))
            exit(0)
        del preset_params["version"]
        for k, v in preset_params.items():
            if k == "specifications":
                params[k] = " ".join([v, params[k]])
            elif params[k] == default_params[k]:
                params[k] = v

    calcs = []
    for name, xyz, file in zip(names, xyzs, files):
        try:
            calc = gen_obj(name=name, xyz=xyz, file=file, **params)
        except CCInputException as e:
            print(f"!!! {str(e)} !!!")
            exit(0)
        else:
            calcs.append(calc)
    return calcs, outputs

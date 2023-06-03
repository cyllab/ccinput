import basis_set_exchange
import numpy as np
import re

from ccinput.utilities import (
    get_method,
    get_solvent,
    get_basis_set,
    clean_xyz,
    get_distance,
    get_angle,
    get_dihedral,
    get_npxyz,
    warn,
)
from ccinput.constants import (
    CalcType,
    ATOMIC_NUMBER,
    LOWERCASE_ATOMIC_SYMBOLS,
    SOFTWARE_MULTIPLICITY,
    SYN_METHODS,
)
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation


class NWChemCalculation:
    TEMPLATE = """TITLE "{}"
    start {}
    memory total {}
    charge {}

    geometry units angstroms
    {}end

    basis
    {}
    end

    {}{}{}
    {}
    """
    # Header
    # Name
    # Amount of memory
    # Charge
    # Geometry
    # Basis set
    # Method Block
    # Calculation Block
    # Additional Block
    # Tasks

    KEYWORDS = {
        CalcType.OPT: ["optimize"],
        #    CalcType.CONSTR_OPT: ["opt"],
        CalcType.TS: ["saddle"],
        CalcType.FREQ: ["freq"],
        CalcType.NMR: ["property"],
        CalcType.SP: ["energy"],
        #    CalcType.UVVIS: ["td"],
        #    CalcType.UVVIS_TDA: ["tda"],
        CalcType.OPTFREQ: ["optimize", "freq"],
    }
    BLOCK_KEYWORDS = {
        "method": [
            "sym",
            "adapt",
            "tol2e",
            "vectors",
            "tresh",
            "maxiter",
            "profile",
            "diis",
            "direct",
            "semidirect",
            "nr",
            "level",
            "convergence",
            "grid",
            "tolerances",
            "decomp",
            "odft",
            "incore",
            "iterations",
            "max_ovl",
            "cgmin",
            "rodft",
            "disp",
            "xdm",
            "noio",
            "print",
            "noprint",
            "fon",
            "nodisk",
            "fukui",
            "noscf",
        ],
        "calculation": [
            "a",  # TO DO
        ],
    }
    BLOCK_NAMES = {
        CalcType.OPT: "driver",
        CalcType.NMR: "property",
        CalcType.FREQ: "freq",
        CalcType.OPTFREQ: "driver",
        CalcType.TS: "driver",
    }

    def __init__(self, calc):
        self.calc = calc
        if self.calc.parameters.theory_level == "hf":  # Name of the block for HF is scf
            self.calc.parameters.theory_level = "scf"
        self.calc.mem = f"{self.calc.mem} mb"
        self.has_scan = False
        self.appendix = []
        self.command_line = ""
        self.method_block = f"{self.calc.parameters.theory_level}"
        self.calculation_block = ""
        self.additional_block = ""
        self.commands = {}
        self.solvation_radii = {}
        self.xyz_structure = ""
        self.tasks = ""
        self.input_file = ""

        if self.calc.type not in self.KEYWORDS:
            raise ImpossibleCalculation(
                f"NWChem does not support calculations of type {self.calc.type}"
            )
        self.handle_tasks()
        self.handle_specifications()
        self.handle_xyz()
        self.handle_basis_sets()
        self.handle_solvation()
        self.close_blocks()
        self.create_input_file()

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. "
        )
        return "".join([c for c in s if c in WHITELIST])

    def separate_lines(self, text):
        lines = text.split(";")
        for line in lines:
            line = self.clean(line)
        return "\n".join(lines)

    def handle_tasks(self):
        for word in self.KEYWORDS[self.calc.type]:
            self.tasks += f"task {self.calc.parameters.theory_level} {word} \n"
        # handle levels of theory
        if (
            self.calc.parameters.method == "hf"
            or self.calc.parameters.method in SYN_METHODS["hf"]
        ):
            scf_block = f"""
            {SOFTWARE_MULTIPLICITY['nwchem'][self.calc.multiplicity]}
            """
            self.method_block += scf_block
        if self.calc.parameters.theory_level == "dft":
            dft_block = f"""
            xc {self.calc.parameters.method}
            mult {self.calc.multiplicity}
            """
            self.method_block += dft_block
            if self.calc.parameters.d3:
                self.method_block += "disp vdw 3 \n"
            elif self.calc.parameters.d3bj:
                self.method_block += "disp vdw 4 \n"

    def handle_basis_sets(self):
        basis_set = get_basis_set(self.calc.parameters.basis_set, "nwchem")
        if basis_set != "":
            self.basis_set = f"* library {basis_set}"

    def handle_specifications(self):
        if self.calc.parameters.specifications != "":
            s = self.separate_lines(self.calc.parameters.specifications)
            for spec in s.split("\n"):
                matched = re.search(r".*\((.*)\)", spec)
                if matched == None:
                    self.additional_block += f"{spec} \n"
                else:
                    command = matched.group(1)
                    block_name = spec[: matched.span(1)[0] - 1]
                    if block_name == "scf":
                        self.method_block += f"{command} \n"
                    elif block_name == "opt" or block_name == "ts":
                        if self.calculation_block == "":
                            self.calculation_block += f"\n driver \n"
                        self.calculation_block += f"{command} \n"
                    elif block_name == "nmr":
                        if self.calculation_block == "":
                            self.calculation_block += f" \n property \n"
                        self.calculation_block += f"{command} \n"
        if self.additional_block != "":
            self.additional_block = "\n" + self.additional_block

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def handle_solvation(self):
        return

    def close_blocks(self):
        if self.method_block != "":
            self.method_block += " end \n"
        if self.calculation_block != "":
            self.calculation_block += " end \n"

    def create_input_file(self):
        raw = self.TEMPLATE.format(
            self.calc.header,
            self.calc.name,
            self.calc.mem,
            self.calc.charge,
            self.xyz_structure,
            self.basis_set,
            self.method_block,
            self.calculation_block,
            self.additional_block,
            self.tasks,
        )
        self.input_file = "\n".join([i.strip() for i in raw.split("\n")])

    @property
    def output(self):
        return self.input_file

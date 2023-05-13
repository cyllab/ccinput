import basis_set_exchange
import numpy as np

from ccinput.utilities import (
    get_method,
    get_solvent,
    get_basis_set,
    clean_xyz,
    get_distance,
    get_angle,
    get_dihedral,
    get_npxyz,
)
from ccinput.constants import CalcType, ATOMIC_NUMBER, LOWERCASE_ATOMIC_SYMBOLS, SOFTWARE_MULTIPLICITY, SYN_METHODS
multiplicity_dict = SOFTWARE_MULTIPLICITY['nwchem']
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
    {}
    {}
    """
    # Header
    # Name
    # Amount of memory
    # Charge
    # Geometry
    # Basis set
    # Additional blocks
    # Tasks

    KEYWORDS = {
        CalcType.OPT: ["optimize"],
    #    CalcType.CONSTR_OPT: ["opt"],
        CalcType.TS: ["saddle"],
        CalcType.FREQ: ["freq"],
    #    CalcType.NMR: ["nmr"],
        CalcType.SP: "energy",
    #    CalcType.UVVIS: ["td"],
    #    CalcType.UVVIS_TDA: ["tda"],
        CalcType.OPTFREQ: ["opt", "freq"],
    }

    # Get a set of all unique calculation keywords
    KEYWORD_LIST = set([j for i in KEYWORDS.values() for j in i])

    CALC_TYPES = [
        CalcType.SP,
        # CalcType.OPT,
        # CalcType.CONSTR_OPT,
        # CalcType.FREQ,
        # CalcType.TS,
        # CalcType.OPTFREQ,
    ]

    def __init__(self, calc):
        self.calc = calc
        if self.calc.parameters.theory_level == 'hf' :
            self.calc.parameters.theory_level = 'scf'
        self.calc.mem = f"{self.calc.mem} mb"
        self.has_scan = False
        self.appendix = []
        self.command_line = ""
        self.Blocks=f"{self.calc.parameters.theory_level}"
        self.commands = {}
        self.solvation_radii = {}

        self.confirmed_specifications = ""
        self.xyz_structure = ""
        self.input_file = ""

        if self.calc.type not in self.KEYWORDS:
            raise ImpossibleCalculation(
                f"NWChem does not support calculations of type {self.calc.type}"
            )
        self.handle_command()
        self.handle_specifications()
        self.handle_xyz()
        self.handle_solvation()
        self.create_input_file()

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. "
        )
        return "".join([c for c in s if c in WHITELIST])

    def add_option(self, key, option):
        _option = option.strip()
        if key in self.commands:
            if _option not in self.commands[key]:
                self.commands[key].append(_option)
        else:
            self.commands[key] = [_option]

    def add_options(self, key, options):
        for option in options:
            self.add_option(key, option)

    def handle_specifications(self):
        s = ""
        specs = self.calc.parameters.specifications.split(';')
        for spec in specs :
            spec = self.clean(spec)
        s = ";".join(specs)
        if(s != '') :
            self.Blocks+='\n'
        self.Blocks+=f"""{s}
        end
        """
        return

    def handle_command(self):
        self.tasks = f'task {self.calc.parameters.theory_level} {self.KEYWORDS[self.calc.type]}'
        basis_set = get_basis_set(self.calc.parameters.basis_set, "nwchem")
        if(basis_set != ''):
            self.basis_set = f"* library {basis_set}"
        if(self.calc.parameters.method == 'hf' or self.calc.parameters.method in SYN_METHODS['hf']) :
            scf_block = f"""
            {multiplicity_dict[self.calc.multiplicity]}
            """
            self.Blocks += scf_block
        if(self.calc.parameters.theory_level == 'dft') :
            dft_block = f"""
            xc {self.calc.parameters.method}
            """
            self.Blocks += dft_block
            if(self.calc.parameters.d3 == True) :
                self.Blocks += "disp vdw 3"
            elif(self.calc.parameters.d3bj == True) :
                self.Blocks += "disp vdw 4"
        return

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def handle_solvation(self):
        return

    def create_input_file(self):
        raw = self.TEMPLATE.format(
            self.calc.header,
            self.calc.name,
            self.calc.mem,
            self.calc.charge,
            self.xyz_structure,
            self.basis_set,
            self.Blocks,
            self.tasks,
        )
        self.input_file = "\n".join([i.strip() for i in raw.split("\n")]).replace(
            "\n\n\n", "\n\n"
        )
    
    @property
    def output(self):
        return self.input_file

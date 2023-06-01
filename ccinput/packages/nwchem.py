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
from ccinput.constants import (
    CalcType,
    ATOMIC_NUMBER,
    LOWERCASE_ATOMIC_SYMBOLS,
    SOFTWARE_MULTIPLICITY,
    SYN_METHODS
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
    {}
    {}
    {}"""
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
        CalcType.SP: ["energy"],
    #    CalcType.UVVIS: ["td"],
    #    CalcType.UVVIS_TDA: ["tda"],
        CalcType.OPTFREQ: ["optimize", "freq"],
    }

    def __init__(self, calc):
        self.calc = calc
        if self.calc.parameters.theory_level == 'hf' : # Name of the block for HF is scf
            self.calc.parameters.theory_level = 'scf'
        self.calc.mem = f"{self.calc.mem} mb"
        self.has_scan = False
        self.appendix = []
        self.command_line = ""
        self.Blocks=f"{self.calc.parameters.theory_level}"
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
        self.handle_additional()
        self.create_input_file()

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. "
        )
        return "".join([c for c in s if c in WHITELIST])

    def separate_lines(self,text):
        lines = text.split(';')
        for line in lines :
            line = self.clean(line)
        return "\n".join(lines)
    
    def handle_tasks(self):
        for word in self.KEYWORDS[self.calc.type]:
            self.tasks += f'task {self.calc.parameters.theory_level} {word} \n'
        if(self.calc.parameters.method == 'hf' or self.calc.parameters.method in SYN_METHODS['hf']) :
            scf_block = f"""
            {SOFTWARE_MULTIPLICITY['nwchem'][self.calc.multiplicity]}
            """
            self.Blocks += scf_block
        if(self.calc.parameters.theory_level == 'dft') :
            dft_block = f"""
            xc {self.calc.parameters.method}
            mult {self.calc.multiplicity}
            """
            self.Blocks += dft_block
            if(self.calc.parameters.d3 == True) :
                self.Blocks += "disp vdw 3"
            elif(self.calc.parameters.d3bj == True) :
                self.Blocks += "disp vdw 4"
    def handle_basis_sets(self):
        basis_set = get_basis_set(self.calc.parameters.basis_set, "nwchem")
        if(basis_set != ''):
            self.basis_set = f"* library {basis_set}"
    def handle_specifications(self):
            s = self.separate_lines(self.calc.parameters.specifications)
            if(s != '') :
                s+= '\n'
                self.Blocks+='\n'
            self.Blocks+=f"{s}end"

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def handle_solvation(self):
        return
    
    def handle_additional(self):
        self.additional = self.separate_lines(self.calc.additional)

    def create_input_file(self):
        raw = self.TEMPLATE.format(
            self.calc.header,
            self.calc.name,
            self.calc.mem,
            self.calc.charge,
            self.xyz_structure,
            self.basis_set,
            self.Blocks,
            self.additional,
            self.tasks,
        )
        self.input_file = "\n".join([i.strip() for i in raw.split("\n")])
    
    @property
    def output(self):
        return self.input_file

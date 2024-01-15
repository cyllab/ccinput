__author__ = "Zarko Ivkovic, zivkoviv7@alumnes.ub.edu"
import basis_set_exchange
import re

from ccinput.utilities import (
    get_solvent,
    get_basis_set,
    clean_xyz,
    warn,
)
from ccinput.constants import (
    CalcType,
    ATOMIC_NUMBER,
    LOWERCASE_ATOMIC_SYMBOLS,
    SOFTWARE_MULTIPLICITY,
)
from ccinput.exceptions import (
    InvalidParameter,
    ImpossibleCalculation,
    UnimplementedError,
    MissingParameter,
)


class NWChemCalculation:
    TEMPLATE = """TITLE "{}"
    start {}
    memory total {}
    charge {}

    geometry units angstroms noautosym
    {}end

    {}

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
        CalcType.CONSTR_OPT: ["optimize"],
        CalcType.TS: ["saddle"],
        CalcType.FREQ: ["freq"],
        CalcType.NMR: ["property"],
        CalcType.SP: ["energy"],
        CalcType.OPTFREQ: ["optimize", "freq"],
        CalcType.MEP: ["neb ignore"],
    }

    def __init__(self, calc):
        self.calc = calc
        self.calc.mem = f"{int(self.calc.mem/self.calc.nproc)} mb"
        self.has_scan = False
        self.appendix = []
        self.command_line = ""
        self.method_block = ""
        self.calculation_block = ""
        self.additional_block = ""
        self.commands = {}
        self.solvation_radii = {}
        self.xyz_structure = ""
        self.tasks = ""
        self.radii_parameters = ""
        self.input_file = ""
        # Some specific syntax processing related only to nwchem
        if (
            self.calc.parameters.theory_level == "hf"
            or self.calc.parameters.method == "uhf"
            or self.calc.parameters.method == "rhf"
        ):  # Name of the block for HF is scf
            self.calc.parameters.theory_level = "scf"
        elif self.calc.parameters.theory_level == "cc":
            # Name of any coupled cluster block is ccsd
            self.calc.parameters.theory_level = "ccsd"
        self.method_block = f"{self.calc.parameters.theory_level}"

        if self.calc.type not in self.KEYWORDS:
            raise ImpossibleCalculation(
                f"NWChem does not support calculations of type {self.calc.type}"
            )
        self.handle_tasks()
        self.handle_solvation()
        self.handle_specifications()
        self.handle_xyz()
        self.handle_basis_sets()
        self.close_blocks()
        self.create_input_file()

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,._ "
        )
        return "".join([c for c in s if c in WHITELIST])

    def separate_lines(self, text):
        text = text.replace(")", ");")
        lines = text.split(";")
        clean = []
        for line in lines:
            if line.strip() != "":
                clean.append(self.clean(line.lower().strip()))
        return "\n".join(clean)

    def handle_tasks(self):
        for word in self.KEYWORDS[self.calc.type]:
            # different cc methods are distinguished in tasks section
            if self.calc.parameters.theory_level in ["ccsd", "mp2"]:
                keyword = self.calc.parameters.method
            else:
                keyword = self.calc.parameters.theory_level
            self.tasks += f"task {keyword} {word} \n"
        # handle levels of theory
        if self.calc.parameters.theory_level == "scf":
            scf_block = "\n"
            if self.calc.parameters.method != "hf":
                scf_block += f"{self.calc.parameters.method} \n"
            scf_block += f"{SOFTWARE_MULTIPLICITY['nwchem'][self.calc.multiplicity]} \n"
            self.method_block += scf_block
        elif self.calc.parameters.theory_level == "dft":
            dft_block = f"""
            xc {self.calc.parameters.method}
            mult {self.calc.multiplicity}
            """
            self.method_block += dft_block
            # dispersion
            if self.calc.parameters.d3:
                self.method_block += "disp vdw 3 \n"
            elif self.calc.parameters.d3bj:
                self.method_block += "disp vdw 4 \n"
        elif self.calc.parameters.theory_level == "mcscf":
            self.method_block += f" \n multiplicity {self.calc.multiplicity} \n"
        elif self.calc.parameters.theory_level in ["mp2", "ccsd"]:
            self.method_block += "\n"
        if self.calc.type == CalcType.NMR:
            self.calculation_block += f" \n property \n shielding \n"
        if (
            self.calc.parameters.method == "rimp2"
            and self.calc.parameters.density_fitting == ""
        ):
            raise InvalidParameter("RI-MP2 calculation requires auxilary basis set")

    def handle_basis_sets(self):
        basis_set = get_basis_set(self.calc.parameters.basis_set, "nwchem")
        if basis_set != "" and len(self.calc.parameters.custom_basis_sets) == 0:
            self.basis_set = f"basis \n * library {basis_set} \n end"
        elif self.calc.parameters.custom_basis_sets != "":
            self.parse_custom_basis_set(basis_set)
        else:
            raise MissingParameter(
                "You must specify a basis set for a nwchem calculation"
            )
        # Density fitting
        if self.calc.parameters.density_fitting != "":
            density_fitting = get_basis_set(
                self.calc.parameters.density_fitting, "nwchem"
            )
            self.basis_set += (
                f'\n \n basis "cd basis" \n * library {density_fitting} \n end'
            )

    def parse_custom_basis_set(self, base_bs):
        custom_basis_sets = self.calc.parameters.custom_basis_sets
        to_append_bs = []
        to_append_ecp = []
        not_recoginzed_bs = {}
        unique_atoms = []
        normal_atoms = []
        for line in self.calc.xyz.split("\n"):
            if line.strip() == "":
                continue
            a, *_ = line.split()
            if a not in unique_atoms:
                unique_atoms.append(a)
                if a not in normal_atoms and a not in custom_basis_sets:
                    normal_atoms.append(a)

        custom_atoms = []
        for el, bs_keyword in custom_basis_sets.items():
            if el not in unique_atoms:
                continue

            custom_atoms.append(el)
            try:
                el_num = ATOMIC_NUMBER[el]
            except KeyError:
                raise InvalidParameter(
                    f"Invalid atom in custom basis set string: '{el}'"
                )

            try:
                bs = basis_set_exchange.get_basis(
                    bs_keyword,
                    fmt="nwchem",
                    elements=[el_num],
                    header=False,
                    optimize_general=True,
                    uncontract_general=True,
                )

            except:
                # Some basis sets are built-in, but use different names as the BSE does (e.g., SDD)
                # In this case, just feed the user keyword in and hope it works.
                # The basis set string has been recognized by ccinput, so it should exist in the program.
                # ECP is NOT added in this case, and the user will recieve a warning
                warn(
                    f"Basis set {bs_keyword} couldn't be pulled from basis set exchange. The ECP will not be added to this basis set."
                )
                not_recoginzed_bs[el] = bs_keyword
            else:
                matched_ECP = re.search(r"ECP\n(.*?)END", bs, re.DOTALL)
                matched_bs = re.search(
                    r'BASIS "ao basis" SPHERICAL PRINT\n(.*?)END', bs, re.DOTALL
                )
                if matched_ECP != None:
                    to_append_ecp.append(matched_ECP.group(1))
                to_append_bs.append(matched_bs.group(1))
        if len(custom_atoms) > 0:
            custom_bs = "\n".join(to_append_bs)
            self.basis_set = f"basis spherical \n {custom_bs} \n "
            if len(normal_atoms) > 0:
                self.basis_set += (
                    f"* library {base_bs} except {' '.join(custom_atoms)} \n"
                )
        if len(not_recoginzed_bs) > 0:
            if self.basis_set == None:
                self.basis_set = "basis\n"
            for element in not_recoginzed_bs:
                self.basis_set += f"{element} library {not_recoginzed_bs[element]} \n"
        try:
            self.basis_set += "end"
        except:
            self.basis_set = f"basis \n * library {base_bs} \n end"
        if len(to_append_ecp) > 0:
            self.basis_set += "\n \n ecp \n"
            self.basis_set += "\n".join(to_append_ecp)
            self.basis_set += " end"

    def handle_specifications(self):
        if self.clean(self.calc.parameters.specifications).strip() != "":
            temp = "\n"  # Here we will store frequency related specifiations in case of FREQOPT calculations
            s = self.separate_lines(self.calc.parameters.specifications)
            # Could be more sophisticated to catch other incorrect specifications
            if s.count("(") != s.count(")"):
                raise InvalidParameter(
                    "Invalid specifications: parenthesis not matching"
                )
            for spec in s.split("\n"):
                # format of the specifications is BLOCK_NAME1(command1);BLOCK_NAME2(command2);...
                matched = re.search(r".*\((.*)\)", spec)
                if matched == None:
                    # To make a difference between neb(defualt mep method) and freezing string method
                    # User has to put some of the following keyword as specification, independant of what calculation was specified in input
                    if spec in [
                        "string",
                        "freezing string sethod",
                        "fsm",
                        "freezing string",
                    ]:
                        self.tasks = self.tasks.replace("neb", "string")
                    else:
                        self.additional_block += f"{spec} \n"
                else:
                    command = matched.group(1)
                    if command.find(",") != -1:
                        command = command.replace(",", "\n").strip()
                    block_name = spec[: matched.span(1)[0] - 1]
                    if block_name == "scf" or block_name == "dft" or block_name == "hf":
                        if (
                            command == "adft"
                            and self.calc.parameters.density_fitting == ""
                        ):
                            raise InvalidParameter(
                                "adft keyword requires auxilary basis set"
                            )
                        self.method_block += f"{command} \n"
                    elif (block_name == "opt" or block_name == "ts") and (
                        self.calc.type
                        in [
                            CalcType.CONSTR_OPT,
                            CalcType.OPT,
                            CalcType.TS,
                            CalcType.OPTFREQ,
                        ]
                    ):
                        if self.calculation_block == "":
                            self.calculation_block += f"\n driver \n"
                        self.calculation_block += f"{command} \n"
                    elif block_name == "nmr" and self.calc.type == CalcType.NMR:
                        self.calculation_block += f"{command} \n"
                    elif block_name == "freq" and self.calc.type == CalcType.FREQ:
                        if self.calculation_block == "":
                            self.calculation_block += f"\n freq \n"
                        self.calculation_block += f"{command} \n"
                    elif block_name == "freq" and self.calc.type == CalcType.OPTFREQ:
                        temp += f"{command} \n"
                    elif (
                        block_name in ["neb", "string", "fsm", "mep"]
                        and self.calc.type == CalcType.MEP
                    ):
                        if self.calculation_block == "":
                            self.calculation_block += f"\n neb \n"
                        self.calculation_block += f"{command} \n"
                    elif (
                        block_name == "sol"
                        or block_name == "cosmo"
                        or block_name == "smd"
                    ):
                        self.additional_block = self.additional_block.replace(
                            "cosmo \n", f"cosmo \n {command} \n"
                        )
                    elif (
                        block_name == "mp2"
                        and self.calc.parameters.theory_level == "mp2"
                    ):
                        self.method_block += f"{command} \n"
                    elif (
                        block_name == "cc"
                        and self.calc.parameters.theory_level == "ccsd"
                    ):
                        self.method_block += f"{command} \n"
                    elif (
                        block_name in ["mcscf", "casscf"]
                        and self.calc.type == CalcType.SP
                    ):
                        self.method_block += f"{command} \n"
            if temp != "\n":
                self.additional_block += f"\n freq {temp} end \n"
        if self.tasks.find("string") != -1:
            self.calculation_block = self.calculation_block.replace("neb", "string")
        # Check if there are necessary specifications for mcscf calculation:
        if (
            self.method_block.find("active") == -1
            or self.method_block.find("actelec") == -1
        ) and self.calc.parameters.theory_level == "mcscf":
            raise MissingParameter(
                "You must specify the number of active orbitals and active electrons in CASSCF calculation"
            )
        # Check if there is density fitting for DFT but adft keyword is not specified
        if (
            self.calc.parameters.theory_level == "dft"
            and self.method_block.find("adft") == -1
            and self.calc.parameters.density_fitting != ""
        ):
            self.method_block += f"adft \n"
        # Handle contraints
        if self.calc.type == CalcType.CONSTR_OPT:
            if len(self.calc.constraints) == 0:
                raise InvalidParameter("No constraint in constrained optimisation mode")
            constraints = ""
            for constraint in self.calc.constraints:
                constraints += constraint.to_nwchem()
            if constraints != "":
                self.additional_block += (
                    f"\n geometry adjust \n zcoord \n {constraints} end \n end \n"
                )

        # Handle scans (TO_DO)
        if self.calc.type == CalcType.CONSTR_OPT:
            for constraint in self.calc.constraints:
                if constraint.scan:
                    self.additional_block += constraint.to_nwchem()

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def handle_solvation(self):
        """Default radii used in nwchem are complex combination of different sources.
        More info can be found: https://nwchemgit.github.io/Solvation-Models
        """
        if self.calc.parameters.solvent.lower() not in ["", "vacuum"]:
            solvation_model = self.calc.parameters.solvation_model.lower()
            solvent_keyword = get_solvent(
                self.calc.parameters.solvent,
                self.calc.parameters.software,
                solvation_model=self.calc.parameters.solvation_model,
            )
            self.additional_block += (
                f"\n cosmo \n minbem 3 \n ificos 1 \n solvent {solvent_keyword} \n"
            )
            """ Best grids recommended by Marenich, A. V.; Cramer, C. J.; Truhlar, D. G.
            The Journal of Physical Chemistry B 2009, 113 (18), 6378-6396. https://doi.org/10.1021/jp810292n."""
            if solvation_model in ["cosmo"]:
                pass
            elif solvation_model in ["smd", "smd18"]:
                self.additional_block += "do_cosmo_smd \n"
                if self.calc.parameters.theory_level != "dft":
                    raise UnimplementedError(
                        "SMD model is only available with DFT in nwchem."
                    )
            else:
                raise UnimplementedError(
                    f"Solvation model {solvation_model} is not implemented in nwchem."
                )
            if self.calc.parameters.solvation_radii not in ["", "default"]:
                raise UnimplementedError(
                    """Only default solvation radii are supported for nwchem.
                                           For manual specification of other solvation radii, use custom solvation radii"""
                )
            if self.calc.parameters.custom_solvation_radii != "":
                self.parse_custom_solvation_radii()
                self.radii_parameters = "\n".join(
                    [
                        f"{element} {self.solvation_radii[element]}"
                        for element in self.solvation_radii
                    ]
                )
                self.additional_block += (
                    f"parameters {self.calc.name}_sol.parameters \n"
                )
                warn(
                    f"Addtitional file {self.calc.name}_sol.parameters was generated. This file will be needed for calculation to run."
                )

    def parse_custom_solvation_radii(self):
        for radius in self.calc.parameters.custom_solvation_radii.split(";"):
            if radius.strip() == "":
                continue
            sradius = radius.split("=")
            if len(sradius) != 2:
                raise InvalidParameter(
                    f"Invalid custom solvation radius specification: '{radius}': must follow the pattern '<atom1>=<radius1>;...'"
                )

            element, rad = sradius
            if element not in LOWERCASE_ATOMIC_SYMBOLS:
                raise InvalidParameter(
                    f"Invalid element in custom solvation radius specifications: '{element}'"
                )

            _element = LOWERCASE_ATOMIC_SYMBOLS[element]  # Add the proper case back

            try:
                _rad = float(rad)
            except ValueError:
                raise InvalidParameter(
                    f"Invalid custom solvation radius for element {element}: '{rad}'"
                )
            self.solvation_radii[_element] = _rad

    def close_blocks(self):
        if self.calc.parameters.theory_level == "mp2" and self.method_block == "mp2\n":
            self.method_block = ""
        if (
            self.calc.parameters.theory_level == "ccsd"
            and self.method_block == "ccsd\n"
        ):
            self.method_block = ""
        if "end \n" not in self.method_block and self.method_block != "":
            self.method_block += " end \n"
        if "end \n" not in self.calculation_block and self.calculation_block != "":
            self.calculation_block += " end \n"
        if (
            "end \n" not in self.additional_block
            and self.clean(self.additional_block) != ""
        ):
            self.additional_block += "end \n"

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
        if self.radii_parameters != "":
            return f"---Input file---\n{self.input_file}\n\n---Solvation parameters file---\n{self.radii_parameters}"
        else:
            return self.input_file

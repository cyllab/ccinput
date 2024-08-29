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
    check_fragments,
    add_fragments_xyz,
    parse_specifications,
)
from ccinput.constants import CalcType, ATOMIC_NUMBER, LOWERCASE_ATOMIC_SYMBOLS
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation


class GaussianCalculation:
    TEMPLATE = """%chk={}.chk
    %nproc={}
    %mem={}MB
    #p {}

    {}

    {} {}
    {}
    {}
    """

    KEYWORDS = {
        CalcType.OPT: ["opt"],
        CalcType.CONSTR_OPT: ["opt"],
        CalcType.TS: ["opt"],
        CalcType.FREQ: ["freq"],
        CalcType.NMR: ["nmr"],
        CalcType.SP: ["sp"],
        CalcType.UVVIS: ["td"],
        CalcType.UVVIS_TDA: ["tda"],
        CalcType.OPTFREQ: ["opt", "freq"],
    }

    # Get a set of all unique calculation keywords
    KEYWORD_LIST = set([j for i in KEYWORDS.values() for j in i])

    # Number of processors
    # Amount of memory
    # Command line
    # Charge
    # Multiplicity
    # XYZ structure
    # Appendix

    def __init__(self, calc):
        self.calc = calc
        self.has_scan = False
        self.appendix = []
        self.command_line = ""

        self.commands = {}
        self.solvation_radii = {}

        self.confirmed_specifications = ""
        self.xyz_structure = ""
        self.input_file = ""
        if self.calc.type not in self.KEYWORDS:
            raise ImpossibleCalculation(
                f"Gaussian 16 does not support calculations of type {self.calc.type}"
            )

        self.handle_specifications()
        self.handle_command()
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

    def add_commands(self, commands):
        for cmd in commands:
            if cmd not in self.commands:
                self.commands[cmd] = []

    def handle_specifications(self):
        s = self.clean(self.calc.parameters.specifications.lower())

        parse_specifications(s, self.add_option)
        if self.calc.parameters.d3:
            self.add_option("EmpiricalDispersion", "GD3")
        elif self.calc.parameters.d3bj:
            self.add_option("EmpiricalDispersion", "GD3BJ")

    def filter_commands(self, commands):
        """Removes all commands which are not relevant"""
        keys = list(self.commands.keys())
        for key in keys:
            if key in self.KEYWORD_LIST and key not in commands:
                del self.commands[key]

    def handle_command(self):
        self.filter_commands(self.KEYWORDS[self.calc.type])
        self.add_commands(self.KEYWORDS[self.calc.type])

        if self.calc.type == CalcType.TS:
            self.add_options("opt", ["ts", "NoEigenTest", "CalcFC"])
        elif self.calc.type == CalcType.CONSTR_OPT:
            self.add_options("opt", ["modredundant"])

            xyz = get_npxyz(self.calc.xyz)
            gaussian_constraints = ""

            has_scan = False

            if len(self.calc.constraints) == 0:
                raise InvalidParameter("No constraint in constrained optimisation mode")

            for constr in self.calc.constraints:
                if constr.scan:
                    has_scan = True
                gaussian_constraints += constr.to_gaussian()

            self.has_scan = has_scan
            self.appendix.append(gaussian_constraints)

        method = get_method(self.calc.parameters.method, "gaussian")
        basis_set = get_basis_set(self.calc.parameters.basis_set, "gaussian")
        custom_basis_set = self.calc.parameters.custom_basis_sets

        if method == "":
            if self.calc.parameters.theory_level == "hf":
                method = "HF"
            else:
                raise InvalidParameter("No method")

        if basis_set != "":
            if len(custom_basis_set) == 0:
                if self.calc.parameters.density_fitting != "":
                    self.command_line += (
                        f"{method}/{basis_set}/{self.calc.parameters.density_fitting} "
                    )
                else:
                    self.command_line += f"{method}/{basis_set} "
            else:
                gen_keyword, to_append = self.parse_custom_basis_set(basis_set)
                if to_append.strip() == "":
                    self.command_line += f"{method}/{basis_set} "
                else:
                    self.appendix.append(to_append)
                    self.command_line += f"{method}/{gen_keyword} "
        else:
            self.command_line += f"{method} "
        # Counterpoise related commands processing
        if "counterpoise" in self.commands.keys():
            check_fragments(
                self.commands["counterpoise"][0], self.calc.fragments, self.calc.xyz
            )

    def parse_custom_basis_set(self, base_bs):
        custom_basis_sets = self.calc.parameters.custom_basis_sets
        to_append_gen = []
        to_append_ecp = []

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
        ecp = False
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
                    bs_keyword, fmt="gaussian94", elements=[el_num], header=False
                )
            except:
                # Some basis sets are built-in, but use different names as the BSE does (e.g., SDD)
                # In this case, just feed the user keyword in and hope it works.
                # The basis set string has been recognized by ccinput, so it should exist in the program.
                # ECP is added if Z > 18 (Ar)
                to_append_gen.append(f"{el} 0\n{bs_keyword}\n****\n")
                if el_num > 18:
                    ecp = True
                    to_append_ecp.append(f"{el} 0\n{bs_keyword}\n")

            else:
                if bs.find("-ECP") != -1:
                    ecp = True
                    sbs = bs.split("\n")
                    ecp_ind = -1
                    for ind, line in enumerate(sbs):
                        if sbs[ind].find("-ECP") != -1:
                            ecp_ind = ind
                            break
                    bs_gen = "\n".join(sbs[: ecp_ind - 2]) + "\n"
                    bs_ecp = "\n".join(sbs[ecp_ind - 2 :])
                    to_append_gen.append(bs_gen)
                    to_append_ecp.append(bs_ecp)
                else:
                    to_append_gen.append(bs)

        if len(custom_atoms) > 0:
            if ecp:
                gen_keyword = "GenECP"
            else:
                gen_keyword = "Gen"

            custom_bs = ""

            if len(normal_atoms) > 0:
                custom_bs += " ".join(normal_atoms) + " 0\n"
                custom_bs += base_bs + "\n"
                custom_bs += "****\n"

            custom_bs += "".join(to_append_gen) + "\n"
            custom_bs += "".join(to_append_ecp).replace("\n\n", "\n")

            return gen_keyword, custom_bs
        else:
            return self.calc.parameters.basis_set, ""

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        # If counterpoise correction is the option, modify xyz corresponding to fragments
        if self.calc.fragments != None:
            lines = add_fragments_xyz(lines, self.calc.fragments)
        self.xyz_structure = "".join(lines)

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

    def get_radii_appendix(self):
        radii_appendix = ""
        for el, rad in self.solvation_radii.items():
            radii_appendix += f"{el} {rad:.2f}\n"
        return radii_appendix

    def handle_solvation(self):
        if self.calc.parameters.solvent.lower() not in ["", "vacuum"]:
            solvent_keyword = get_solvent(
                self.calc.parameters.solvent,
                self.calc.parameters.software,
                solvation_model=self.calc.parameters.solvation_model,
            )

            model = self.calc.parameters.solvation_model
            radii_set = self.calc.parameters.solvation_radii
            custom_radii = self.calc.parameters.custom_solvation_radii

            DEFAULT_RADII_SETS = {
                "smd": ["", "default"],
                "pcm": ["", "uff", "default"],
                "cpcm": ["", "uff", "default"],
            }

            self.add_options("SCRF", [model.upper(), f"Solvent={solvent_keyword}"])

            if radii_set not in DEFAULT_RADII_SETS[model] or custom_radii != "":
                self.add_option("SCRF", "Read")
                self.parse_custom_solvation_radii()

            solv_appendix = ""

            if model == "smd":
                if radii_set == "smd18":
                    # Refined solvation radii
                    # E. Engelage, N. Schulz, F. Heinen, S. M. Huber, D. G. Truhlar,
                    # C. J. Cramer, Chem. Eur. J. 2018, 24, 15983-15987.
                    if "Br" not in self.solvation_radii:
                        self.solvation_radii["Br"] = 2.60
                    if "I" not in self.solvation_radii:
                        self.solvation_radii["I"] = 2.74
            elif model in ["pcm", "cpcm"]:
                if radii_set not in DEFAULT_RADII_SETS[model]:
                    solv_appendix += f"Radii={radii_set}\n"
            else:
                raise InvalidParameter(
                    f"Invalid solvation method for Gaussian: '{self.calc.parameters.solvation_model}'"
                )

            if len(self.solvation_radii) > 0:
                solv_appendix += "modifysph\n\n"
                solv_appendix += self.get_radii_appendix()

            self.appendix.append(solv_appendix)

    def create_input_file(self):
        for cmd, options in self.commands.items():
            option_str = ", ".join(options)
            if option_str != "":
                cmd_formatted = f"{cmd}({option_str}) "
            else:
                cmd_formatted = f"{cmd} "

            if cmd.lower() != "scrf":
                spec_options = []
                for o in options:
                    if o.strip() == "":
                        continue
                    if self.calc.type == CalcType.TS and o in [
                        "ts",
                        "NoEigenTest",
                        "CalcFC",
                    ]:
                        continue
                    elif self.calc.type == CalcType.CONSTR_OPT and o == "modredundant":
                        continue
                    spec_options.append(o)
                if len(spec_options) > 0:
                    spec_option_str = ", ".join(spec_options)

                    confirmed_spec = f"{cmd}({spec_option_str}) "
                    if cmd in self.KEYWORD_LIST:
                        self.confirmed_specifications = (
                            confirmed_spec + self.confirmed_specifications
                        )
                    else:
                        self.confirmed_specifications += confirmed_spec

            # This ensures that the command line follows this pattern:
            # CMD1 <CMD2> METHOD/BASIS_SET <ADDITIONAL_OPTION1> ...
            if cmd in self.KEYWORD_LIST:
                self.command_line = cmd_formatted + self.command_line
            else:
                self.command_line = self.command_line + cmd_formatted

        self.confirmed_specifications = self.confirmed_specifications.strip()

        raw = self.TEMPLATE.format(
            self.calc.name,
            self.calc.nproc,
            self.calc.mem,
            self.command_line.strip(),
            self.calc.header,
            self.calc.charge,
            self.calc.multiplicity,
            self.xyz_structure,
            "\n".join(self.appendix),
        )
        self.input_file = "\n".join([i.strip() for i in raw.split("\n")]).replace(
            "\n\n\n", "\n\n"
        )

    @property
    def output(self):
        return self.input_file

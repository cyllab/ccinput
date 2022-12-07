import numpy as np
import copy

from ccinput.utilities import (
    get_method,
    get_solvent,
    get_basis_set,
    clean_xyz,
)
from ccinput.constants import CalcType, ATOMIC_NUMBER, LOWERCASE_ATOMIC_SYMBOLS
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation


class PysisDriver:

    SUPPORTED_PACKAGES = ["xtb"]

    SUPPORTED_KEYWORDS = {
        CalcType.TS: "tsopt",
    }

    DEFAULT_SECTION_PARAMS = {
        "tsopt": {
            "type": "rsirfo",
            "do_hess": False,
            "hessian_recalc": 5,
        },
        "calc": {
            "type": "xtb",
            "pal": 1,
            "charge": 0,
            "mult": 1,
        },
        "geom": {
            "type": "dlc",
            "fn": "calc.xyz",
        },
    }

    """
        Keys of self.sections are sections of the input file (e.g., geom, calc, tsopt).
        The corresponding values are dictionaries containing the key-value pairs within that section.

        Specifications/parameters are given through the specifications field in a pseudo-Gaussian format:
        $ [...] --specifications "geom(type=redund) tsopt(do_hess=False)"
    """

    def __init__(self, calc):
        self.calc = calc
        self.command = "pysis calc.inp"

        self.input_file = ""

        if self.calc.type not in self.SUPPORTED_KEYWORDS:
            raise ImpossibleCalculation(
                f"Pysis does not support calculations of type {self.calc.type} or its support is not implemented in ccinput"
            )

        self.sections = {}

        self.add_section("geom")
        self.handle_main_parameters()

        self.handle_specifications()

        self.create_input_file()

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. _"
        )
        return "".join([c for c in s if c in WHITELIST])

    def add_option(self, section, opt, val):
        if section not in self.sections:
            self.add_section(section)

        self.sections[section][opt] = val

    def add_options(self, sec, options):
        for opt, val in options:
            self.add_option(sec, opt, val)

    def add_section(self, section):
        if section not in self.sections:
            if section in self.DEFAULT_SECTION_PARAMS:
                self.sections[section] = copy.deepcopy(
                    self.DEFAULT_SECTION_PARAMS[section]
                )
            else:
                self.sections[section] = {}
        else:
            logger.warning(f"Duplicate specifications for section {section} were given")

    def handle_specifications(self):

        s = self.clean(self.calc.parameters.specifications.lower())

        # Duplicate code from ccinput/packages/Gaussian.py
        # TODO: refactor
        if s.count("(") != s.count(")"):
            raise InvalidParameter("Invalid specifications: parenthesis not matching")

        _specifications = ""
        remove = False
        for c in s:
            if c == " " and remove:
                continue
            _specifications += c
            if c == "(":
                remove = True
            elif c == ")":
                remove = False

        for spec in _specifications.split(" "):
            if spec.strip() == "":
                continue
            if spec.find("(") != -1:
                key, options = spec.split("(")
                options = options.replace(")", "")
                for option in options.split(","):
                    if option.find("=") == -1:
                        raise InvalidParameter(
                            f"Invalid specification {option}: no key/value delimiter (equals sign)"
                        )
                    if option.strip() != "":
                        param, val = option.split("=")
                        self.add_option(key, param, val)
            else:
                raise InvalidParameter(
                    f"Invalid specification {spec}: should be in the form 'section(key1=value1, ...)'"
                )

    def handle_main_parameters(self):

        if self.calc.parameters.software not in self.SUPPORTED_PACKAGES:
            raise InvalidParameter(
                f"{self.calc.parameters} is not currently supported with the pysis driver"
            )

        self.add_option("calc", "type", self.calc.parameters.software)
        self.add_option("calc", "charge", self.calc.charge)
        self.add_option("calc", "mult", self.calc.multiplicity)
        self.add_option("calc", "pal", self.calc.nproc)

        if self.calc.type == CalcType.TS:
            self.add_section("tsopt")
        else:
            raise InvalidParameter(
                f"Calculation type {self.calc.type} is not currently implemented for pysis in ccinput"
            )

        """
        For other computation packages:
            method
            basis set
            custom basis set
            constraints
        """

        if self.calc.parameters.solvent != "":
            if self.calc.parameters.software == "xtb":
                if self.calc.parameters.solvation_model not in ["gbsa", "alpb"]:
                    raise InvalidParameter(
                        f"Invalid solvation method for xtb: {self.calc.parameters.solvation_model}"
                    )

                try:
                    solvent_keyword = get_solvent(
                        self.calc.parameters.solvent, self.calc.parameters.software
                    )
                except KeyError:
                    raise InvalidParameter("Invalid solvent")

                self.add_option(
                    "calc", self.calc.parameters.solvation_model, solvent_keyword
                )
            else:
                raise InvalidParameter(
                    f"Cannot apply solvation for {self.calc.parameters.software}"
                )

    def create_input_file(self):
        self.input_file = ""
        for section, content in sorted(self.sections.items(), key=lambda i: i[0]):
            self.input_file += f"{section}:\n"
            for key, val in sorted(content.items(), key=lambda i: i[0]):
                self.input_file += f"    {key}: {val}\n"

    @property
    def output(self):
        return self.input_file

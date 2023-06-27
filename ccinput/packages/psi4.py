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
    warn,
)
from ccinput.calculation import Calculation
from ccinput.constants import CalcType, ATOMIC_NUMBER, LOWERCASE_ATOMIC_SYMBOLS
from ccinput.exceptions import UnimplementedError, InvalidParameter


class Psi4Calculation:
    TEMPLATE = """#{header}

{memory_line}

molecule {{
{charge} {multiplicity}  
{coordinates}}}

{command_line}"""

    CALC_TYPES = {
        CalcType.SP: ["sp", "energy"],
        CalcType.OPT: ["opt", "optimize"],
        # CalcType.CONSTR_OPT,
        CalcType.FREQ: ["freq", "frequency", "frequencies"],
        # CalcType.TS,
        CalcType.OPTFREQ: ["Opt+Freq Calculation", "optfreq", "opt+freq"],
    }

    def __init__(self, calc: Calculation):
        self.calc = calc
        self.has_scan = False
        self.appendix = []
        self.command_line = ""
        self.commands = {}
        self.solvation_radii = {}
        self.confirmed_specifications = ""
        self.xyz_structure = ""
        self.input_file = ""

        if self.calc.type not in self.CALC_TYPES:
            raise UnimplementedError(
                f"Calculation Type {self.calc.type} not implemented yet for psi4"
            )

        self.type_method = ""

        # self.handle_specifications()
        self.handle_command()
        self.handle_xyz()

        self.handle_solvation()

        self.create_input_file()

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,._ "
        )
        return "".join([c for c in s if c in WHITELIST])

    def add_option(self, key, option):
        _option = option.strip()
        if key in self.commands:
            if _option not in self.commands[key]:
                self.commands[key].append(_option)  # remove.append and add =
        else:
            self.commands[key] = [_option]  # remove the brakets of the list

    def add_options(self, key, options):
        for option in options:
            self.add_option(key, option)

    def handle_specifications(self, _specifications):
        """Returns a string of specifications separated by comma."""

        # _specifications = self.clean(self.calc.parameters.specifications)

        if _specifications != "":
            sspecs = _specifications.split()
            if "," in str(sspecs):
                specs_no_comma = [s.replace(",", "") for s in sspecs]
                specifications_list = ", ".join(specs_no_comma)
            else:
                specifications_list = ", ".join(sspecs)
            return specifications_list

    def extract_string_between_delimiter(self, s, start_delimiter, end_delimiter):
        start = start_delimiter
        end = end_delimiter
        return s[s.find(start) + len(start) : s.rfind(end)]

    def extract_string_before_delimeter(self, s, delimiter):
        end = delimiter
        s = s[0 : s.rfind(end)]
        return s

    def check_multispec_syntax(self, multi_specs):
        """Performs simple syntax check for specifications of combined calculations.
        Returns: dict {True: 'wrong specification'}
        """
        # simple syntax check for combined calculations, for now is for optfreq. It can be more sophisitcated.
        # checks that multi specifications have the syntax "calctype(specifications)". If wrong it appends True to an empty list and the worng value to another empty list.

        cheked_specs = []
        wrong_spec_list = []

        if multi_specs.find("(") == -1:
            cheked_specs.append(True)
            wrong_spec_list.append(multi_specs)
        if multi_specs.count("(") != multi_specs.count(")"):
            cheked_specs.append(True)
            wrong_spec_list.append(multi_specs)
        if multi_specs.find("(") != -1:
            _specifications = multi_specs.strip().split(")")
            for spec in _specifications:
                if spec != "":
                    if spec.find("(") == -1:
                        cheked_specs.append(True)
                        wrong_spec_list.append(spec)
                    # checks if the string calctype is less than 5 chars. This is the case for opt or freq.
                    # This can be more sophisticated for other type of combined calculation types.
                    elif (
                        spec.find("(") != -1
                        and len(self.extract_string_before_delimeter(spec, "(")) > 5
                    ):
                        cheked_specs.append(True)
                        wrong_spec_list.append(spec)
                    elif self.calc.type == CalcType.OPTFREQ:
                        calctype_string = self.extract_string_before_delimeter(
                            spec, "("
                        ).strip()
                        if calctype_string not in set(
                            ("opt", "freq")
                        ):  # could this be extended to calctype synonyms?
                            cheked_specs.append(True)
                            wrong_spec_list.append(calctype_string)
        wrong_specs_bool_dict = dict(zip(cheked_specs, wrong_spec_list))
        return wrong_specs_bool_dict

    def convert_multispecs_to_dictionary(self):
        """Returns a dictionary of calculation call and specifications when specifications are necessary for combined claculations like OPTFREQ (optimization followed by frequency analysis).
        If specifications are needed for both calculation calls use: "calctype1(spec1 spec2 ...) calctype2(spec1 spec2 ...)" otherwhise use "calctype(specs1 ...)"

        ### Example:
        #### Case for combined calculations as optfreq calculation.
        - Specifications for both calls : "opt(spec1 spec2) freq(spec1 spec2)"
        - Specifications for one call : "opt(spec1 spec2)"
        - Returns: {'opt':'spec1 spec2', 'freq': 'spec1 spec2'} or {'opt': 'spec1 spec2'}
        """
        _specifications = self.clean(self.calc.parameters.specifications)
        keys = []
        values = []

        # If there are specifications, it verifies if their syntax is OK.
        checked_spec_value = self.check_multispec_syntax(_specifications)
        if _specifications != "" and any(checked_spec_value.keys()) == True:
            raise InvalidParameter(
                (
                    f"Invalid specification: '{checked_spec_value[True]}'. Perhaps you forgot to write the calctype name?\nAre you using the syntax 'calctype1(spec1 spec2 ...) calctype2(spec1 spec2 ...) or just calctype(specs)', where specs are space separeted?"
                )
            )

        if _specifications.find("(") != -1:
            multi_specifications = _specifications.split(")")
            for specs in multi_specifications:
                if specs != "":
                    _key = self.extract_string_before_delimeter(specs, "(")
                    keys.append(_key.strip())
                    _value = self.extract_string_between_delimiter(specs, "(", "")
                    values.append(_value.strip())

        specs_dictionary = dict(zip(keys, values))
        return specs_dictionary

    def create_command_line_string(self, calc_type, method, basis_set, specifications):
        """ "Returns a string containing input parameters as per calculation type"""

        if isinstance(specifications, str) and str != "":
            command_line_string = (
                f"{calc_type}('{method}/{basis_set}', {specifications})"
            )
        if specifications is None:
            command_line_string = f"{calc_type}('{method}/{basis_set}')"
        return command_line_string

    def handle_command(self):
        """ "Handles parameters to generate the input file commands as per calcuation type.
        Returns the {comand_line} string"""

        # method label, can be used to specify quantum method or functional.
        method = get_method(self.calc.parameters.method, "psi4")
        basis_set = get_basis_set(self.calc.parameters.basis_set, "psi4")
        _specifications = self.clean(self.calc.parameters.specifications)
        sspecs = self.handle_specifications(_specifications)

        if self.calc.type == CalcType.SP:
            calc_type = "energy"
            cmls_sp = self.create_command_line_string(
                calc_type, method, basis_set, specifications=sspecs
            )
            self.command_line += cmls_sp
        elif self.calc.type == CalcType.OPT:
            calc_type = "optimize"
            cmls_opt = self.create_command_line_string(
                calc_type, method, basis_set, specifications=sspecs
            )
            self.command_line += cmls_opt
        elif self.calc.type == CalcType.FREQ:
            calc_type = "frequencies"
            cmls_freqs = self.create_command_line_string(
                calc_type, method, basis_set, sspecs
            )
            self.command_line += cmls_freqs
        elif self.calc.type == CalcType.OPTFREQ:
            multi_specs = self.convert_multispecs_to_dictionary()
            calc_type = {"freq": "frequencies", "opt": "optimize"}
            opt_specs = None
            freq_specs = None

            if "opt" in multi_specs:
                opt_specs = str("")
                opt_specs += self.handle_specifications(multi_specs["opt"])
            if "freq" in multi_specs:
                freq_specs = str("")
                freq_specs += self.handle_specifications(multi_specs["freq"])
            elif "opt" not in multi_specs:
                opt_specs = None
            elif "freq" not in multi_specs:
                freq_specs = None
            else:
                raise InvalidParameter(
                    f"Keys in {multi_specs.keys()} not valid for {self.calc.type} type of calculation"
                )

            cmls_opt = self.create_command_line_string(
                calc_type["opt"], method, basis_set, specifications=opt_specs
            )
            cmls_freqs = self.create_command_line_string(
                calc_type["freq"], method, basis_set, specifications=freq_specs
            )
            self.command_line += f"{cmls_opt}\n{cmls_freqs}"
        else:
            raise UnimplementedError(
                f"Calculation Type {self.calc.type} not implemented yet for psi4"
            )

    def handle_memory(self):
        """Returns a str for the input memory line"""
        memory_line = ""
        mem = self.calc.mem
        self.commands["mem"] = mem
        if bool(self.commands["mem"]) == False:
            memory_line
        elif bool(self.commands["mem"]) == True:
            memory_line += f"memory {mem} mb"
        return memory_line

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def handle_solvation(self):
        return

    def create_input_file(self) -> None:
        self.input_file = self.TEMPLATE.format(
            header=self.calc.header,
            memory_line=self.handle_memory(),
            charge=self.calc.charge,
            multiplicity=self.calc.multiplicity,
            coordinates=self.xyz_structure,
            command_line=self.command_line,
        )
        return

    @property
    def output(self):
        return self.input_file

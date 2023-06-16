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

        #self.handle_specifications()
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
                self.commands[key].append(_option) #remove.append and add =
        else:
            self.commands[key] = [_option] # remove the brakets of the list

    def add_options(self, key, options):
        for option in options:
            self.add_option(key, option)


    def handle_specifications(self, _specifications):
        """ Returns a string of specifications separated by comma.
        """

        #_specifications = self.clean(self.calc.parameters.specifications)
        
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
        return s[s.find(start)+len(start):s.rfind(end)]

    def extract_string_before_delimeter(self, s, delimiter):
        end = delimiter
        s = s[0:s.rfind(end)]
        return s
    
    def convert_multispecs_to_dictionary(self):
        """Returns a dictionary of specifications when specifications are defined for more than one keyword such as more than one caulation type.
        For example: specifications: "opt(spec1 spec2) freq(spec1 spec2)
        Retruns {'opt':'spec1 spec2', 'freq': 'spec1 spec2'}
        """
        _specifications = self.clean(self.calc.parameters.specifications)
        keys = []
        values = []
        
        if _specifications.find("(") != -1:
            multi_specifications = _specifications.split(")")
    
            for specs in multi_specifications:
                if specs !="":
                    _key = self.extract_string_before_delimeter(specs, "(") # strip removes undesirable chars default = " "
                    keys.append(_key.strip())
                    _value = self.extract_string_between_delimiter(specs,"(","")
                    values.append(_value.strip())
            
        specs_dictionary = dict(zip(keys, values))
        return specs_dictionary
    
        
    def create_command_line_string(self, calc_type, method, basis_set, specifications):
        """"Returns a string containing input parameters as per calculation type"""

        if isinstance(specifications, str) and str != "":
            
            # if self.calc.type == CalcType.OPTFREQ:
            #     warn("Warning: Specifications for opt+freq only implemented for the optimization call. Not implemented for the frequencies call yet.")
            #     #command_line_string = f"{calc_type}('{method}/{basis_set}', {specifications})"
            #     command_line_string = f"{calc_type}('{method}/{basis_set}', {specifications})"
            # else:
            #     command_line_string = f"{calc_type}('{method}/{basis_set}', {specifications})"
            command_line_string = f"{calc_type}('{method}/{basis_set}', {specifications})"
        if specifications is None:
            command_line_string = f"{calc_type}('{method}/{basis_set}')"
        return command_line_string

    def handle_command(self):
        """"Handles parameters to generate the input file commands as per calcuation type. 
        Returns the {comand_line} string """
        
        # method label, can be used to specify quantum method or functional.
        method = get_method(self.calc.parameters.method, "psi4")
        basis_set = get_basis_set(self.calc.parameters.basis_set, "psi4")
        _specifications = self.clean(self.calc.parameters.specifications)
        sspecs = self.handle_specifications(_specifications)
        
        if self.calc.type == CalcType.SP:
            calc_type = "energy"
            cmls_sp = self.create_command_line_string(calc_type, method, basis_set, specifications=sspecs)
            self.command_line += cmls_sp
        elif self.calc.type == CalcType.OPT:
            calc_type = "optimize"
            cmls_opt = self.create_command_line_string(calc_type, method, basis_set, specifications=sspecs)
            self.command_line += cmls_opt
        elif self.calc.type == CalcType.FREQ:
            calc_type = "frequencies"
            cmls_freqs = self.create_command_line_string(calc_type, method, basis_set,sspecs)
            self.command_line += cmls_freqs
        elif self.calc.type == CalcType.OPTFREQ:
            multi_specs = self.convert_multispecs_to_dictionary()
            calc_type= {"freq":"frequencies", "opt":"optimize"}
            opt_specs = None
            freq_specs = None
            
            if "opt" in multi_specs:
                opt_specs = str("")
                opt_specs += self.handle_specifications(multi_specs['opt'])
            if "freq" in multi_specs:
                freq_specs = str("")
                freq_specs += self.handle_specifications(multi_specs['freq'])
            elif "opt" not in multi_specs:
                opt_specs = None
            elif "freq" not in multi_specs:
                freq_specs = None
            else:
                raise InvalidParameter( f"Keys in {multi_specs.keys()} not valid for {self.calc.type} type of calculation")
            
            cmls_opt = self.create_command_line_string(calc_type["opt"], method, basis_set, specifications=opt_specs)
            cmls_freqs = self.create_command_line_string(calc_type["freq"], method, basis_set, specifications=freq_specs)
            self.command_line += f"{cmls_opt}\n{cmls_freqs}"
        else:
            raise UnimplementedError(
                f"Calculation Type {self.calc.type} not implemented yet for psi4"
            )
    
    
    def handle_memory(self):
        """Retruns a str for the input memory line"""
        memory_line = ""
        mem = self.calc.mem
        self.commands['mem'] = mem
        if bool(self.commands['mem']) == False:
            memory_line
        elif bool(self.commands['mem']) == True:
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

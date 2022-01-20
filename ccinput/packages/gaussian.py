import basis_set_exchange
import numpy as np

from ccinput.utilities import get_method, get_solvent, get_basis_set, clean_xyz, \
                              get_distance, get_angle, get_dihedral, get_npxyz
from ccinput.constants import CalcType, ATOMIC_NUMBER
from ccinput.exceptions import InvalidParameter

class GaussianCalculation:

    # Refined solvation radii
    # E. Engelage, N. Schulz, F. Heinen, S. M. Huber, D. G. Truhlar,
    # C. J. Cramer, Chem. Eur. J. 2018, 24, 15983–15987.
    SMD18_APPENDIX = """modifysph

    Br 2.60
    I 2.74
    """

    TEMPLATE = """%chk={}.chk
    %nproc={}
    %mem={}MB
    #p {} {}

    {}

    {} {}
    {}
    {}
    """

    KEYWORDS = {
                CalcType.OPT: 'opt',
                CalcType.CONSTR_OPT: 'opt',
                CalcType.TS: 'opt',
                CalcType.FREQ: 'freq',
                CalcType.NMR: 'nmr',
                CalcType.SP: 'sp',
                CalcType.UVVIS: 'td',
            }

    #Number of processors
    #Amount of memory
    #Command line
    #Additional commands
    #Charge
    #Multiplicity
    #XYZ structure
    #Appendix

    def __init__(self, calc):
        self.calc = calc
        self.has_scan = False
        self.appendix = []
        self.command_line = ""
        self.additional_commands = []
        self.command_specifications = []
        self.confirmed_specifications = ""
        self.xyz_structure = ""
        self.input_file = ""

        self.handle_specifications()
        self.handle_command()
        self.handle_xyz()

        self.handle_solvation()

        self.create_input_file()

    def clean(self, s):
        WHITELIST = set("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. ")
        return ''.join([c for c in s if c in WHITELIST])

    def handle_specifications(self):
        specs = {}
        def add_spec(key, option):
            if option == '':
                return
            if key in specs:
                if option not in specs[key]:
                    specs[key].append(option)
            else:
                specs[key] = [option]

        s = self.clean(self.calc.parameters.specifications.lower())

        #Could be more sophisticated to catch other incorrect specifications
        if s.count('(') != s.count(')'):
            raise InvalidParameter("Invalid specifications: parenthesis not matching")

        _specifications = ""
        remove = False
        for c in s:
            if c == ' ' and remove:
                continue
            _specifications += c
            if c == '(':
                remove = True
            elif c == ')':
                remove = False

        for spec in _specifications.split(' '):
            if spec.strip() == '':
                continue
            if spec.find("(") != -1:
                key, options = spec.split('(')
                options = options.replace(')', '')
                if key == self.KEYWORDS[self.calc.type]:
                    for option in options.split(','):
                        if option not in self.command_specifications:
                            self.command_specifications.append(option.strip())
                else:
                    if key in self.KEYWORDS.values():
                        continue#Invalid specification
                    for option in options.split(','):
                        add_spec(key, option.strip())
            else:
                if spec not in self.additional_commands:
                    self.additional_commands.append(spec)

        for spec, option in specs.items():
            specs_str = ','.join(option)
            spec_formatted = f'{spec}({specs_str}) '
            self.additional_commands.append(spec_formatted)

    def handle_command(self):
        cmd = ""
        base_specs = []
        if self.calc.type == CalcType.NMR:
            cmd = "nmr"
        elif self.calc.type == CalcType.UVVIS:
            cmd = "td"
        elif self.calc.type == CalcType.OPT:
            cmd = "opt"
        elif self.calc.type == CalcType.TS:
            cmd = "opt"
            base_specs = ['ts', 'NoEigenTest', 'CalcFC']
        elif self.calc.type == CalcType.FREQ:
            cmd = "freq"
        elif self.calc.type == CalcType.CONSTR_OPT:
            cmd = "opt"
            base_specs = ['modredundant']

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
        elif self.calc.type == CalcType.SP:
            cmd = "sp"

        if len(self.command_specifications) > 0:
            self.confirmed_specifications += f"{cmd}({', '.join(self.command_specifications)}) "
        full_specs = base_specs + self.command_specifications

        if len(full_specs) == 0:
            cmd_full = cmd + ' '
        else:
            cmd_specifications_str = f"({', '.join(full_specs)})"
            cmd_full = f"{cmd}{cmd_specifications_str} "

        self.command_line += cmd_full

        method = get_method(self.calc.parameters.method, "gaussian")
        basis_set = get_basis_set(self.calc.parameters.basis_set, "gaussian")
        custom_basis_set = self.calc.parameters.custom_basis_sets

        if method == "":
            if self.calc.parameters.theory_level == "hf":
                method = "HF"
            else:
                raise InvalidParameter("No method")

        if basis_set != "":
            if custom_basis_set == "":
                if self.calc.parameters.density_fitting != '':
                    self.command_line += f"{method}/{basis_set}/{self.calc.parameters.density_fitting} "
                else:
                    self.command_line += f"{method}/{basis_set} "
            else:
                gen_keyword, to_append = self.parse_custom_basis_set(basis_set)
                if to_append.strip() == '':
                    self.command_line += f"{method}/{basis_set} "
                else:
                    self.appendix.append(to_append)
                    self.command_line += f"{method}/{gen_keyword} "
        else:
            self.command_line += f"{method} "

    def parse_custom_basis_set(self, base_bs):
        custom_basis_set = self.calc.parameters.custom_basis_sets
        entries = [i.strip() for i in custom_basis_set.split(';') if i.strip() != ""]
        to_append_gen = []
        to_append_ecp = []

        custom_atoms_requested = []
        for entry in entries:
            sentry = entry.split('=')

            if len(sentry) != 2:
                raise InvalidParameter(f"Invalid custom basis set string: '{entry}'")

            el, bs_keyword = sentry
            custom_atoms_requested.append(el)


        unique_atoms = []
        normal_atoms = []
        for line in self.calc.xyz.split('\n'):
            if line.strip() == "":
                continue
            a, *_ = line.split()
            if a not in unique_atoms:
                unique_atoms.append(a)
                if a not in normal_atoms and a not in custom_atoms_requested:
                    normal_atoms.append(a)

        custom_atoms = []
        ecp = False
        for entry in entries:
            sentry = entry.split('=')

            el, bs_keyword = sentry

            if el not in unique_atoms:
                continue

            custom_atoms.append(el)
            try:
                el_num = ATOMIC_NUMBER[el]
            except KeyError:
                raise InvalidParameter(f"Invalid atom in custom basis set string: '{el}'")

            bs = basis_set_exchange.get_basis(bs_keyword, fmt='gaussian94', elements=[el_num], header=False)
            if bs.find('-ECP') != -1:
                ecp = True
                sbs = bs.split('\n')
                ecp_ind = -1
                for ind, line in enumerate(sbs):
                    if sbs[ind].find("-ECP") != -1:
                        ecp_ind = ind
                        break
                bs_gen = '\n'.join(sbs[:ecp_ind-2]) + '\n'
                bs_ecp = '\n'.join(sbs[ecp_ind-2:])
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
                custom_bs += ' '.join(normal_atoms) + ' 0\n'
                custom_bs += base_bs + '\n'
                custom_bs += '****\n'

            custom_bs += ''.join(to_append_gen)
            custom_bs += ''.join(to_append_ecp).replace('\n\n', '\n')

            return gen_keyword, custom_bs
        else:
            return self.calc.parameters.basis_set, ''

    def handle_xyz(self):
        lines = [i + '\n' for i in clean_xyz(self.calc.xyz).split('\n') if i != '' ]
        self.xyz_structure = ''.join(lines)

    def handle_solvation(self):
        if self.calc.parameters.solvent.lower() not in ["", "vacuum"]:
            solvent_keyword = get_solvent(self.calc.parameters.solvent, self.calc.parameters.software, solvation_model=self.calc.parameters.solvation_model)
            if self.calc.parameters.solvation_model == "smd":
                if self.calc.parameters.solvation_radii in ["",  "default"]:
                    self.command_line += f"SCRF(SMD, Solvent={solvent_keyword}) "
                else:
                    self.command_line += f"SCRF(SMD, Solvent={solvent_keyword}, Read) "
                    self.appendix.append(self.SMD18_APPENDIX)
            elif self.calc.parameters.solvation_model == "pcm":
                if self.calc.parameters.solvation_radii in ["uff", ""]:
                    self.command_line += f"SCRF(PCM, Solvent={solvent_keyword}) "
                else:
                    self.command_line += f"SCRF(PCM, Solvent={solvent_keyword}, Read) "
                    self.appendix.append(f"Radii={self.calc.parameters.solvation_radii}\n")
            elif self.calc.parameters.solvation_model == "cpcm":
                if self.calc.parameters.solvation_radii in ["uff", ""]:
                    self.command_line += f"SCRF(CPCM, Solvent={solvent_keyword}) "
                else:
                    self.command_line += f"SCRF(CPCM, Solvent={solvent_keyword}, Read) "
                    self.appendix.append(f"Radii={self.calc.parameters.solvation_radii}\n")
            else:
                raise InvalidParameter(f"Invalid solvation method for Gaussian: '{self.calc.parameters.solvation_model}'")

    def create_input_file(self):
        additional_commands = " ".join([i.strip() for i in self.additional_commands]).strip()
        self.confirmed_specifications += additional_commands
        raw = self.TEMPLATE.format(self.calc.name, self.calc.nproc, self.calc.mem, self.command_line.strip(), additional_commands, self.calc.header, self.calc.charge, self.calc.multiplicity, self.xyz_structure, '\n'.join(self.appendix))
        self.input_file = '\n'.join([i.strip() for i in raw.split('\n')]).replace('\n\n\n', '\n\n')


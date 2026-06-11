from ccinput.constants import (
    CalcType,
    ATOMIC_NUMBER,
    LOWERCASE_ATOMIC_SYMBOLS,
    SOFTWARE_BASIS_SETS,
)
from ccinput.utilities import (
    get_method,
    get_solvent,
    clean_xyz,
    warn,
    parse_specifications,
)
from ccinput.exceptions import (
    InvalidParameter,
    UnimplementedError,
    ImpossibleCalculation,
)


class PySCFCalculation:
    calc = None

    has_scan = False
    pal = 0
    mem = 0
    blocks = []

    TEMPLATE = """!{}
    *xyz {} {}
    {}*
    {}"""
    # Command Line
    # Charge
    # Multiplicity
    # XYZ structure
    # Options blocks

    command_line = ""
    xyz_structure = ""
    block_lines = ""

    input_file = ""

    CALC_TYPES = [
        CalcType.SP,
        CalcType.OPT,
        CalcType.CONSTR_OPT,
        CalcType.FREQ,
        CalcType.TS,
        CalcType.OPTFREQ,
    ]

    def __init__(self, calc):
        self.calc = calc
        self.has_scan = False
        self.pal = 0
        self.mem = 0

        ###

        self.mol_specifications = {
            "verbose": 5,
        }

        self.command_line = ""
        self.output_elements = {}
        self.mf_lines = []

        self.additional_commands = ""
        self.xyz_structure = ""
        self.input_file = ""

        self.specifications = {}
        self.specifications_lines = []
        self.import_lines = []

        self.solvation_radii = {}
        self.aux_basis_set = ""

        if self.calc.type not in self.CALC_TYPES:
            raise ImpossibleCalculation(
                f"PySCF does not support calculations of type {self.calc.type}"
            )

        self.handle_specifications()

        self.handle_command()
        self.handle_xyz()

        self.handle_pal_mem()
        self.parse_custom_solvation_radii()
        self.handle_solvation()

        self.create_input_file()

    @property
    def confirmed_specifications(self):
        """Returns the effective additional commands; saved in the Parameters object"""
        return self.additional_commands

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. {}_[]"
        )
        return "".join([c for c in s if c in WHITELIST])

    def add_to_block(self, block, lines):
        if block in self.blocks:
            self.blocks[block] += lines
        else:
            self.blocks[block] = lines

    def add_option(self, key, option):
        if option == "":
            # Maybe not rigorously correct, but should do the expected thing
            if key[-2:] == "/c":
                self.aux_basis_set = key[:-2]
            elif key[-3:] == "/jk":
                self.aux_basis_set = key[:-3]
            elif key[-2:] == "/j":
                self.aux_basis_set = key[:-2]
            elif key not in self.specifications_list:
                self.specifications_list.append(key)
        else:
            try:
                option = int(option)
            except ValueError:
                pass
            self.mol_specifications[key] = option

    def handle_specifications(self):
        _specifications = (
            self.clean(self.calc.parameters.specifications).lower().strip()
        )

        parse_specifications(_specifications, self.add_option, condense=False)

        if self.calc.parameters.d3:
            raise Exception("Not implemented")
        elif self.calc.parameters.d3bj:
            raise Exception("Not implemented")

    def handle_command(self):
        if self.calc.type == CalcType.OPT:
            self.command_line = "OPT "
            # imports
        elif self.calc.type == CalcType.OPTFREQ:
            self.command_line = "OPT FREQ "
            # imports
        elif self.calc.type == CalcType.TS:
            self.command_line = """
                mf.kernel()
                params = {'transition': True, 'trust': 0.02, 'tmax': 0.06, "hessian": "file:hessian.txt"}
                mol_ts = mf.Gradients().optimizer(solver='geomeTRIC').kernel(params)
                """.strip()
            self.output_elements = {
                "xyz": "mol_ts.tostring()",
            }
        elif self.calc.type == CalcType.FREQ:
            self.command_line = """
                mf.kernel()
                ohess = mf.Hessian()
                hess = ohess.kernel()
                res_vib = harmonic_analysis(mf.mol, hess, imaginary_freq=False)
                res_thermo = thermo(mf, res_vib["freq_au"])
                """.strip()
            self.import_lines.append(
                "from pyscf.hessian.thermo import thermo, harmonic_analysis"
            )
            self.output_elements = {
                "SCF": "mf.e_tot",
                "H": 'res_thermo["H_tot"][0]',
                "G": 'res_thermo["G_tot"][0]',
                "hessian": "hess.tolist()",
                "freqs": 'res_vib["freq_wavenumber"].tolist()',
                "modes": 'res_vib["norm_mode"].tolist()',
            }
        elif self.calc.type == CalcType.CONSTR_OPT:
            ###
            # imports
            if len(self.calc.constraints) == 0:
                raise InvalidParameter("No constraints for constrained optimisation")

            scans = []
            freeze = []

            for constr in self.calc.constraints:
                if constr.scan:
                    scans.append(constr.to_orca())
                else:
                    freeze.append(constr.to_orca())

            if len(scans) > 0:
                scan_block = """Scan
                {}
                end
                """
                self.add_to_block(
                    "geom", scan_block.format("".join(scans).strip()).split("\n")
                )

            if len(freeze) > 0:
                freeze_block = """Constraints
                {}
                end
                """
                self.add_to_block(
                    "geom", freeze_block.format("".join(freeze).strip()).split("\n")
                )
        elif self.calc.type == CalcType.SP:
            self.command_line = "mf.kernel()"
            self.output_elements = {
                "SCF": "mf.e_tot",
            }

        self.method = get_method(self.calc.parameters.method, "orca")  #### pyscf

        # if DFT

        self.mf_lines = ['mf = dft.rks.RKS(mol, xc="{}")'.format(self.method)]
        self.import_lines.append("from pyscf import dft")

        if self.aux_basis_set != "":
            self.mf_lines.append(
                f'mf = mf.density_fit(auxbasis="{self.aux_basis_set}")'
            )

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def handle_pal_mem(self):
        self.pal = self.calc.nproc
        self.mem_per_core = int(self.calc.mem / self.calc.nproc)
        self.mem = self.calc.mem

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

            _element = ATOMIC_NUMBER[LOWERCASE_ATOMIC_SYMBOLS[element]]

            try:
                _rad = float(rad)
            except ValueError:
                raise InvalidParameter(
                    f"Invalid custom solvation radius for element {element}: '{rad}'"
                )
            self.solvation_radii[_element] = _rad

    def handle_solvation(self):
        return  ####

        if self.calc.parameters.solvent.lower() not in ["vacuum", ""]:
            solvent_keyword = get_solvent(
                self.calc.parameters.solvent,
                self.calc.parameters.software,
                solvation_model=self.calc.parameters.solvation_model,
            )
            model = self.calc.parameters.solvation_model
            radii_set = self.calc.parameters.solvation_radii
            custom_radii = self.calc.parameters.custom_solvation_radii

            if model == "smd":
                # Refined solvation radii
                # E. Engelage, N. Schulz, F. Heinen, S. M. Huber, D. G. Truhlar,
                # C. J. Cramer, Chem. Eur. J. 2018, 24, 15983–15987.
                if radii_set == "smd18":
                    if 53 not in self.solvation_radii:
                        self.solvation_radii[53] = 2.74
                    if 35 not in self.solvation_radii:
                        self.solvation_radii[35] = 2.60

                # self.add_to_block(
                #    "cpcm", ["smd true", f'SMDsolvent "{solvent_keyword}"', radii_specs]
                # )
            else:
                raise InvalidParameter(
                    f"Invalid solvation model for PySCF: "
                    f"'{self.calc.parameters.solvation_model}'"
                )

    def create_input_file(self):
        self.block_lines = ""

        mol_options = 'atom="""{}""",\n'.format(self.xyz_structure.strip())
        mol_options += 'basis="{}",\n'.format(self.calc.parameters.basis_set)
        mol_options += "charge={},\n".format(self.calc.charge)
        mol_options += "spin={},\n".format(self.calc.multiplicity - 1)

        for k, v in self.mol_specifications.items():
            if isinstance(v, int):
                mol_options += "{}={},\n".format(k, v)
            else:
                mol_options += '{}="{}",\n'.format(k, v)

        # from pyscf.geomopt.geometric_solver import optimize

        output = """
        import os, json, pyscf
        from pyscf.lib import num_threads
        {imports}

        os.environ["OMP_NUM_THREADS"] = "{pal}"
        os.environ["MKL_NUM_THREADS"] = "{pal}"
        os.environ["OPENBLAS_NUM_THREADS"] = "{pal}"
        num_threads({pal})

        mol = pyscf.gto.M(
        {mol_options})

        mol.max_memory = {mem}
        {mf_line}
        {specifications}{command_line}
        print(json.dumps({{
            {output_elements}
        }}))
        """

        extra_imports = "\n".join(self.import_lines)

        output_elements_str = ""

        for k, v in self.output_elements.items():
            output_elements_str += f'"{k}": {v},\n'

        raw = output.format(
            imports=extra_imports,
            pal=self.pal,
            mol_options=mol_options,
            mem=self.mem,
            mf_line="\n".join(self.mf_lines),
            specifications="\n".join(self.specifications_lines),
            command_line=self.command_line,
            output_elements=output_elements_str.strip(),
        )

        self.input_file = "\n".join([i.strip() for i in raw.split("\n")]).replace(
            "\n\n\n", "\n\n"
        )

    @property
    def output(self):
        return self.input_file

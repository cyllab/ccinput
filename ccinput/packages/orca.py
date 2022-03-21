import basis_set_exchange

from ccinput.constants import CalcType, ATOMIC_NUMBER, LOWERCASE_ATOMIC_SYMBOLS
from ccinput.utilities import get_method, get_basis_set, get_solvent, clean_xyz
from ccinput.exceptions import InvalidParameter, UnimplementedError


class OrcaCalculation:

    calc = None

    has_scan = False
    pal = 0
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

    def __init__(self, calc):
        self.calc = calc
        self.has_scan = False
        self.pal = 0
        self.blocks = []
        self.command_line = ""
        self.additional_commands = ""
        self.xyz_structure = ""
        self.block_lines = ""
        self.input_file = ""
        self.specifications = {}
        self.solvation_radii = {}

        self.handle_specifications()

        self.handle_command()
        self.handle_custom_basis_sets()
        self.handle_xyz()

        self.handle_pal()
        self.parse_custom_solvation_radii()
        self.handle_solvation()

        self.create_input_file()

    @property
    def confirmed_specifications(self):
        """Returns the effective additional commands; saved in the Parameters object"""
        return self.additional_commands

    def clean(self, s):
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. "
        )
        return "".join([c for c in s if c in WHITELIST])

    def handle_specifications(self):
        _specifications = (
            self.clean(self.calc.parameters.specifications).lower().strip()
        )

        specifications_list = []

        if _specifications != "":
            sspecs = _specifications.split()
            ind = 0
            while ind < len(sspecs):
                spec = sspecs[ind]
                if spec == "--phirshfeld":
                    HIRSHFELD_BLOCK = """%output
                    Print[ P_Hirshfeld] 1
                    end"""
                    self.blocks.append(HIRSHFELD_BLOCK)
                elif spec == "--nimages":
                    nimages = sspecs[ind + 1]
                    try:
                        nimages = int(nimages)
                    except ValueError:
                        raise InvalidParameter("Invalid specifications")
                    self.specifications["nimages"] = nimages
                    ind += 1
                elif spec not in specifications_list:
                    specifications_list.append(spec)

                ind += 1

        if self.calc.parameters.d3:
            specifications_list.append("d3zero")
        elif self.calc.parameters.d3bj:
            specifications_list.append("d3bj")

        if len(specifications_list) > 0:
            self.additional_commands = " ".join(specifications_list)

    def handle_command(self):
        if self.calc.type == CalcType.NMR:
            self.command_line = "NMR "
        elif self.calc.type == CalcType.OPT:
            self.command_line = "OPT "
        elif self.calc.type == CalcType.OPTFREQ:
            self.command_line = "OPT FREQ "
        elif self.calc.type == CalcType.TS:
            self.command_line = "OPTTS "
        elif self.calc.type == CalcType.MO:
            self.command_line = "SP "
            struct = clean_xyz(self.calc.xyz)

            electrons = 0
            for line in struct.split("\n"):
                if line.strip() == "":
                    continue
                el = line.split()[0]
                electrons += ATOMIC_NUMBER[el]

            electrons -= self.calc.charge

            if self.calc.multiplicity != 1:
                raise InvalidParameter("Unimplemented multiplicity")

            n_HOMO = int(electrons / 2) - 1
            n_LUMO = int(electrons / 2)
            n_LUMO1 = int(electrons / 2) + 1
            n_LUMO2 = int(electrons / 2) + 2

            mo_block = f"""%plots
                        dim1 45
                        dim2 45
                        dim3 45
                        min1 0
                        max1 0
                        min2 0
                        max2 0
                        min3 0
                        max3 0
                        Format Gaussian_Cube
                        MO("in-HOMO.cube",{n_HOMO},0);
                        MO("in-LUMO.cube",{n_LUMO},0);
                        MO("in-LUMOA.cube",{n_LUMO1},0);
                        MO("in-LUMOB.cube",{n_LUMO2},0);
                        end
                        """
            self.blocks.append(mo_block)
        elif self.calc.type == CalcType.FREQ:
            self.command_line = "FREQ "
        elif self.calc.type == CalcType.CONSTR_OPT:
            self.command_line = "OPT "

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
                scan_block = """%geom Scan
                {}
                end
                end"""
                self.blocks.append(scan_block.format("".join(scans).strip()))

            if len(freeze) > 0:
                freeze_block = """%geom Constraints
                {}
                end
                end"""
                self.blocks.append(freeze_block.format("".join(freeze).strip()))
        elif self.calc.type == CalcType.SP:
            self.command_line = "SP "
        elif self.calc.type == CalcType.MEP:  #### Second structure to handle
            self.command_line = "NEB "
            neb_block = """%neb
                        product "{}.xyz"
                        nimages {}
                        end"""
            if "nimages" in self.specifications:
                nimages = self.specifications["nimages"]
            else:
                nimages = 8
            self.blocks.append(neb_block.format(self.calc.aux_name, nimages))

        method = get_method(self.calc.parameters.method, "orca")
        if self.calc.parameters.theory_level not in [
            "xtb",
            "semiempirical",
            "special",
        ]:
            basis_set = get_basis_set(self.calc.parameters.basis_set, "orca")
            self.command_line += f"{method} {basis_set} "
        else:
            self.command_line += f"{method} "

    def handle_custom_basis_sets(self):
        if self.calc.parameters.custom_basis_sets == "":
            return

        entries = [
            i.strip()
            for i in self.calc.parameters.custom_basis_sets.split(";")
            if i.strip() != ""
        ]

        unique_atoms = []
        for line in self.calc.xyz.split("\n"):
            if line.strip() == "":
                continue
            a, *_ = line.split()
            if a not in unique_atoms:
                unique_atoms.append(a)

        BS_TEMPLATE = """%basis
        {}
        end"""

        custom_bs = ""
        for entry in entries:
            sentry = entry.split("=")

            if len(sentry) != 2:
                raise InvalidParameter("Invalid custom basis set string")

            el, bs_keyword = sentry

            if el not in unique_atoms:
                continue

            try:
                el_num = ATOMIC_NUMBER[el]
            except KeyError:
                raise InvalidParameter("Invalid atom in custom basis set string")

            bs = basis_set_exchange.get_basis(
                bs_keyword, fmt="ORCA", elements=[el_num], header=False
            ).strip()
            sbs = bs.split("\n")
            if bs.find("ECP") != -1:
                clean_bs = "\n".join(sbs[3:]).strip() + "\n"
                clean_bs = clean_bs.replace("\n$END", "$END").replace("$END", "end")
                custom_bs += f"newgto {el}\n"
                custom_bs += clean_bs
            else:
                clean_bs = "\n".join(sbs[3:-1]).strip() + "\n"
                custom_bs += f"newgto {el}\n"
                custom_bs += clean_bs
                custom_bs += "end"

        if custom_bs != "":
            self.blocks.append(BS_TEMPLATE.format(custom_bs))

    def handle_xyz(self):
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def handle_pal(self):
        if self.calc.parameters.theory_level == "semiempirical":
            self.pal = 1
        else:
            self.pal = self.calc.nproc

        pal_block = f"""%pal
        nprocs {self.pal}
        end"""

        self.blocks.append(pal_block)

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

    def get_radii_specs(self):
        radii_specs = ""
        for el, rad in self.solvation_radii.items():
            radii_specs += f"radius[{el}] {rad:.2f}\n"
        return radii_specs

    def handle_solvation(self):
        if self.calc.parameters.solvent.lower() not in ["vacuum", ""]:
            solvent_keyword = get_solvent(
                self.calc.parameters.solvent,
                self.calc.parameters.software,
                solvation_model=self.calc.parameters.solvation_model,
            )
            model = self.calc.parameters.solvation_model
            radii_set = self.calc.parameters.solvation_radii
            custom_radii = self.calc.parameters.custom_solvation_radii

            solv_block = ""

            if self.calc.parameters.method[:3].lower() == "xtb":
                self.command_line += f" ALPB({solvent_keyword})"
            elif model == "smd":
                # Refined solvation radii
                # E. Engelage, N. Schulz, F. Heinen, S. M. Huber, D. G. Truhlar,
                # C. J. Cramer, Chem. Eur. J. 2018, 24, 15983–15987.
                if radii_set == "smd18":
                    if 53 not in self.solvation_radii:
                        self.solvation_radii[53] = 2.74
                    if 35 not in self.solvation_radii:
                        self.solvation_radii[35] = 2.60

                radii_specs = self.get_radii_specs()
                solv_block = f"""%cpcm
                smd true
                SMDsolvent "{solvent_keyword}"
                {radii_specs}end"""
                self.blocks.append(solv_block)
            elif model == "cpcm":
                self.command_line += f"CPCM({solvent_keyword}) "
                radii_specs = self.get_radii_specs()
                if radii_specs != "":
                    solv_block = f"""%cpcm
                    {radii_specs}end"""
                    self.blocks.append(solv_block)
                if radii_set not in ["", "default", "bondi"]:
                    raise UnimplementedError(
                        "As of now, ccinput does not know how to request "
                        "other solvation radii than the default ones in ORCA with CPCM"
                    )
            else:
                raise InvalidParameter(
                    f"Invalid solvation model for ORCA: "
                    f"'{self.calc.parameters.solvation_model}'"
                )

    def create_input_file(self):
        self.block_lines = "\n".join(self.blocks)
        cmd = f"{self.command_line} {self.additional_commands}".replace("  ", " ")
        raw = self.TEMPLATE.format(
            cmd,
            self.calc.charge,
            self.calc.multiplicity,
            self.xyz_structure,
            self.block_lines,
        )
        self.input_file = "\n".join([i.strip() for i in raw.split("\n")])

    @property
    def output(self):
        return self.input_file

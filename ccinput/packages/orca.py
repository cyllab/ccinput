import os
import periodictable
import basis_set_exchange

from ccinput.constants import *
from ccinput.utilities import get_method, get_basis_set, get_solvent, clean_xyz

class OrcaCalculation:

    calc = None

    has_scan = False
    pal = 0
    blocks = []

    TEMPLATE = """!{}
    *xyz {} {}
    {}*
    {}"""
    #Command Line
    #Charge
    #Multiplicity
    #XYZ structure
    #Options blocks

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

        self.handle_specifications()

        self.handle_command()
        self.handle_custom_basis_sets()
        self.handle_xyz()

        self.handle_pal()
        self.handle_solvation()

        self.create_input_file()

    def clean(self, s):
        WHITELIST = set("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. ")
        return ''.join([c for c in s if c in WHITELIST])

    def handle_specifications(self):
        _specifications = self.clean(self.calc.parameters.specifications).lower().strip()
        if _specifications == '':
            return

        specifications_list = []
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
                nimages = sspecs[ind+1]
                try:
                    nimages = int(nimages)
                except ValueError:
                    raise Exception("Invalid specifications")
                self.specifications['nimages'] = nimages
                ind += 1
            elif spec not in specifications_list:
                specifications_list.append(spec)

            ind += 1

        if len(specifications_list) > 0:
            self.additional_commands = " ".join(specifications_list)

    def handle_command(self):
        if self.calc.type == CalcType.NMR:
            self.command_line = "NMR "
        elif self.calc.type == CalcType.OPT:
            self.command_line = "OPT "
        elif self.calc.type == CalcType.TS:
            self.command_line = "OPTTS "
        elif self.calc.type == CalcType.MO:
            self.command_line = "SP "
            struct = clean_xyz(self.calc.xyz)

            electrons = 0
            for line in struct.split('\n'):
                if line.strip() == "":
                    continue
                el = line.split()[0]
                electrons += ATOMIC_NUMBER[el]

            electrons -= self.calc.charge

            if self.calc.multiplicity != 1:
                raise Exception("Unimplemented multiplicity")

            n_HOMO = int(electrons/2)-1
            n_LUMO = int(electrons/2)
            n_LUMO1 = int(electrons/2)+1
            n_LUMO2 = int(electrons/2)+2

            mo_block = """%plots
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
                        MO("in-HOMO.cube",{},0);
                        MO("in-LUMO.cube",{},0);
                        MO("in-LUMOA.cube",{},0);
                        MO("in-LUMOB.cube",{},0);
                        end
                        """.format(n_HOMO, n_LUMO, n_LUMO1, n_LUMO2)
            self.blocks.append(mo_block)
        elif self.calc.type == CalcType.FREQ:
            self.command_line = "FREQ "
        elif self.calc.type == CalcType.CONSTR_OPT:
            self.command_line = "OPT "

            orca_constraints = ""
            has_scan = False
            scans = []
            freeze = []

            if self.calc.constraints.strip() == '':
                raise Exception("No constraints for constrained optimisation")

            for cmd in self.calc.constraints.split(';'):
                if cmd.strip() == '':
                    continue
                _cmd, ids = cmd.split('/')
                _cmd = _cmd.split('_')
                ids = ids.split('_')
                ids = [int(i)-1 for i in ids]
                type = len(ids)
                if _cmd[0] == "Scan":
                    has_scan = True
                else:
                    if type == 2:
                        freeze.append("{{ B {} {} C }}\n".format(*ids))
                    if type == 3:
                        freeze.append("{{ A {} {} {} C }}\n".format(*ids))
                    if type == 4:
                        freeze.append("{{ D {} {} {} {} C }}\n".format(*ids))
            if has_scan:
                for cmd in self.calc.constraints.split(';'):
                    if cmd.strip() == '':
                        continue
                    _cmd, ids = cmd.split('/')
                    ids = ids.split('_')
                    _cmd = _cmd.split('_')
                    ids_str = "{}".format(int(ids[0])-1)
                    for i in ids[1:]:
                        ids_str += " {}".format(int(i)-1)
                    if len(ids) == 2:
                        type = "B"
                    if len(ids) == 3:
                        type = "A"
                    if len(ids) == 4:
                        type = "D"
                    if _cmd[0] == "Scan":
                        scans.append("{} {} = {}, {}, {}\n".format(type, ids_str, *_cmd[1:]))

            self.has_scan = has_scan

            if len(scans) > 0:
                scan_block = """%geom Scan
                {}
                end
                end"""
                self.blocks.append(scan_block.format(''.join(scans).strip()))

            if len(freeze) > 0:
                freeze_block = """%geom Constraints
                {}
                end
                end"""
                self.blocks.append(freeze_block.format(''.join(freeze).strip()))
        elif self.calc.type == CalcType.SP:
            self.command_line = "SP "
        elif self.calc.type == CalcType.MEP:#### Second structure to handle
            self.command_line = "NEB "
            neb_block = """%neb
                        product "struct2.xyz"
                        nimages {}
                        end"""
            if 'nimages' in self.specifications:
                nimages = self.specifications['nimages']
            else:
                nimages = 8
            self.blocks.append(neb_block.format(nimages))

        method = get_method(self.calc.parameters.method, "orca")
        if self.calc.parameters.theory_level not in ['xtb', 'semi-empirical', 'special']:
            basis_set = get_basis_set(self.calc.parameters.basis_set, "orca")
            self.command_line += "{} {} ".format(method, basis_set)
        else:
            self.command_line += "{} ".format(method)


    def handle_custom_basis_sets(self):
        if self.calc.parameters.custom_basis_sets == "":
            return

        entries = [i.strip() for i in self.calc.parameters.custom_basis_sets.split(';') if i.strip() != ""]

        unique_atoms = []
        for line in self.calc.xyz.split('\n'):
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
            sentry = entry.split('=')

            if len(sentry) != 2:
                raise Exception("Invalid custom basis set string")

            el, bs_keyword = sentry

            if el not in unique_atoms:
                continue

            try:
                el_num = ATOMIC_NUMBER[el]
            except KeyError:
                raise Exception("Invalid atom in custom basis set string")

            bs = basis_set_exchange.get_basis(bs_keyword, fmt='ORCA', elements=[el_num], header=False).strip()
            sbs = bs.split('\n')
            if bs.find("ECP") != -1:
                clean_bs = '\n'.join(sbs[3:]).strip() + '\n'
                clean_bs = clean_bs.replace("\n$END", '$END').replace('$END', 'end')
                custom_bs += "newgto {}\n".format(el)
                custom_bs += clean_bs
            else:
                clean_bs = '\n'.join(sbs[3:-1]).strip() + '\n'
                custom_bs += "newgto {}\n".format(el)
                custom_bs += clean_bs
                custom_bs += "end"

        if custom_bs != "":
            self.blocks.append(BS_TEMPLATE.format(custom_bs))

    def handle_xyz(self):
        lines = [i + '\n' for i in clean_xyz(self.calc.xyz).split('\n') if i != '' ]
        self.xyz_structure = ''.join(lines)

    def handle_pal(self):
        if self.calc.parameters.theory_level == "semi-empirical":
            self.pal = 1
        else:
            self.pal = self.calc.nproc

        pal_block = """%pal
        nprocs {}
        end""".format(self.pal)

        self.blocks.append(pal_block)

    def handle_solvation(self):
        if self.calc.parameters.solvent.lower() not in ["vacuum", ""]:
            solvent_keyword = get_solvent(self.calc.parameters.solvent, self.calc.parameters.software, solvation_model=self.calc.parameters.solvation_model)

            if self.calc.parameters.method[:3] == 'gfn':
                self.command_line += " ALPB({})".format(solvent_keyword)
            elif self.calc.parameters.solvation_model == "smd":
                if self.calc.parameters.solvation_radii in ["default", ""]:
                    smd_block = '''%cpcm
                    smd true
                    SMDsolvent "{}"
                    end'''.format(solvent_keyword)
                    self.blocks.append(smd_block)
                elif self.calc.parameters.solvation_radii == "smd18":
                    smd_block = '''%cpcm
                    smd true
                    SMDsolvent "{}"
                    radius[53] 2.74
                    radius[35] 2.60
                    end'''.format(solvent_keyword)
                    self.blocks.append(smd_block)
            elif self.calc.parameters.solvation_model == "cpcm":
                self.command_line += "CPCM({}) ".format(solvent_keyword)
                ###CPCM radii
            else:
                raise Exception("Invalid solvation model for ORCA: '{}'".format(self.calc.parameters.solvation_model))

    def create_input_file(self):
        self.block_lines = '\n'.join(self.blocks)
        cmd = "{} {}".format(self.command_line, self.additional_commands).replace('  ', ' ')
        raw = self.TEMPLATE.format(cmd, self.calc.charge, self.calc.multiplicity, self.xyz_structure, self.block_lines)
        self.input_file = '\n'.join([i.strip() for i in raw.split('\n')])


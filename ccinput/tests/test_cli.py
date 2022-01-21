import shlex
import os
from unittest import TestCase

from ccinput.calculation import Calculation, Parameters
from ccinput.wrapper import gen_input, get_input_from_args, get_parser

class ManualCliTests(TestCase):
    def are_equivalent(self, api_args, cmd_line):
        ref = gen_input(**api_args)

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))
        inp = get_input_from_args(args)
        return ref == inp

    def struct(self, name):
        return os.path.join('/'.join(__file__.split('/')[:-1]), "structures/", name + '.xyz')

    def test_basic(self):
        args = {
            'software': "gaussian",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_no_nproc(self):
        args = {
            'software': "gaussian",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' --mem 1G -c -1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_basic_file(self):
        args = {
            'software': "gaussian",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
        }
        line = f"gaussian sp HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G"
        self.assertTrue(self.are_equivalent(args, line))

    def test_freeze_distance(self):
        args = {
            'software': "gaussian",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Freeze/1_2;",
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 1 2"
        self.assertTrue(self.are_equivalent(args, line))

    def test_freeze_angle(self):
        args = {
            'software': "gaussian",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Freeze/3_2_4;",
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 3 2 4"
        self.assertTrue(self.are_equivalent(args, line))

    def test_freeze_dihedral(self):
        args = {
            'software': "gaussian",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Freeze/6_3_1_2;",
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 6 3 1 2"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance_multiple(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --nsteps 10 --scan 3 4 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance_multiple_step(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --step -0.1 --scan 3 4 --from 2.0 --to 1.0 --step -0.1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance_step_wrong_sign(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --step 0.1 --scan 3 4 --from 2.0 --to 1.0 --step 0.1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_angle(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2_3;",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 3 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_dihedral(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'in_file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2_3_4;",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 3 4 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_name(self):
        args = {
            'software': "gaussian",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
            'name': "Chloride in vacuum",
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --name 'Chloride in vacuum'"
        self.assertTrue(self.are_equivalent(args, line))

    def test_custom_basis_set(self):
        args = {
            'software': "gaussian",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'custom_basis_sets': "Cl=Def2-SVPD;",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 -cbs 'Cl=Def2-SVPD;'"
        self.assertTrue(self.are_equivalent(args, line))

    def test_solvation(self):
        args = {
            'software': "gaussian",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
            'solvation_model': 'SMD',
            'solvent': 'Chloroform',
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --solvent chloroform --solvation_model smd"
        self.assertTrue(self.are_equivalent(args, line))

    def test_solvation_smd18(self):
        args = {
            'software': "gaussian",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
            'solvation_model': 'SMD',
            'solvent': 'Chloroform',
            'solvation_radii': 'SMD18',
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --solvent chloroform --solvation_model smd --solvation_radii SMD18"
        self.assertTrue(self.are_equivalent(args, line))


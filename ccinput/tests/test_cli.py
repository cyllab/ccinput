import shlex
import os

from ccinput.calculation import Calculation, Parameters
from ccinput.wrapper import gen_input, get_input_from_args, get_parser
from ccinput.tests.testing_utilities import InputTests

class CliEquivalenceTests(InputTests):
    def are_equivalent(self, api_args, cmd_line):
        ref = gen_input(**api_args)

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))
        calcs, outputs = get_input_from_args(args)
        inp = calcs[0].input_file
        return self.is_equivalent(ref, inp)

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
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'name': 'ethanol',
        }
        line = f"gaussian sp HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G"
        self.assertTrue(self.are_equivalent(args, line))

    def test_freeze_distance(self):
        args = {
            'software': "gaussian",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Freeze/1_2;",
            'name': 'ethanol',
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 1 2"
        self.assertTrue(self.are_equivalent(args, line))

    def test_freeze_angle(self):
        args = {
            'software': "gaussian",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Freeze/3_2_4;",
            'name': 'ethanol',
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 3 2 4"
        self.assertTrue(self.are_equivalent(args, line))

    def test_freeze_dihedral(self):
        args = {
            'software': "gaussian",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Freeze/6_3_1_2;",
            'name': 'ethanol',
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 6 3 1 2"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;",
            'name': 'ethanol',
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance_multiple(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            'name': 'ethanol',
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --nsteps 10 --scan 3 4 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance_multiple_no_from(self):
        args = {
            'software': "gaussian",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            'name': 'ethanol',
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --to 1.0 --nsteps 10 --scan 3 4 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance_multiple_step(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            'name': 'ethanol',
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --step -0.1 --scan 3 4 --from 2.0 --to 1.0 --step -0.1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_distance_step_wrong_sign(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            'name': 'ethanol',
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --step 0.1 --scan 3 4 --from 2.0 --to 1.0 --step 0.1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_angle(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2_3;",
            'name': 'ethanol',
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 3 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.are_equivalent(args, line))

    def test_scan_dihedral(self):
        args = {
            'software': "orca",
            'type': "constr_opt",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct('ethanol'),
            'nproc': 1,
            'mem': "1G",
            'constraints': "Scan_2.0_1.0_10/1_2_3_4;",
            'name': 'ethanol',
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

    def test_opt_freq(self):
        args = {
            'software': "gaussian",
            'type': "opt+freq",
            'method': "HF",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
        }
        line = "gaussian opt+freq HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_opt_freq2(self):
        args = {
            'software': "gaussian",
            'type': "opt+freq",
            'method': "HF",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
        }
        line = "gaussian opt-freq HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
        self.assertTrue(self.are_equivalent(args, line))

    def test_d3_gaussian(self):
        args = {
            'software': "gaussian",
            'type': "opt",
            'method': "M06",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
            'd3': True,
        }
        line = "gaussian opt M06 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3"
        self.assertTrue(self.are_equivalent(args, line))

    def test_d3bj_gaussian(self):
        args = {
            'software': "gaussian",
            'type': "opt",
            'method': "PBE0",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
            'd3bj': True,
        }
        line = "gaussian opt PBE0 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3bj"
        self.assertTrue(self.are_equivalent(args, line))

    def test_d3_orca(self):
        args = {
            'software': "ORCA",
            'type': "opt",
            'method': "M06",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
            'd3': True,
        }
        line = "orca opt M06 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3"
        self.assertTrue(self.are_equivalent(args, line))

    def test_d3bj_orca(self):
        args = {
            'software': "ORCA",
            'type': "opt",
            'method': "PBE0",
            'basis_set': "Def2SVP",
            'xyz': "Cl 0 0 0\n",
            'nproc': 1,
            'mem': "1G",
            'charge': -1,
            'd3bj': True,
        }
        line = "orca opt PBE0 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3bj"
        self.assertTrue(self.are_equivalent(args, line))

class ManualCliTests(InputTests):
    def test_multiple_files_no_output(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        assert len(objs) == 2
        assert len(outputs) == 0

        args1 = {
            'software': "ORCA",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct("ethanol"),
            'nproc': 1,
            'mem': "1G",
        }
        args2 = {
            'software': "ORCA",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct("CH4"),
            'nproc': 1,
            'mem': "1G",
        }

        inp1 = gen_input(**args1)
        inp2 = gen_input(**args2)

        assert self.is_equivalent(inp1, objs[0].input_file)
        assert self.is_equivalent(inp2, objs[1].input_file)

    def test_multiple_files_output(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o test.inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        assert len(objs) == 2
        assert len(outputs) == 2

        assert os.path.basename(outputs[0]) == "test_ethanol.inp"
        assert os.path.basename(outputs[1]) == "test_CH4.inp"

        args1 = {
            'software': "ORCA",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct("ethanol"),
            'nproc': 1,
            'mem': "1G",
        }
        args2 = {
            'software': "ORCA",
            'type': "sp",
            'method': "HF",
            'basis_set': "Def2SVP",
            'file': self.struct("CH4"),
            'nproc': 1,
            'mem': "1G",
        }

        inp1 = gen_input(**args1)
        inp2 = gen_input(**args2)

        assert self.is_equivalent(inp1, objs[0].input_file)
        assert self.is_equivalent(inp2, objs[1].input_file)

    def test_multiple_files_output_no_name(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o .inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        assert len(objs) == 2
        assert len(outputs) == 2

        assert outputs[0] == "ethanol.inp"
        assert outputs[1] == "CH4.inp"

    def test_multiple_files_output_name_no_override(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o .inp -n 1 --mem 1G --name test"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        assert len(objs) == 2
        assert len(outputs) == 2

        assert outputs[0] == "ethanol.inp"
        assert outputs[1] == "CH4.inp"

    def test_multiple_files_output_directory(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o calc_dir/.inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        assert len(objs) == 2
        assert len(outputs) == 2

        assert outputs[0] == "calc_dir/ethanol.inp"
        assert outputs[1] == "calc_dir/CH4.inp"

    def test_single_file_output(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} -o calc_dir/ethanol.inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        assert len(objs) == 1
        assert len(outputs) == 1

        assert outputs[0] == "calc_dir/ethanol.inp"


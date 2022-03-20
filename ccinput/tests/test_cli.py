import shlex
import os
from mock import patch
import tempfile
import json

from ccinput.calculation import Calculation, Parameters
from ccinput.wrapper import gen_input, get_input_from_args, get_parser, cmd
from ccinput.tests.testing_utilities import InputTests
from ccinput import presets

from contextlib import contextmanager, redirect_stdout
from os import devnull

PRESET_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "presets")


@contextmanager
def hide_cmd_output():
    with open(devnull, "w") as gobble:
        with redirect_stdout(gobble) as out:
            yield out


class CliEquivalenceTests(InputTests):
    def test_basic(self):
        args = {
            "software": "gaussian",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_no_nproc(self):
        args = {
            "software": "gaussian",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' --mem 1G -c -1"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_basic_file(self):
        args = {
            "software": "gaussian",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "name": "ethanol",
        }
        line = f"gaussian sp HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_freeze_distance(self):
        args = {
            "software": "gaussian",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Freeze/1_2;",
            "name": "ethanol",
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 1 2"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_freeze_angle(self):
        args = {
            "software": "gaussian",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Freeze/3_2_4;",
            "name": "ethanol",
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 3 2 4"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_freeze_dihedral(self):
        args = {
            "software": "gaussian",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Freeze/6_3_1_2;",
            "name": "ethanol",
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --freeze 6 3 1 2"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_scan_distance(self):
        args = {
            "software": "orca",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Scan_2.0_1.0_10/1_2;",
            "name": "ethanol",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_scan_distance_multiple(self):
        args = {
            "software": "orca",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            "name": "ethanol",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --nsteps 10 --scan 3 4 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_scan_distance_multiple_no_from(self):
        args = {
            "software": "gaussian",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            "name": "ethanol",
        }
        line = f"gaussian constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --to 1.0 --nsteps 10 --scan 3 4 --to 1.0 --nsteps 10"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_scan_distance_multiple_step(self):
        args = {
            "software": "orca",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            "name": "ethanol",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --step -0.1 --scan 3 4 --from 2.0 --to 1.0 --step -0.1"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_scan_distance_step_wrong_sign(self):
        args = {
            "software": "orca",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Scan_2.0_1.0_10/1_2;Scan_2.0_1.0_10/3_4;",
            "name": "ethanol",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 --from 2.0 --to 1.0 --step 0.1 --scan 3 4 --from 2.0 --to 1.0 --step 0.1"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_scan_angle(self):
        args = {
            "software": "orca",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Scan_2.0_1.0_10/1_2_3;",
            "name": "ethanol",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 3 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_scan_dihedral(self):
        args = {
            "software": "orca",
            "type": "constr_opt",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
            "constraints": "Scan_2.0_1.0_10/1_2_3_4;",
            "name": "ethanol",
        }
        line = f"orca constr_opt HF -bs Def2SVP -f {self.struct('ethanol')} -n 1 --mem 1G --scan 1 2 3 4 --from 2.0 --to 1.0 --nsteps 10"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_name(self):
        args = {
            "software": "gaussian",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "name": "Chloride in vacuum",
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --name 'Chloride in vacuum'"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_custom_basis_set(self):
        args = {
            "software": "gaussian",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "custom_basis_sets": "Cl=Def2-SVPD;",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 -cbs 'Cl=Def2-SVPD;'"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_solvation(self):
        args = {
            "software": "gaussian",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "solvation_model": "SMD",
            "solvent": "Chloroform",
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --solvent chloroform --solvation_model smd"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_solvation_smd18(self):
        args = {
            "software": "gaussian",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "solvation_model": "SMD",
            "solvent": "Chloroform",
            "solvation_radii": "SMD18",
        }
        line = "gaussian sp HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --solvent chloroform --solvation_model smd --solvation_radii SMD18"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_opt_freq(self):
        args = {
            "software": "gaussian",
            "type": "opt+freq",
            "method": "HF",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
        }
        line = "gaussian opt+freq HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_opt_freq2(self):
        args = {
            "software": "gaussian",
            "type": "opt+freq",
            "method": "HF",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
        }
        line = "gaussian opt-freq HF -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_d3_gaussian(self):
        args = {
            "software": "gaussian",
            "type": "opt",
            "method": "M06",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "d3": True,
        }
        line = "gaussian opt M06 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_d3bj_gaussian(self):
        args = {
            "software": "gaussian",
            "type": "opt",
            "method": "PBE0",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "d3bj": True,
        }
        line = (
            "gaussian opt PBE0 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3bj"
        )
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_d3_orca(self):
        args = {
            "software": "ORCA",
            "type": "opt",
            "method": "M06",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "d3": True,
        }
        line = "orca opt M06 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_d3bj_orca(self):
        args = {
            "software": "ORCA",
            "type": "opt",
            "method": "PBE0",
            "basis_set": "Def2SVP",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "d3bj": True,
        }
        line = "orca opt PBE0 -bs Def2SVP --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 --d3bj"
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_preset_winchars(self):
        presets.data_dir = PRESET_DIR
        args = {
            "software": "gaussian",
            "type": "opt",
            "method": "HF",
            "basis_set": "def2tzvp",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
        }

        line = '--preset winchars --xyz "Cl 0 0 0" -c -1'
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_preset_whitespace(self):
        presets.data_dir = PRESET_DIR
        args = {
            "software": "gaussian",
            "type": "opt",
            "method": "HF",
            "basis_set": "def2tzvp",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
        }

        line = '--preset whitespace --xyz "Cl 0 0 0" -c -1'
        self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_preset_additive_specifications(self):
        presets.data_dir = PRESET_DIR
        args = {
            "software": "gaussian",
            "type": "opt",
            "method": "HF",
            "basis_set": "def2tzvp",
            "xyz": "Cl 0 0 0\n",
            "nproc": 1,
            "mem": "1G",
            "charge": -1,
            "specifications": "nosymm 5d",
        }

        line = '--preset specification --xyz "Cl 0 0 0" -c -1 --specification "5d"'
        self.assertTrue(self.args_cmd_equivalent(args, line))


class ManualCliTests(InputTests):
    def setUp(self):
        self.warnings = []

    def get_warn(self, msg):
        self.warnings.append(msg)

    def test_multiple_files_no_output(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 2)
        self.assertEqual(len(outputs), 0)

        args1 = {
            "software": "ORCA",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
        }
        args2 = {
            "software": "ORCA",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("CH4"),
            "nproc": 1,
            "mem": "1G",
        }

        inp1 = gen_input(**args1)
        inp2 = gen_input(**args2)

        self.assertTrue(self.is_equivalent(inp1, objs[0].input_file))
        self.assertTrue(self.is_equivalent(inp2, objs[1].input_file))

    def test_multiple_files_output(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o test.inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)

        self.assertEqual(len(objs), 2)
        self.assertEqual(len(outputs), 2)

        self.assertEqual(os.path.basename(outputs[0]), "test_ethanol.inp")
        self.assertEqual(os.path.basename(outputs[1]), "test_CH4.inp")

        args1 = {
            "software": "ORCA",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("ethanol"),
            "nproc": 1,
            "mem": "1G",
        }
        args2 = {
            "software": "ORCA",
            "type": "sp",
            "method": "HF",
            "basis_set": "Def2SVP",
            "file": self.struct("CH4"),
            "nproc": 1,
            "mem": "1G",
        }

        inp1 = gen_input(**args1)
        inp2 = gen_input(**args2)

        self.assertTrue(self.is_equivalent(inp1, objs[0].input_file))
        self.assertTrue(self.is_equivalent(inp2, objs[1].input_file))

    def test_multiple_files_output_no_name(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o .inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 2)
        self.assertEqual(len(outputs), 2)

        self.assertEqual(outputs[0], "ethanol.inp")
        self.assertEqual(outputs[1], "CH4.inp")

    def test_multiple_files_output_name_no_override(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o .inp -n 1 --mem 1G --name test"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 2)
        self.assertEqual(len(outputs), 2)

        self.assertEqual(outputs[0], "ethanol.inp")
        self.assertEqual(outputs[1], "CH4.inp")

    def test_multiple_files_output_directory(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} {self.struct('CH4')} -o calc_dir/.inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 2)
        self.assertEqual(len(outputs), 2)

        self.assertEqual(outputs[0], "calc_dir/ethanol.inp")
        self.assertEqual(outputs[1], "calc_dir/CH4.inp")

    def test_single_file_output(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} -o calc_dir/ethanol.inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 1)
        self.assertEqual(len(outputs), 1)

        self.assertEqual(outputs[0], "calc_dir/ethanol.inp")

    def test_single_file_output_change_filename(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} -o calc_dir/sp_HF.inp -n 1 --mem 1G"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 1)
        self.assertEqual(len(outputs), 1)

        self.assertEqual(outputs[0], "calc_dir/sp_HF.inp")

    def test_single_file_output_change_name(self):
        cmd_line = f"orca sp HF -bs Def2SVP -f {self.struct('ethanol')} -o calc_dir/sp_HF.inp -n 1 --mem 1G --name my_sp"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 1)
        self.assertEqual(len(outputs), 1)

        self.assertEqual(outputs[0], "calc_dir/sp_HF.inp")

    def test_xtb_commandline(self):
        cmd_line = f"xtb sp gfn2-xtb -f {self.struct('ethanol')}"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 1)
        self.assertEqual(len(outputs), 0)

        self.assertEqual(objs[0].input_file, "")
        self.assertEqual(objs[0].command, "xtb ethanol.xyz")

    def test_crest_commandline(self):
        cmd_line = f"xtb 'constr conf search' gfn2-xtb -f {self.struct('ethanol')} --freeze 1 3"

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))

        objs, outputs = get_input_from_args(args)
        self.assertEqual(len(objs), 1)
        self.assertEqual(len(outputs), 0)

        INPUT_REF = """
        $constrain
        force constant=1.0
        reference=ethanol.xyz
        distance: 1, 3, auto
        atoms: 1,3
        $metadyn
        atoms: 2,4-9"""

        self.assertTrue(self.is_equivalent(INPUT_REF, objs[0].input_file))
        self.assertEqual(
            objs[0].command, "crest ethanol.xyz -cinp input -rthr 0.6 -ewin 6"
        )

    @patch("ccinput.utilities.warn")
    def test_warn_unknown(self, warn_fn):
        warn_fn.side_effect = self.get_warn
        line = 'Gaussian sp abcd -bs cc-pvdz --xyz "Cl 0 0 0" -c -1'
        parser = get_parser()
        args = parser.parse_args(shlex.split(line))

        objs, outputs = get_input_from_args(args)

        self.assertNotEqual(len(self.warnings), 0)

    @patch("ccinput.utilities.warn")
    def test_no_warn_utpsstpss(self, warn_fn):
        warn_fn.side_effect = self.get_warn
        line = 'Gaussian sp utpsstpss -bs cc-pvdz --xyz "Cl 0 0 0" -c -1'
        parser = get_parser()
        args = parser.parse_args(shlex.split(line))

        objs, outputs = get_input_from_args(args)

        self.assertEqual(len(self.warnings), 0)

    @patch("ccinput.utilities.warn")
    def test_no_warn_rtpsstpss(self, warn_fn):
        warn_fn.side_effect = self.get_warn
        line = 'Gaussian sp rtpsstpss -bs cc-pvdz --xyz "Cl 0 0 0" -c -1'
        parser = get_parser()
        args = parser.parse_args(shlex.split(line))

        objs, outputs = get_input_from_args(args)

        self.assertEqual(len(self.warnings), 0)

    @patch("ccinput.utilities.warn")
    def test_no_warn_tpsstpss(self, warn_fn):
        warn_fn.side_effect = self.get_warn
        line = 'Gaussian sp tpsstpss -bs cc-pvdz --xyz "Cl 0 0 0" -c -1'
        parser = get_parser()
        args = parser.parse_args(shlex.split(line))

        objs, outputs = get_input_from_args(args)

        self.assertEqual(len(self.warnings), 0)

    def test_synonym_basis_set(self):
        line = 'Gaussian sp tpsstpss -bs ccpvdz --xyz "Cl 0 0 0" -c -1'
        parser = get_parser()
        args = parser.parse_args(shlex.split(line))

        objs, outputs = get_input_from_args(args)


class CliPresetTests(InputTests):
    def test_create_preset(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "Gaussian sp tpsstpss -bs ccpvdz --save my_preset"
            cmd(cmd_line=line)

            content = os.listdir(tmp_dir)
            self.assertEqual(len(content), 1)
            self.assertEqual(content[0], "my_preset.preset")

    def test_create_partial_preset(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "Gaussian sp -bs ccpvdz --save my_preset"
            cmd(cmd_line=line)

            content = os.listdir(tmp_dir)
            self.assertEqual(len(content), 1)
            self.assertEqual(content[0], "my_preset.preset")

    def test_integrity_preset1(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "Gaussian sp tpsstpss -bs ccpvdz --save my_preset"
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertEqual(preset["software"], "Gaussian")
            self.assertEqual(preset["type"], "sp")
            self.assertEqual(preset["method"], "tpsstpss")
            self.assertEqual(preset["basis_set"], "ccpvdz")
            self.assertEqual(len(preset.keys()), 5)  # 4 parameters + version

    def test_integrity_preset2(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = 'orca opt m062x -bs def2tzvp --specifications "tightscf" -sm smd -sr smd18 -s methanol --save my_preset'
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertEqual(preset["software"], "orca")
            self.assertEqual(preset["type"], "opt")
            self.assertEqual(preset["method"], "m062x")
            self.assertEqual(preset["basis_set"], "def2tzvp")
            self.assertEqual(preset["specifications"], "tightscf")
            self.assertEqual(preset["solvation_model"], "smd")
            self.assertEqual(preset["solvation_radii"], "smd18")
            self.assertEqual(preset["solvent"], "methanol")

            self.assertEqual(len(preset.keys()), 9)

    def test_integrity_partial_preset(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = 'orca opt -bs def2tzvp --specifications "tightscf" -sm smd -sr smd18 -s methanol --save my_preset'
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertEqual(preset["software"], "orca")
            self.assertEqual(preset["type"], "opt")
            self.assertEqual(preset["basis_set"], "def2tzvp")
            self.assertEqual(preset["specifications"], "tightscf")
            self.assertEqual(preset["solvation_model"], "smd")
            self.assertEqual(preset["solvation_radii"], "smd18")
            self.assertEqual(preset["solvent"], "methanol")

            self.assertNotIn("method", preset)

            self.assertEqual(len(preset.keys()), 8)

    def test_structure_not_saved(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = 'orca opt m062x -bs def2tzvp --specifications "tightscf" -sm smd -sr smd18 -s methanol --xyz "Cl 0 0 0" -c -1 --save my_preset'
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertNotIn("xyz", preset)
            self.assertEqual(len(preset.keys()), 10)

    def test_override_partially(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = 'orca opt m062x -bs def2tzvp --specifications "tightscf" -sm smd -sr smd18 -s methanol --save my_preset'
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertEqual(preset["basis_set"], "def2tzvp")
            self.assertEqual(len(preset.keys()), 9)

            line = "orca -bs def2qzvp --save my_preset"
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertEqual(preset["software"], "orca")
            self.assertEqual(preset["type"], "opt")
            self.assertEqual(preset["method"], "m062x")
            self.assertEqual(preset["specifications"], "tightscf")
            self.assertEqual(preset["solvation_model"], "smd")
            self.assertEqual(preset["solvation_radii"], "smd18")
            self.assertEqual(preset["solvent"], "methanol")

            self.assertEqual(preset["basis_set"], "def2qzvp")
            self.assertEqual(len(preset.keys()), 9)

    def test_override_completely(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = 'orca opt m062x -bs def2tzvp --specifications "tightscf" -sm smd -sr smd18 -s methanol --save my_preset'
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertEqual(preset["basis_set"], "def2tzvp")
            self.assertEqual(len(preset.keys()), 9)

            line = 'g16 sp m06 -bs def2qzvp --specifications "opt(maxstep=10) symm" -sm pcm -sr bondi -s acetonitrile --save my_preset'
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertEqual(preset["software"], "g16")
            self.assertEqual(preset["type"], "sp")
            self.assertEqual(preset["method"], "m06")
            self.assertEqual(preset["basis_set"], "def2qzvp")
            self.assertEqual(preset["specifications"], "opt(maxstep=10) symm")
            self.assertEqual(preset["solvation_model"], "pcm")
            self.assertEqual(preset["solvation_radii"], "bondi")
            self.assertEqual(preset["solvent"], "acetonitrile")

            self.assertEqual(len(preset.keys()), 9)

    def test_load_preset(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt m062x -bs def2tzvp --save my_preset"
            cmd(cmd_line=line)

            args = {
                "software": "orca",
                "type": "opt",
                "method": "m062x",
                "basis_set": "def2tzvp",
                "xyz": "Cl 0 0 0\n",
                "nproc": 1,
                "mem": "1G",
                "charge": -1,
            }
            line = "--preset my_preset --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
            self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_preset_override1(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt m062x -bs def2tzvp --save my_preset"
            cmd(cmd_line=line)

            args = {
                "software": "gaussian",
                "type": "opt",
                "method": "m062x",
                "basis_set": "def2tzvp",
                "xyz": "Cl 0 0 0\n",
                "nproc": 1,
                "mem": "1G",
                "charge": -1,
            }
            line = "gaussian --preset my_preset --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
            self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_preset_override2(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt m062x -bs def2tzvp --save my_preset"
            cmd(cmd_line=line)

            args = {
                "software": "orca",
                "type": "sp",
                "method": "m06",
                "basis_set": "def2tzvp",
                "xyz": "Cl 0 0 0\n",
                "nproc": 1,
                "mem": "1G",
                "charge": -1,
            }
            line = "orca sp m06 --preset my_preset --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1"
            self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_preset_override3(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt m062x -bs def2tzvp --save my_preset"
            cmd(cmd_line=line)

            args = {
                "software": "orca",
                "type": "opt",
                "method": "m062x",
                "basis_set": "def2svp",
                "xyz": "Cl 0 0 0\n",
                "nproc": 1,
                "mem": "1G",
                "charge": -1,
            }
            line = "--preset my_preset --xyz 'Cl 0 0 0' -n 1 --mem 1G -c -1 -bs def2svp"
            self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_preset_override4(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt m062x -bs def2tzvp --save my_preset"
            cmd(cmd_line=line)

            args = {
                "software": "orca",
                "type": "opt",
                "method": "m062x",
                "basis_set": "def2tzvp",
                "xyz": "Cl 0 0 0\n",
                "nproc": 1,
                "mem": "1G",
                "charge": -1,
                "d3": True,
            }
            line = "--xyz 'Cl 0 0 0' -n 1 --preset my_preset --mem 1G -c -1 --d3"
            self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_preset_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt -bs def2tzvp --save my_preset"
            cmd(cmd_line=line)

            line = "--xyz 'Cl 0 0 0' -c -1 -n 1 --mem 1G --preset my_preset"

            with self.assertRaises(SystemExit):
                cmd(cmd_line=line)

    def test_load_preset_complete1(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt -bs def2tzvp --save my_preset"
            cmd(cmd_line=line)

            line = (
                "orca opt m062x --xyz 'Cl 0 0 0' -c -1 -n 1 --mem 1G --preset my_preset"
            )

            cmd(cmd_line=line)

    def test_load_preset_complete2(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt m062x --save my_preset"
            cmd(cmd_line=line)

            args = {
                "software": "orca",
                "type": "opt",
                "method": "m062x",
                "basis_set": "def2tzvp",
                "xyz": "Cl 0 0 0\n",
                "nproc": 1,
                "mem": "1G",
                "charge": -1,
                "d3": True,
            }
            line = "--xyz 'Cl 0 0 0' -n 1 -bs def2tzvp --preset my_preset --mem 1G -c -1 --d3"
            self.assertTrue(self.args_cmd_equivalent(args, line))

    def test_load_unknown_preset(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "orca opt m062x --save my_preset"
            cmd(cmd_line=line)

            line = "--xyz 'Cl 0 0 0' -n 1 -bs def2tzvp --preset other_preset --mem 1G -c -1 --d3"

            with self.assertRaises(SystemExit):
                cmd(cmd_line=line)

    def test_structures_not_saved(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "Gaussian sp tpsstpss -bs ccpvdz -f mystruct.xyz --save my_preset"
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertNotIn("file", preset)

    def test_output_not_saved(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "Gaussian sp tpsstpss -bs ccpvdz -o calc.com --save my_preset"
            cmd(cmd_line=line)

            with open(os.path.join(tmp_dir, "my_preset.preset")) as f:
                preset = json.load(f)

            self.assertNotIn("output", preset)

    def test_save_and_load(self):
        with tempfile.TemporaryDirectory() as tmp_dir, hide_cmd_output():
            presets.data_dir = tmp_dir
            line = "Gaussian sp tpsstpss -bs ccpvdz -o calc.com --save my_preset --preset my_preset"
            with self.assertRaises(SystemExit):
                cmd(cmd_line=line)

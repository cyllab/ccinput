import os
import shlex
import subprocess
import tempfile
from unittest import TestCase

from ccinput.wrapper import gen_obj, gen_input, get_input_from_args, get_parser


class InputTests(TestCase):
    def generate_calculation(self, **params):
        params["file"] = os.path.join(
            "/".join(__file__.split("/")[:-1]), "structures/", params["file"]
        )
        return gen_obj(**params)

    def args_cmd_equivalent(self, api_args, cmd_line):
        ref = gen_input(**api_args)

        parser = get_parser()
        args = parser.parse_args(shlex.split(cmd_line))
        default_params = vars(parser.parse_args([]))
        calcs, outputs = get_input_from_args(args, default_params=default_params)
        inp = calcs[0].input_file
        return self.is_equivalent(ref, inp)

    def is_equivalent(self, ref, res):
        ref_lines = [i.strip() for i in ref.strip().split("\n")]
        res_lines = [i.strip() for i in res.strip().split("\n")]

        if len(ref_lines) != len(res_lines):
            print(f"Different number of lines: {len(ref_lines)} and {len(res_lines)}")
            print("----START REFERENCE---")
            print("\n".join(ref_lines))
            print("----END REFERENCE---")
            print("----START RESULT---")
            print("\n".join(res_lines))
            print("----END RESULT---")
            return False

        for line1, line2 in zip(ref_lines, res_lines):
            if line1 != line2:
                if self.xyz_line_equivalent(line1, line2):
                    continue
                print("")
                print("Difference found:")
                print("----REFERENCE---")
                print(line1)
                print()
                print("----RESULT---")
                print(line2)
                print("")
                return False

        return True

    def xyz_line_equivalent(self, line1, line2):
        sline1 = line1.strip().split()
        sline2 = line2.strip().split()

        if len(sline1) != 4 or len(sline2) != 4:
            return False

        if sline1[0] != sline2[0]:
            return False

        for c1, c2 in zip(sline1[1:], sline2[1:]):
            try:
                _c1 = float(c1)
                _c2 = float(c2)
            except ValueError:
                return False
            if abs(_c1 - _c2) > 1e-7:
                return False
        return True

    def struct(self, name):
        return os.path.join(
            "/".join(__file__.split("/")[:-1]), "structures/", name + ".xyz"
        )


class CalculationTests(InputTests):
    def run_calc(self, **params):
        obj = self.generate_calculation(**params)
        with tempfile.TemporaryDirectory() as tmpdir:
            self.run_software(obj, tmpdir)
            return self.get_energy(os.path.join(tmpdir))

    energies = []

    def known_energy(self, E, params, fail=False):
        if E == -1:
            raise Exception("Invalid calculation")
        for entry in self.energies:
            if entry[1] == E:
                if not fail:
                    print("")
                    print("Clash detected:")
                    print(entry[0])
                    print(params)
                    print("")
                return True

        self.energies.append([params, E])
        return False


class GaussianCalculationTests(CalculationTests):
    def run_software(self, obj, tmpdir):
        with open(os.path.join(tmpdir, "calc.com"), "w") as out:
            out.write(obj.input_file)

        with open(os.path.join(tmpdir, "gaussian.log"), "w") as out:
            ret = subprocess.run(
                shlex.split("g16 calc.com"), cwd=tmpdir, stdout=out, stderr=out
            )

        if ret.returncode != 0:
            print(f"Calculation ended with return code {ret.returncode}")

    def get_energy(self, tmpdir):
        path = os.path.join(tmpdir, "calc.log")
        with open(path) as f:
            lines = f.readlines()
            ind = len(lines) - 1
        while lines[ind].find("SCF Done") == -1:
            ind -= 1
        return float(lines[ind].split()[4])


class OrcaCalculationTests(CalculationTests):
    def run_software(self, obj, tmpdir):
        with open(os.path.join(tmpdir, "calc.inp"), "w") as out:
            out.write(obj.input_file)

        with open(os.path.join(tmpdir, "calc.out"), "w") as out:
            ret = subprocess.run(
                shlex.split("orca calc.inp"), cwd=tmpdir, stdout=out, stderr=out
            )

        if ret.returncode != 0:
            print(f"Calculation ended with return code {ret.returncode}")

    def get_energy(self, tmpdir):
        path = os.path.join(tmpdir, "calc.out")
        with open(path) as f:
            lines = f.readlines()
            ind = len(lines) - 1

        for line in lines:
            print(line)
        while lines[ind].find("FINAL SINGLE POINT ENERGY") == -1:
            ind -= 1

        return float(lines[ind].split()[4].strip())


class XtbCalculationTests(CalculationTests):
    def run_software(self, obj, tmpdir):
        os.system(f"cp {obj.calc.file} {tmpdir}/{obj.get_output_name()}")

        if obj.input_file != "":
            with open(os.path.join(tmpdir, "input"), "w") as out:
                out.write(obj.input_file)

        with open(os.path.join(tmpdir, "calc.out"), "w") as out:
            ret = subprocess.run(
                shlex.split(obj.command), cwd=tmpdir, stdout=out, stderr=out
            )

        if ret.returncode != 0:
            print(f"Calculation ended with return code {ret.returncode}")

    def get_energy(self, tmpdir):
        path = os.path.join(tmpdir, "calc.out")
        with open(path) as f:
            lines = f.readlines()
            ind = len(lines) - 1

        for line in lines:
            print(line)
        while lines[ind].find("TOTAL ENERGY") == -1:
            ind -= 1

        return float(lines[ind].split()[3].strip())

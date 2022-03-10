import os
import shlex
from unittest import TestCase

from ccinput.wrapper import gen_obj
from ccinput.wrapper import gen_input, get_input_from_args, get_parser


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

import os
from unittest import TestCase

from ccinput.wrapper import generate_calculation

class InputTests(TestCase):
    def generate_calculation(self, **params):
        params['in_file'] = os.path.join('/'.join(__file__.split('/')[:-1]), "structures/", params['in_file'])
        return generate_calculation(**params)

    def is_equivalent(self, ref, res):
        ref_lines = [i.strip() for i in ref.strip().split('\n')]
        res_lines = [i.strip() for i in res.strip().split('\n')]

        if len(ref_lines) != len(res_lines):
            print("Different number of lines: {} and {}".format(len(ref_lines), len(res_lines)))
            print("----START REFERENCE---")
            print('\n'.join(ref_lines))
            print("----END REFERENCE---")
            print("----START RESULT---")
            print('\n'.join(res_lines))
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
            if abs(_c1-_c2) > 1e-7:
                return False
        return True


from ccinput.tests.testing_utilities import InputTests
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation


class PysisTests(InputTests):
    def test_xtb_charge_mult(self):
        params = {
            "charge": 1,
            "multiplicity": 2,
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 1
            mult: 2
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts_parallel(self):
        params = {
            "nproc": 8,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 8
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts_do_hess(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "specifications": "tsopt(do_hess=true)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: true
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts_filename(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "specifications": "geom(fn=my_file.xyz)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: my_file.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_internal_coordinates(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "specifications": "geom(type=cart)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: cart
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_accuracy(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "specifications": "calc(acc=0.1)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            acc: 0.1
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts_explicit_vacuum(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "solvent": "vAcUUm",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts_alpb(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "solvent": "dichloromethane",
            "solvation_model": "alpb",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            alpb: ch2cl2
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts_gbsa(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "solvent": "dichloromethane",
            "solvation_model": "gbsa",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            gbsa: ch2cl2
            mult: 1
            pal: 1
            type: xtb
        geom:
            fn: calc.xyz
            type: dlc
        tsopt:
            do_hess: False
            hessian_recalc: 5
            type: rsirfo
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_ts_solvent_requires_model(self):
        params = {
            "nproc": 1,
            "type": "TS Optimization",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "solvent": "dichloromethane",
            "solvation_model": "",
        }

        with self.assertRaises(InvalidParameter):
            self.generate_calculation(**params)

    def test_xtb_gsm_simple(self):
        params = {
            "type": "Minimum Energy Path",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "charge": 0,
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        cos:
            climb: True
            max_nodes: 15
            reparam_every: 2
            reparam_every_full: 3
            type: gs
        geom:
            fn: [mini_ts.xyz, calc2.xyz]
            type: dlc
        opt:
            align: False
            stop_in_when_full: -1
            type: string
        preopt:
            max_cycles: 10
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_gsm_aux_name(self):
        params = {
            "type": "Minimum Energy Path",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "charge": 0,
            "aux_name": "end",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        cos:
            climb: True
            max_nodes: 15
            reparam_every: 2
            reparam_every_full: 3
            type: gs
        geom:
            fn: [mini_ts.xyz, end.xyz]
            type: dlc
        opt:
            align: False
            stop_in_when_full: -1
            type: string
        preopt:
            max_cycles: 10
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_xtb_gsm_node_num(self):
        params = {
            "type": "Minimum Energy Path",
            "file": "mini_ts.xyz",
            "software": "xtb",
            "driver": "pysis",
            "charge": 0,
            "specifications": "cos(max_nodes=21)",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            mult: 1
            pal: 1
            type: xtb
        cos:
            climb: True
            max_nodes: 21
            reparam_every: 2
            reparam_every_full: 3
            type: gs
        geom:
            fn: [mini_ts.xyz, calc2.xyz]
            type: dlc
        opt:
            align: False
            stop_in_when_full: -1
            type: string
        preopt:
            max_cycles: 10
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_orca_gsm(self):
        params = {
            "type": "Minimum Energy Path",
            "file": "mini_ts.xyz",
            "software": "orca",
            "driver": "pysis",
            "charge": 0,
            "method": "M062x",
            "basis_set": "def2svp",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            keywords: M062X Def2-SVP
            mem: 1000
            mult: 1
            pal: 1
            type: orca
        cos:
            climb: True
            max_nodes: 15
            reparam_every: 2
            reparam_every_full: 3
            type: gs
        geom:
            fn: [mini_ts.xyz, calc2.xyz]
            type: dlc
        opt:
            align: False
            stop_in_when_full: -1
            type: string
        preopt:
            max_cycles: 10
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_orca_gsm_cpcm(self):
        params = {
            "type": "Minimum Energy Path",
            "file": "mini_ts.xyz",
            "software": "orca",
            "driver": "pysis",
            "charge": 0,
            "method": "M062x",
            "basis_set": "def2svp",
            "solvent": "DCM",
            "solvation_model": "CPCM",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            charge: 0
            keywords: M062X Def2-SVP CPCM(ch2cl2)
            mem: 1000
            mult: 1
            pal: 1
            type: orca
        cos:
            climb: True
            max_nodes: 15
            reparam_every: 2
            reparam_every_full: 3
            type: gs
        geom:
            fn: [mini_ts.xyz, calc2.xyz]
            type: dlc
        opt:
            align: False
            stop_in_when_full: -1
            type: string
        preopt:
            max_cycles: 10
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_orca_gsm_smd18(self):
        params = {
            "type": "Minimum Energy Path",
            "file": "I.xyz",
            "software": "orca",
            "driver": "pysis",
            "charge": -1,
            "method": "M062x",
            "basis_set": "def2svp",
            "solvent": "DCM",
            "solvation_model": "SMD",
            "solvation_radii": "SMD18",
        }

        inp = self.generate_calculation(**params)

        REF = """
        calc:
            blocks: "%cpcm
            smd true
            SMDsolvent \\"ch2cl2\\"
            radius[53] 2.74
            radius[35] 2.60
            end"
            charge: -1
            keywords: M062X Def2-SVP
            mem: 1000
            mult: 1
            pal: 1
            type: orca
        cos:
            climb: True
            max_nodes: 15
            reparam_every: 2
            reparam_every_full: 3
            type: gs
        geom:
            fn: [I.xyz, calc2.xyz]
            type: dlc
        opt:
            align: False
            stop_in_when_full: -1
            type: string
        preopt:
            max_cycles: 10
        """

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

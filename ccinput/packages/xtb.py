import os

from ccinput.constants import CalcType, ATOMIC_NUMBER, LOWERCASE_ATOMIC_SYMBOLS
from ccinput.utilities import get_solvent
from ccinput.exceptions import InvalidParameter


class XtbCalculation:

    EXECUTABLES = {
        CalcType.OPT: "xtb",
        CalcType.CONSTR_OPT: "xtb",
        CalcType.FREQ: "xtb",
        CalcType.SP: "xtb",
        CalcType.UVVIS_TDA: "stda",
        CalcType.OPTFREQ: "xtb",
        CalcType.CONF_SEARCH: "crest",
        CalcType.CONSTR_CONF_SEARCH: "crest",
    }

    def __init__(self, calc):
        self.calc = calc
        self.program = ""
        self.main_command = ""
        self.cmd_arguments = ""
        self.input_file = ""
        self.specifications = ""
        self.force_constant = 1.0
        self.confirmed_specifications = ""

        self.handle_command()
        self.handle_specifications()

        if self.calc.type == CalcType.CONSTR_CONF_SEARCH:
            self.handle_constraints_crest()
        elif self.calc.type == CalcType.CONSTR_OPT:
            self.handle_constraints_scan()

        self.handle_parameters()

        self.create_command()

    def handle_parameters(self):
        if self.calc.parameters.solvent != "":
            try:
                solvent_keyword = get_solvent(
                    self.calc.parameters.solvent, self.calc.parameters.software
                )
            except KeyError:
                raise InvalidParameter("Invalid solvent")

            if self.calc.parameters.solvation_model == "gbsa":
                self.main_command += f"-g {solvent_keyword} "
            elif self.calc.parameters.solvation_model == "alpb":
                self.main_command += f"--alpb {solvent_keyword} "
            else:
                raise InvalidParameter(
                    "Invalid solvation method for xtb: {}".format(
                        self.calc.parameters.solvation_model
                    )
                )

        if self.calc.charge != 0:
            self.main_command += f"--chrg {self.calc.charge} "

        if self.calc.multiplicity != 1:
            self.main_command += f"--uhf {self.calc.multiplicity} "

    def handle_constraints_scan(self):
        if len(self.calc.constraints) == 0:
            raise InvalidParameter("No constraint in constrained optimisation mode")

        self.input_file += "$constrain\n"
        self.input_file += f"force constant={self.force_constant}\n"
        self.has_scan = False

        for cmd in self.calc.constraints:
            self.input_file += cmd.to_xtb()
            if cmd.scan:
                self.has_scan = True

        if self.has_scan:
            self.input_file += "$scan\n"
            for counter, cmd in enumerate(self.calc.constraints):
                if cmd.scan:
                    self.input_file += (
                        f"{counter+1}: {cmd.start_d}, {cmd.end_d}, {cmd.num_steps}\n"
                    )

    def compress_indices(self, arr):
        comp = []

        def add_to_str(curr):
            if len(curr) == 0:
                return ""
            elif len(curr) == 1:
                return f"{curr[0]}"
            else:
                return f"{curr[0]}-{curr[-1]}"

        _arr = sorted(set(arr))
        curr_atoms = []

        for a in _arr:
            if len(curr_atoms) == 0:
                curr_atoms.append(a)
            else:
                if a == curr_atoms[-1] + 1:
                    curr_atoms.append(a)
                else:
                    comp.append(add_to_str(curr_atoms))
                    curr_atoms = [a]

        comp.append(add_to_str(curr_atoms))
        return ",".join(comp)

    def handle_constraints_crest(self):
        if len(self.calc.constraints) == 0:
            raise InvalidParameter("No constraint in constrained optimisation mode")

        num_atoms = len(self.calc.xyz.split("\n"))
        input_file_name = os.path.basename(self.calc.file)

        self.input_file += "$constrain\n"
        self.input_file += f"force constant={self.force_constant}\n"
        self.input_file += f"reference={input_file_name}\n"
        constr_atoms = []
        for cmd in self.calc.constraints:
            self.input_file += cmd.to_xtb()
            constr_atoms += cmd.ids

        self.input_file += f"atoms: {self.compress_indices(constr_atoms)}\n"

        mtd_atoms = list(range(1, num_atoms))
        for a in constr_atoms:
            if int(a) in mtd_atoms:
                mtd_atoms.remove(int(a))

        self.input_file += "$metadyn\n"
        self.input_file += f"atoms: {self.compress_indices(mtd_atoms)}\n"

    def handle_specifications(self):
        accuracy = -1
        iterations = -1
        method = "gfn2-xtb"
        opt_level = "tight"
        rthr = 0.6
        ewin = 6
        cmd_arguments = ""

        ALLOWED = "qwertyuiopasdfghjklzxcvbnm-1234567890./= "
        clean_specs = "".join(
            [
                i
                for i in self.specifications
                + self.calc.parameters.specifications.lower()
                if i in ALLOWED
            ]
        )
        clean_specs = clean_specs.replace("=", " ").replace("  ", " ")

        specs = clean_specs.strip().split("--")

        for spec in specs:
            if spec.strip() == "":
                continue
            ss = spec.strip().split()
            if len(ss) == 1:
                if ss[0] in ["gfn2", "gfn1", "gfn0", "gfnff", "gfn2//gfnff"]:
                    if ss[0] == "gfn2//gfnff" and self.calc.type not in [
                        CalcType.CONF_SEARCH,
                        CalcType.CONSTR_CONF_SEARCH,
                    ]:
                        raise InvalidParameter(
                            f"Invalid method for calculation type: {ss[0]}"
                        )
                    if ss[0] in ["gfn2", "gfn1", "gfn0"]:
                        method = f"{ss[0][:-1]} {ss[0][-1]}"
                    else:
                        method = ss[0]
                elif ss[0] == "nci":
                    self.cmd_arguments += "--nci "
                elif ss[0] == "quick":
                    self.cmd_arguments += "--quick "
                elif ss[0] == "squick":
                    self.cmd_arguments += "--squick "
                elif ss[0] == "mquick":
                    self.cmd_arguments += "--mquick "
                else:
                    raise InvalidParameter("Invalid specification")
            elif len(ss) == 2:
                if ss[0] == "o" or ss[0] == "opt":
                    if ss[1] not in [
                        "crude",
                        "sloppy",
                        "loose",
                        "lax",
                        "normal",
                        "tight",
                        "vtight",
                        "extreme",
                    ]:
                        raise InvalidParameter("Invalid optimization specification")
                    opt_level = ss[1]
                elif ss[0] == "rthr":
                    if self.calc.type not in [
                        CalcType.CONF_SEARCH,
                        CalcType.CONSTR_CONF_SEARCH,
                    ]:
                        raise InvalidParameter(
                            "Invalid specification for calculation type: rthr"
                        )
                    try:
                        float(ss[1])
                    except ValueError:
                        raise InvalidParameter(
                            "Parameter rthr must be a floating point value"
                        )
                    else:
                        rthr = ss[1]

                elif ss[0] == "ewin":
                    if self.calc.type not in [
                        CalcType.CONF_SEARCH,
                        CalcType.CONSTR_CONF_SEARCH,
                    ]:
                        raise InvalidParameter(
                            "Invalid specification for calculation type: ewin"
                        )
                    try:
                        float(ss[1])
                    except ValueError:
                        raise InvalidParameter(
                            "Parameter ewin must be a floating point value"
                        )
                    else:
                        ewin = ss[1]
                elif ss[0] == "acc":
                    accuracy = float(ss[1])
                elif ss[0] == "iterations":
                    try:
                        iterations = int(ss[1])
                    except ValueError:
                        raise InvalidParameter(
                            "Invalid number of iterations: must be an integer"
                        )
                elif ss[0] == "forceconstant":
                    try:
                        self.force_constant = float(ss[1])
                    except ValueError:
                        raise InvalidParameter(
                            "Invalid force constant: must be a floating point value"
                        )
                elif ss[0] == "gfn":
                    if ss[1] not in ["0", "1", "2"]:
                        raise InvalidParameter("Invalid GFN version")
                    method = f"{ss[0]} {ss[1]}"
                else:
                    raise InvalidParameter(f"Unknown specification: {ss[0]}")
            else:
                raise InvalidParameter(f"Invalid specification: {ss}")

        if accuracy != -1:
            self.cmd_arguments += f"--acc {accuracy:.2f} "
        if iterations != -1:
            self.cmd_arguments += f"--iterations {iterations} "
        if method != "gfn2-xtb" and method != "gfn 2":
            self.cmd_arguments += f"--{method} "
        if opt_level != "normal":
            self.main_command = self.main_command.replace(
                "--opt ", f"--opt {opt_level} "
            )
            self.confirmed_specifications += f"--opt {opt_level} "

        if self.calc.type in [CalcType.CONF_SEARCH, CalcType.CONSTR_CONF_SEARCH]:
            self.cmd_arguments += f"--rthr {rthr} --ewin {ewin}"
            self.confirmed_specifications += self.cmd_arguments.strip()

            self.cmd_arguments = self.cmd_arguments.replace(
                "--", "-"
            )  # Crest 2.10.2 does not read arguments with double dashes
        else:
            self.confirmed_specifications += self.cmd_arguments.strip()

    def handle_command(self):
        self.program = self.EXECUTABLES[self.calc.type]

        if self.calc.type == CalcType.OPT:
            self.specifications = "--opt tight "
            self.main_command += "--opt "
        elif self.calc.type == CalcType.OPTFREQ:
            self.main_command = "--ohess "  # Not sure if the tightness will be parsed
        elif self.calc.type == CalcType.CONSTR_CONF_SEARCH:
            self.main_command += "-cinp input "
        elif self.calc.type == CalcType.CONSTR_OPT:
            self.main_command += "--opt --input input "
        elif self.calc.type == CalcType.FREQ:
            self.main_command += "--hess "

    def create_command(self):
        if self.calc.file:
            input_file_name = os.path.basename(self.calc.file)
        else:
            input_file_name = self.calc.name + ".xyz"

        if self.calc.type in [CalcType.CONF_SEARCH, CalcType.CONSTR_CONF_SEARCH]:
            self.main_command = self.main_command.replace(
                "--", "-"
            )  # Crest 2.10.2 does not read arguments with double dashes

        if self.main_command != "":
            self.command = f"{self.program} {input_file_name} {self.main_command.strip()} {self.cmd_arguments}".strip()
        else:
            self.command = (
                f"{self.program} {input_file_name} {self.cmd_arguments}".strip()
            )

    @property
    def output(self):
        if self.input_file == "":
            return self.command
        else:
            return f"{self.command}\n\n---input:\n{self.input_file}"

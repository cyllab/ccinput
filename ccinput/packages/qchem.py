from typing import Dict, List

from ccinput.calculation import Calculation
from ccinput.constants import (
    CalcType,
    ATOMIC_NUMBER,
    LOWERCASE_ATOMIC_SYMBOLS,
    THEORY_LEVELS,
)
from ccinput.utilities import clean_xyz, get_basis_set, get_method, get_solvent
from ccinput.exceptions import InvalidParameter

_CALC_TYPE_TO_JOBTYPE = {
    CalcType.SP: "sp",
    CalcType.OPT: "opt",
    # CalcType.CONSTR_OPT: [],
    CalcType.FREQ: "freq",
    CalcType.TS: "ts",
    CalcType.MEP: "fsm",
    # Go through moprop rather than d-scf
    CalcType.NMR: "sp",
    CalcType.UVVIS: "sp",
    CalcType.UVVIS_TDA: "sp",
    # plots?
    # CalcType.MO: [],
    # does combined exist?
    # CalcType.OPTFREQ: [],
}

_TEMPLATE = """$comment
{comment}
$end

$molecule
{charge} {multiplicity}
{xyz_structure}$end

$rem
{rem_lines}
$end
{other_blocks}
"""


class QChemCalculation:
    def __init__(self, calc: Calculation) -> None:
        self.calc = calc
        self.other_blocks: List[str] = []
        self.solvation_radii: Dict[int, float] = {}

        self._rems: Dict[str, str] = {}

        self.xyz_structure = ""
        self.input_file = ""

        self._rems["jobtype"] = _CALC_TYPE_TO_JOBTYPE[self.calc.type]
        if self.calc.type == CalcType.NMR:
            self._rems["moprop"] = "nmr"
        self._rems["method"] = get_method(self.calc.parameters.method, "qchem")
        if self._rems["method"] not in THEORY_LEVELS["semiempirical"]:
            self._rems["basis"] = get_basis_set(self.calc.parameters.basis_set, "qchem")
        else:
            self._rems["semi_empirical"] = "true"
        self._rems["mem_total"] = str(self.calc.mem)

        if self.calc.parameters.d3:
            self._rems["dft_d"] = "d3_zero"
        elif self.calc.parameters.d3bj:
            self._rems["dft_d"] = "d3_bj"

        self.handle_specifications()

        # self.handle_command()
        self.handle_xyz()

        self.parse_custom_solvation_radii()
        self.handle_solvation()

        self.create_input_file()

    @property
    def confirmed_specifications(self):
        pass

    def clean(self, s: str) -> str:
        WHITELIST = set(
            "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. "
        )
        return "".join([c for c in s if c in WHITELIST])

    def add_option(self, key, option) -> None:
        pass

    def adds_option(self, key, options) -> None:
        pass

    def add_commands(self, commands) -> None:
        pass

    def handle_specifications(self) -> None:
        _specifications = self.clean(self.calc.parameters.specifications.lower())

        if _specifications != "":
            sspecs = _specifications.split()
            ind = 0
            while ind < len(sspecs):
                spec = sspecs[ind]
                if spec == "--phirshfeld":
                    self._rems["hirshfeld"] = "true"
                ind += 1

    def handle_command(self):
        ctype = self.calc.type
        if ctype == CalcType.MEP:
            raise RuntimeError
        else:
            raise RuntimeError

    def parse_custom_basis_set(self, base_bs):
        pass

    def handle_custom_basis_sets(self) -> None:
        if len(self.calc.parameters.custom_basis_sets) == 0:
            return

        unique_atoms = _get_unique_atoms(self.calc.xyz)

        for el, bs_keyword in self.calc.parameters.custom_basis_sets.items():
            pass

    def handle_xyz(self) -> None:
        lines = [i + "\n" for i in clean_xyz(self.calc.xyz).split("\n") if i != ""]
        self.xyz_structure = "".join(lines)

    def parse_custom_solvation_radii(self) -> None:
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

    def get_radii_specs(self) -> str:
        radii_specs = []
        if self.solvation_radii:
            radii_specs.append("$van_der_waals")
            radii_specs.append("1")
            for el, rad in self.solvation_radii.items():
                radii_specs.append(f"{el} {rad:.2f}")
            radii_specs.append("$end")
        return "\n".join(radii_specs)

    def handle_solvation(self) -> None:
        if self.calc.parameters.solvent.lower() not in ["vacuum", ""]:
            solvent_keyword = get_solvent(
                self.calc.parameters.solvent,
                self.calc.parameters.software,
                solvation_model=self.calc.parameters.solvation_model,
            )
            model = self.calc.parameters.solvation_model
            radii_set = self.calc.parameters.solvation_radii
            custom_radii = self.calc.parameters.custom_solvation_radii

            if radii_set not in _ALLOWED_RADII_SETS:
                raise InvalidParameter(
                    "invalid argument for solvent radii type: "
                    f"{radii_set} not in {_ALLOWED_RADII_SETS}"
                )

            radii_specs = self.get_radii_specs()
            model_block = []
            solv_block = []

            if model == "kirkwood":
                raise NotImplementedError(
                    f"solvent_method = {model} not implemented yet"
                )
            elif "pcm" in model:
                pcm_model = model
                model = "pcm"
                model_block = ["$pcm", f"theory {pcm_model}"]
                if radii_set != "":
                    model_block.append(f"radii {radii_set}")
                if radii_specs:
                    if radii_set not in ("", "read"):
                        raise InvalidParameter(
                            f"custom radii set but not requested to read: {radii_set}"
                        )
                    if radii_set == "":
                        model_block.append("radii read")
                model_block.append("$end")
                solv_block = ["$solvent", f"solventname {solvent_keyword}", "$end"]
            elif model == "isosvp":
                raise NotImplementedError(
                    f"solvent_method = {model} not implemented yet"
                )
            elif model == "cosmo":
                raise NotImplementedError(
                    f"solvent_method = {model} not implemented yet"
                )
            elif model == "sm8":
                raise NotImplementedError(
                    f"solvent_method = {model} not implemented yet"
                )
            elif model == "sm12":
                raise NotImplementedError(
                    f"solvent_method = {model} not implemented yet"
                )
            elif model == "smd":
                solv_block = ["$smx", f"solvent {solvent_keyword}", "$end"]
            elif model == "chem_sol":
                raise NotImplementedError(
                    f"solvent_method = {model} not implemented yet"
                )
            elif model == "peqs":
                raise NotImplementedError(
                    f"solvent_method = {model} not implemented yet"
                )
            else:
                raise InvalidParameter(
                    f"Invalid solvation model for Q-Chem: "
                    f"'{self.calc.parameters.solvation_model}'"
                )

            self._rems["solvent_method"] = model

            if model_block:
                self.other_blocks.append("")
                self.other_blocks.extend(model_block)
            self.other_blocks.append("")
            self.other_blocks.extend(solv_block)
            if radii_specs:
                self.other_blocks.append("")
                self.other_blocks.append(radii_specs)

    def create_input_file(self) -> None:
        # breakpoint()
        self.input_file = _TEMPLATE.format(
            comment=self.calc.header,
            rem_lines="\n".join(f"{k} = {v}" for k, v in self._rems.items()),
            charge=self.calc.charge,
            multiplicity=self.calc.multiplicity,
            xyz_structure=self.xyz_structure,
            other_blocks="\n".join(self.other_blocks),
        )

    @property
    def output(self):
        return self.input_file


def _get_unique_atoms(xyz: str) -> List[str]:
    unique_atoms = []
    for line in xyz.split("\n"):
        if line.strip() == "":
            continue
        a, *_ = line.split()
        if a not in unique_atoms:
            unique_atoms.append(a)
    return unique_atoms


_ALLOWED_RADII_SETS = ("", "bondi", "ff", "read")

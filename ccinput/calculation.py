import hashlib
from itertools import zip_longest

from ccinput.exceptions import InvalidParameter, InternalError, ImpossibleCalculation, \
                               MissingParameter
from ccinput.utilities import get_abs_software, get_method, get_abs_basis_set, \
                              get_abs_solvent, get_theory_level, standardize_memory, \
                              get_npxyz, get_coord, has_dispersion_parameters
from ccinput.constants import ATOMIC_NUMBER, SYN_SOFTWARE
from ccinput.logging import warn

class Calculation:
    """
        Holds all the data required to generate an input file.
        Its fields are the parameters likely to change (charge, multiplicity,
        xyz, calculation type...). The other parameters are contained in the
        Parameters class (accessed through self.parameters).
    """
    def __init__(self, xyz, parameters, type, constraints=[], nproc=0, mem=0,
            charge=0, multiplicity=1, aux_name="calc2", name="calc",
            header="File created by ccinput", software=""):
        self.xyz = xyz
        self.parameters = parameters
        self.type = type

        if software not in SYN_SOFTWARE:
            raise InvalidParameter(f"Invalid software: '{software}'")

        if nproc == 0:
            raise MissingParameter("Number of cores unspecified")

        try:
            self.nproc = int(nproc)
        except ValueError:
            raise InvalidParameter(f"Invalid number of cores: '{nproc}'")

        if abs(self.nproc - float(nproc)) > 1e-4:
            raise InvalidParameter(f"The number of cores must be an integer (received '{nproc}')")

        if self.nproc < 1:
            raise InvalidParameter(f"The number of cores must at least 1 (received '{nproc}')")

        try:
            self.mem = standardize_memory(mem)
        except InvalidParameter:
            if software in ["gaussian"]:
                raise

        try:
            self.charge = int(charge)
        except ValueError:
            raise InvalidParameter(f"Invalid charge: '{charge}'")

        if abs(self.charge - float(charge)) > 1e-4:
            raise InvalidParameter(f"Charge must be an integer (received '{charge}')")

        try:
            self.multiplicity = int(multiplicity)
        except ValueError:
            raise InvalidParameter(f"Invalid multiplicity: '{multiplicity}'")

        if abs(self.multiplicity - float(multiplicity)) > 1e-4:
            raise InvalidParameter(f"Multiplicity must be an integer (received '{multiplicity}')")
        if self.multiplicity < 1:
            raise InvalidParameter(f"Multiplicity must at least 1 (received '{multiplicity}')")
        self.verify_charge_mult()
        self.constraints = constraints

        self.name = name
        self.aux_name = aux_name
        self.header = header

    def verify_charge_mult(self):
        """ Verifies that the requested charge and multiplicity are possible for the structure """
        electrons = 0
        for line in self.xyz.split('\n'):
            if line.strip() == '':
                continue
            el = line.split()[0]
            electrons += ATOMIC_NUMBER[el]

        electrons -= self.charge
        odd_e = electrons % 2
        odd_m = self.multiplicity % 2

        if odd_e == odd_m:
            raise ImpossibleCalculation(f"This combination of charge ({self.charge}) " +
                    f"and multiplicity ({self.multiplicity}) is impossible")

class Parameters:
    """
        Holds all the parameters about the computational method.
        These parameters do not depend on the particular system
        (e.g. regarding the charge and multiplicity) and can be reused.
    """

    def __init__(self, software, solvent="", solvation_model="", solvation_radii="",
            custom_solvation_radii="", basis_set="", method="", specifications="",
            density_fitting="", custom_basis_sets="", d3=False, d3bj=False, **kwargs):

        self.solvent = get_abs_solvent(solvent)
        if self.solvent != "":
            self.solvation_model = solvation_model.lower()
            self.solvation_radii = solvation_radii.lower()
            self.custom_solvation_radii = custom_solvation_radii.lower()

            if self.solvation_model.strip() == "":
                raise InvalidParameter("No solvation model specified, " +
                            "although solvation is requested")
            if self.solvation_radii.strip() == "":
                warn("No solvation radii specified; using default radii")
        else:
            self.solvation_model = ""
            self.solvation_radii = ""
            self.custom_solvation_radii = ""

        self.software = get_abs_software(software)

        if method == "":
            if 'functional' in kwargs:
                method = get_method(kwargs['functional'], self.software)
            else:
                raise InvalidParameter("No calculation method specified (method='...')")
        else:
            self.method = get_method(method, self.software)

        self.theory_level = get_theory_level(method)

        if self.theory_level not in ["semi-empirical", "xtb", "special"]:
            self.basis_set = get_abs_basis_set(basis_set)
        else:
            self.basis_set = ""

        if d3 and d3bj:
            raise InvalidParameter("Cannot use both D3(0) and D3BJ dispersion corrections")

        self.d3 = d3
        self.d3bj = d3bj

        if d3 and not has_dispersion_parameters(self.method, version='d3'):
            warn(f"Your calculation requests a method ({self.method}) that may not have D3 "
                 "parameters. Be aware that the calculation might result in an error")

        if d3bj and not has_dispersion_parameters(self.method, version='d3bj'):
            warn(f"Your calculation requests a method ({self.method}) that may not have D3BJ "
                 "parameters. Be aware that the calculation might result in an error")

        self.specifications = specifications
        self.density_fitting = density_fitting
        self.custom_basis_sets = custom_basis_sets
        self.kwargs = kwargs

    def __eq__(self, other):
        values = [(k,v) for k,v in self.__dict__.items()]
        other_values = [(k,v) for k,v in other.__dict__.items()]

        return values == other_values

    @property
    def md5(self):
        """
            Returns a hash digest of the parameters to easily compare
            parameters without using objects.
        """
        values = [(k,v) for k,v in self.__dict__.items()]
        params_str = ""
        for k, v in values:
            if isinstance(v, int):
                params_str += f"{k}={v};"
            elif isinstance(v, str):
                params_str += f"{k}={v.lower()};"
            else:
                raise InternalError("Unknown value type")
        hash = hashlib.md5(bytes(params_str, 'UTF-8'))
        return hash.digest()

class Constraint:
    """
        Class to contain a single constraint (freeze, scan and the like) and provide useful functions for it.
    """
    def __init__(self, scan=False, start_d=None, end_d=None, step_size=None,
                 num_steps=None, ids=[], xyz=None, software=""):
        self.scan = scan
        self.ids = ids # One-indexed (Gaussian like, not like ORCA)
        self.software = software

        if len(ids) < 2 or len(ids) > 4:
            raise InvalidParameter(f"Invalid number of atoms: {len(ids)}, needs to be between 2 and 4")

        _xyz = get_npxyz(xyz)

        if max(ids) > len(_xyz):
            raise InvalidParameter(f"Invalid atom index {max(ids)}: larger than the number of atoms")
        if len(ids) != len(set(ids)):
            raise InvalidParameter(f"The provided atom indices are not all unique")

        if self.scan:
            num_params = 0
            for p in [end_d, step_size, num_steps]:
                if p is not None:
                    num_params += 1

            if num_params < 2:
                raise InvalidParameter("Not enough constraint parameters given")
            elif num_params == 3:
                raise InvalidParameter("Too many constraint parameters given, more than one possibility of unique parameters")

            if start_d is not None and end_d is not None and xyz is None:
                raise InvalidParameter("The XYZ structure needs to be specified in order to calculate the initial coordinate value")

        if software != 'gaussian':
            self.start_d = start_d
        else:
            if start_d is not None:
                warn("Gaussian only allows scans that start from the current structure; " +
                     "overriding the specified initial value")
            self.start_d = None

        self.end_d = end_d
        self.step_size = step_size
        self.num_steps = num_steps

        if self.start_d is None:
            self.start_d = get_coord(_xyz, self.ids)

        if self.scan:
            self.complete_parameters()

    def complete_parameters(self):
        if self.end_d is None:
            assert self.step_size is not None
            assert self.num_steps is not None

            self.end_d = self.start_d + self.num_steps*self.step_size

        if self.step_size is None:
            assert self.end_d is not None
            assert self.num_steps is not None

            self.step_size = round((self.end_d-self.start_d)/self.num_steps, 2)

        if self.num_steps is None:
            assert self.end_d is not None
            assert self.step_size is not None

            # Make sure the number of steps is positive overall
            # Also use the absolute step size here, as the start and end are defined
            # An mathematically incorrect sign won't matter, as the step size is recalculated
            self.num_steps = abs(int((self.end_d-self.start_d)/abs(self.step_size)))

            # Adjust the step size to end up exactly at the end point
            self.step_size = round((self.end_d-self.start_d)/self.num_steps, 2)

    def to_orca(self):
        ids_str = ' '.join([str(i-1) for i in self.ids])
        type = len(self.ids)

        if type == 2:
            t = "B"
        elif type == 3:
            t = "A"
        elif type == 4:
            t = "D"

        if self.scan:
            return f"{t} {ids_str} = {self.start_d:.2f}, {self.end_d:.2f}, {self.num_steps}\n"
        else:
            return f"{{ {t} {ids_str} C }}\n"

    def to_gaussian(self):
        ids_str = ' '.join([str(i) for i in self.ids])
        type = len(self.ids)

        if type == 2:
            t = "B"
        elif type == 3:
            t = "A"
        elif type == 4:
            t = "D"
            if self.scan:
                # Different sign convention used by Gaussian for dihedral angles (?)
                self.step_size *= -1

        if self.scan:
            return f"{t} {ids_str} S {self.num_steps} {self.step_size}\n"
        else:
            return f"{t} {ids_str} F\n"


def parse_freeze_constraints(arr, xyz_str, software=""):
    if len(arr) == 0:
        return []
    constr = ""
    for c in arr:
        constr += f'Freeze/{"_".join([str(i) for i in c])};'

    return parse_str_constraints(constr, xyz_str, software=software)

def parse_scan_constraints(arr, sfrom, sto, snsteps, sstep, xyz_str, software=""):
    if len(arr) == 0:
        return []

    scans = []
    for ids, fro, to, nsteps, step in zip_longest(arr, sfrom, sto, snsteps, sstep):
        if ids is None:
            raise InvalidParameter("Not enough sets of atom indices specified for the number of other parameters")
        _ids = [int(i) for i in ids]
        scans.append(gen_constraint(_ids, xyz_str, 'scan', start_str=fro, end_str=to, nsteps_str=nsteps, step_str=step, software=software))

    return scans

def parse_str_constraints(s, xyz_str, software=""):
    if s.strip() == "":
        return []

    if software == "":
        warn("No software specified for the constraints; the behaviour might be incorrect")

    constraints = []
    cs = s.split(';')
    for c in cs:
        if c.strip() == '':
            continue

        if c.count('/') != 1:
            raise InvalidParameter(f"Invalid constraint string: '{c}'")

        specs_str, ids_str = c.lower().split('/')

        try:
            ids = [int(i) for i in ids_str.split('_')]
        except ValueError:
            raise InvalidParameter(f"Could not parse the atom numbers from the string '{ids_str}'")

        constraints.append(gen_constraint(ids, xyz_str, *specs_str.split('_'), software=software))

    return constraints

def gen_constraint(ids, xyz_str, option, start_str=None, end_str=None, nsteps_str=None,
                   step_str=None, software=""):
    """
        Generate a constraint object from arrays of parameters.
        In case of multiple values in the arrays, each parameter with the same index
        is assumed to be of the same command constraint.

        The atom indices must already be integers in ids as array of arrays of ints.

    """
    t = len(ids)

    if option not in ['scan', 'freeze']:
        raise InvalidParameter(f"Invalid type of scan: '{specs[0]}'")

    if option == 'scan':
        if start_str is not None:
            try:
                start_d = float(start_str)
            except ValueError:
                raise InvalidParameter(f"Invalid initial value: '{start_str}'")

            if t == 2 and start_d < 0.01:
                raise InvalidParameter(f"Invalid initial distance: '{start_str}'")
        else:
            start_d = None

        if end_str is not None:
            try:
                end_d = float(end_str)
            except ValueError:
                raise InvalidParameter(f"Invalid final value: '{end_str}'")

            if t == 2 and end_d < 0.01:
                raise InvalidParameter(f"Invalid final distance: '{end_str}'")
        else:
            end_d = None

        if nsteps_str is not None:
            try:
                num_steps = int(nsteps_str)
            except ValueError:
                raise InvalidParameter(f"Invalid number of steps: '{nsteps_str}'")

            if num_steps < 1:
                raise InvalidParameter(f"Invalid number of steps: '{nsteps_str}'")
        else:
            num_steps = None

        if step_str is not None:
            try:
                step_size = float(step_str)
            except ValueError:
                raise InvalidParameter(f"Invalid final value: '{step_str}'")
        else:
            step_size = None

    else:
        start_d = None
        end_d = None
        num_steps = None
        step_size = None

    return Constraint(scan=option == 'scan', start_d=start_d, end_d=end_d, num_steps=num_steps,
                      step_size=step_size, ids=ids, xyz=xyz_str, software=software.lower())

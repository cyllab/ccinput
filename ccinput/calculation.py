import hashlib

from ccinput.exceptions import *
from ccinput.utilities import *
from ccinput.constants import ATOMIC_NUMBER, SYN_SOFTWARE
from ccinput.logging import warn

class Calculation:
    """
        Holds all the data required to generate an input file. Its fields are the parameters likely to change (charge, multiplicity, xyz, calculation type...). The other parameters are contained in the Parameters class (accessed through self.parameters).
    """
    def __init__(self, xyz, parameters, type, constraints="", nproc=0, mem=0, charge=0, multiplicity=1, name="calc", header="File created by ccinput", software=""):
        self.xyz = xyz
        self.parameters = parameters
        self.type = type

        if software not in SYN_SOFTWARE.keys():
            raise InvalidParameter("Invalid software: '{}'".format(software))

        if nproc == 0:
            raise MissingParameter("Number of cores unspecified")

        try:
            self.nproc = int(nproc)
        except ValueError:
            raise InvalidParameter("Invalid number of cores: '{}'".format(nproc))

        if abs(self.nproc - float(nproc)) > 1e-4:
            raise InvalidParameter("The number of cores must be an integer (received '{}')".format(nproc))

        if self.nproc < 1:
            raise InvalidParameter("The number of cores must at least 1 (received '{}')".format(nproc))

        try:
            self.mem = standardize_memory(mem)
        except InvalidParameter:
            if software in ["gaussian"]:
                raise

        try:
            self.charge = int(charge)
        except ValueError:
            raise InvalidParameter("Invalid charge: '{}'".format(charge))

        if abs(self.charge - float(charge)) > 1e-4:
            raise InvalidParameter("Charge must be an integer (received '{}')".format(charge))

        try:
            self.multiplicity = int(multiplicity)
        except ValueError:
            raise InvalidParameter("Invalid multiplicity: '{}'".format(multiplicity))

        if abs(self.multiplicity - float(multiplicity)) > 1e-4:
            raise InvalidParameter("Multiplicity must be an integer (received '{}')".format(multiplicity))
        if self.multiplicity < 1:
            raise InvalidParameter("Multiplicity must at least 1 (received '{}')".format(multiplicity))
        self.verify_charge_mult()
        self.constraints = constraints

        self.name = name
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
            raise ImpossibleCalculation("This combination of charge ({}) and multiplicity ({}) is impossible".format(self.charge, self.multiplicity))

class Parameters:
    """
        Holds all the parameters about the computational method. These parameters do not depend on the particular system (e.g. regarding the charge and multiplicity) and can be reused.
    """

    def __init__(self, software, solvent="", solvation_model="", solvation_radii="", basis_set="", method="", specifications="", density_fitting="", custom_basis_sets="", **kwargs):

        if solvent.strip() != "":
            self.solvent = get_abs_solvent(solvent)
            self.solvation_model = solvation_model.lower()
            self.solvation_radii = solvation_radii.lower()

            if self.solvation_model.strip() == "":
                raise InvalidParameter("No solvation model specified, although solvation is requested")
            if self.solvation_radii.strip() == "":
                warn("No solvation radii specified; using default radii")
        else:
            self.solvent = ""
            self.solvation_model = ""
            self.solvation_radii = ""

        self.software = get_abs_software(software)

        if method == "":
            if 'functional' in kwargs.keys():
                method = get_abs_method(kwargs['functional'])
            else:
                raise InvalidParameter("No calculation method specified (method='...')")
        else:
            self.method = get_abs_method(method)

        self.theory_level = get_theory_level(method)

        if self.theory_level not in ["semi-empirical", "xtb", "special"]:
            self.basis_set = get_abs_basis_set(basis_set)
        else:
            self.basis_set = ""

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
        """ Returns a hash digest of the parameters to easily compare parameters without using objects """
        values = [(k,v) for k,v in self.__dict__.items()]
        params_str = ""
        for k, v in values:
            if isinstance(v, int):
                params_str += "{}={};".format(k, v)
            elif isinstance(v, str):
                params_str += "{}={};".format(k, v.lower())
            else:
                raise InternalError("Unknown value type")
        hash = hashlib.md5(bytes(params_str, 'UTF-8'))
        return hash.digest()


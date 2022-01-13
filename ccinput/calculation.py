import hashlib

from ccinput.exceptions import *
from ccinput.utilities import *
from ccinput.constants import ATOMIC_NUMBER

class Calculation:
    def __init__(self, xyz, parameters, type, constraints="", nproc=0, mem=0, charge=0, multiplicity=1, name="calc", header="File created by ccinput"):
        self.xyz = xyz
        self.parameters = parameters
        self.type = type

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

        self.mem = standardize_memory(mem)

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
    def __init__(self, software, solvent="", solvation_model="", solvation_radii="", basis_set="", method="", specifications="", density_fitting="", custom_basis_sets="", **kwargs):
        if solvent.strip() != "":
            self.solvent = get_abs_solvent(solvent)
        else:
            self.solvent = ""

        self.solvation_model = solvation_model
        self.solvation_radii = solvation_radii
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

    # save as json or something

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


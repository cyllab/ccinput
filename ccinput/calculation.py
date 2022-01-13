import hashlib

from ccinput.exceptions import *
from ccinput.utilities import standardize_memory

class Calculation:
    def __init__(self, xyz, parameters, type, constraints="", nproc=0, mem=0, charge=0, multiplicity=1):
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

        self.constraints = constraints


class Parameters:
    def __init__(self, software, solvent="", solvation_model="", solvation_radii="", basis_set="", theory_level="", method="", specifications="", density_fitting="", custom_basis_sets="", **kwargs):
        self.solvent = solvent
        self.solvation_model = solvation_model
        self.solvation_radii = solvation_radii
        self.software = software
        self.basis_set = basis_set
        self.theory_level = theory_level # necessary?
        self.method = method
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


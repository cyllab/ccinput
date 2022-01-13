import hashlib

class Calculation:
    def __init__(self, xyz, parameters, type, constraints="", nproc=0, mem=0, charge=0, multiplicity=1):
        self.xyz = xyz
        self.parameters = parameters
        self.type = type

        #checks

        self.nproc = nproc
        self.mem = mem # not required explicitly in all packages... to check

        try:
            self.charge = int(charge)
        except ValueError:
            raise Exception("Invalid charge: {}".format(charge))

        try:
            self.multiplicity = int(multiplicity) #>= 1
        except ValueError:
            raise Exception("Invalid charge: {}".format(charge))

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
                raise Exception("Unknown value type")
        hash = hashlib.md5(bytes(params_str, 'UTF-8'))
        return hash.digest()


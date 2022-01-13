import string
import numpy as np

from ccinput.constants import *

def clean(s):
	WHITELIST = set("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ/()=-,. ")
	return ''.join([c for c in s if c in WHITELIST])

def get_abs_type(str_type):
    """
        Converts a string calculation type into the correct CalcType.
        Takes into account different equivalent ways to write the calculation types
    """
    _str_type = str_type.lower()
    if _str_type in STR_TYPES.keys():
        return STR_TYPES[_str_type]
    else:
        raise Exception("Invalid calculation type: {}".format(str_type))

def get_abs_software(software):
    _software = software.lower().strip()
    for s in SYN_SOFTWARE.keys():
        if _software in SYN_SOFTWARE[s] or _software == s:
            return s
    return -1

def get_abs_method(method):
    for m in SYN_METHODS.keys():
        if method.lower() in SYN_METHODS[m] or method.lower() == m:
            return m
    return -1

def get_abs_basis_set(basis_set):
    for bs in SYN_BASIS_SETS.keys():
        if basis_set.lower() in SYN_BASIS_SETS[bs] or basis_set.lower() == bs:
            return bs
    return -1

def get_abs_solvent(solvent):
    for solv in SYN_SOLVENTS.keys():
        if solvent.lower() in SYN_SOLVENTS[solv] or solvent.lower() == solv:
            return solv
    return -1

def get_method(method, software):
    abs_method = get_abs_method(method)
    if abs_method == -1:
        return method
    return SOFTWARE_METHODS[software][abs_method]

def get_basis_set(basis_set, software):
    abs_basis_set = get_abs_basis_set(basis_set)
    if abs_basis_set == -1:
        return basis_set
    return SOFTWARE_BASIS_SETS[software][abs_basis_set]

def get_solvent(solvent, software, solvation_model="SMD"):
    abs_solvent = get_abs_solvent(solvent)

    if abs_solvent == -1:
        return solvent

    if software == "ORCA" and abs_solvent == "n-octanol":
        #Weird exception in ORCA
        if solvation_model == "SMD":
            return "1-octanol"
        elif solvation_model == "CPCM":
            return "octanol"
        #Note that ch2cl2 is a valid keyword for SMD, although not listed in the manual

    return SOFTWARE_SOLVENTS[software][abs_solvent]

def clean_xyz(xyz):
    return ''.join([x if x in string.printable else ' ' for x in xyz])

def get_distance(xyz, a, b):
    return np.linalg.norm(xyz[a-1][1] - xyz[b-1][1])

def get_angle(xyz, a, b, c):
    v1 = xyz[a-1][1] - xyz[b-1][1]
    v2 = xyz[c-1][1] - xyz[b-1][1]

    return np.arccos(v1.dot(v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))*180/np.pi

def get_dihedral(xyz, a, b, c, d):
    v1 = xyz[b-1][1] - xyz[a-1][1]
    v2 = xyz[c-1][1] - xyz[b-1][1]
    v3 = xyz[d-1][1] - xyz[c-1][1]

    n1 = np.cross(v1, v2)
    n1 = n1/np.linalg.norm(n1)

    n2 = np.cross(v2, v3)
    n2 = n2/np.linalg.norm(n2)

    m1 = np.cross(n1, v2/np.linalg.norm(v2))
    x = n1.dot(n2)
    y = m1.dot(n2)

    return np.arctan2(y, x)*180/np.pi


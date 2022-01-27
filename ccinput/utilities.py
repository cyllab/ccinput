import os
import string
import numpy as np

from ccinput.constants import *
from ccinput.exceptions import *
from ccinput.logging import warn

MEMORY_FACTORS = {
            'm': 1,
            'mb': 1,
            'mib': 1.048576, # (1024/1000)^2
            'g': 1000,
            'gb': 1000,
            'gib': 1073.741824,# 1000*(1024/1000)^3
            't': 1000000,
            'tb': 1000000,
            'tib': 1099511.628,# 1000^2*(1024/1000)^4
        }

def verify_memory_valid(mem):
    if mem < 0:
        raise InvalidParameter(f"The amount of memory must be positive, not '{mem}'")

    if mem > 1e8:
        raise InvalidParameter(f"Unreasonable amount of memory requested: '{mem}' MB")

    return round(mem)

def standardize_memory(mem):
    """ Converts a string specifying an amount of memory to an integer number of megabytes """

    if isinstance(mem, (int, float)):
        return verify_memory_valid(mem)

    if isinstance(mem, str):
        ind = 0
        for ind, c in enumerate(mem.lower()):
            if c in string.ascii_lowercase:
                break
        else:
            ind += 1

        _val = mem[:ind]
        unit = mem[ind:].lower()

        if len(_val) == 0:
            raise InvalidParameter(f"Invalid memory specification: '{mem}'")

        if len(unit) == 0:
            warn("The amount of memory does not specify a unit; supposing MB")
            unit = "mb"
        else:
            if unit not in MEMORY_FACTORS:
                raise InvalidParameter(f"Unknown unit used for memory specification: '{unit}'")

        try:
            val = float(_val)
        except ValueError:
            raise InvalidParameter(f"The amount of memory must be a number (received '{mem}')")
        standard_mem = val*MEMORY_FACTORS[unit]
        return verify_memory_valid(standard_mem)
    else:
        raise InvalidParameter(f"Invalid type for the memory specification: {type(mem)}")


def standardize_xyz(xyz):
    """
        Converts variations of the XYZ format into a uniform format for this project.

        Format:
            "El1 X1 Y1 Z2\n El2 X2 Y2 Z2 [...]"

        Example:
            H 0.0 0.0 0.0
            H 1.0 0.0 0.0
    """

    standard_xyz = ""
    if isinstance(xyz, list):
        arr_xyz = xyz
    elif isinstance(xyz, str):
        arr_xyz = xyz.strip().split('\n')
    else:
        raise InvalidParameter(f"Cannot parse xyz from type {type(xyz)}")

    # Check if the xyz contains the two line header
    # Remove it if it does
    try:
        num_atoms = int(arr_xyz[0])
    except ValueError:
        pass
    else:
        if num_atoms == len(arr_xyz)-2:
            arr_xyz = arr_xyz[2:]
        else:
            raise InvalidXYZ(f"Invalid xyz header: {num_atoms} atoms specified, " +
                    "but actually contains {len(arr_xyz) - 2} atoms")

    for el in arr_xyz:
        line_data = []

        if not isinstance(el, str):
            raise InvalidXYZ(f"Could not parse xyz from array: contains element '{el}' " +
                    "and not string")

        if el.strip() == '':
            continue

        sel = el.strip().split()
        if not len(sel) == 4:
            raise InvalidXYZ(f"Invalid xyz: found line '{el}'")

        if sel[0] not in ATOMIC_NUMBER:
            if sel[0].isdigit():
                try:
                    el_Z = int(sel[0])
                except ValueError:
                    raise InvalidXYZ(f"Invalid atomic label: '{sel[0]}'")
                else:
                    if el_Z not in ATOMIC_SYMBOL:
                        raise InvalidXYZ(f"Invalid atomic number: '{el_Z}'")
                    el_symb = ATOMIC_SYMBOL[el_Z]
                    line_data.append(el_symb)
            else:
                if sel[0].lower() in LOWERCASE_ATOMIC_SYMBOLS:
                    line_data.append(LOWERCASE_ATOMIC_SYMBOLS[sel[0].lower()])
                else:
                    raise InvalidXYZ(f"Invalid atomic label: '{sel[0]}'")
        else:
            line_data.append(sel[0])

        for coord in sel[1:]:
            try:
                c = float(coord)
            except ValueError:
                raise InvalidXYZ(f"Invalid atomic coordinate: '{coord}'")
            else:
                line_data.append(c)

        standard_xyz += "{:<2} {:>12.8f} {:>12.8f} {:>12.8f}\n".format(*line_data)

    return standard_xyz

def parse_xyz_from_file(path):
    if not os.path.isfile(path):
        raise InvalidParameter(f"Input file not found: {path}")

    with open(path) as f:
        lines = f.readlines()

    return standardize_xyz(lines)

def get_theory_level(method):
    for level, methods in THEORY_LEVELS.items():
        if method.lower() in methods:
            return level
    return 'dft'

def get_abs_type(str_type):
    """
        Converts a string calculation type into the correct CalcType.
        Takes into account different equivalent ways to write the calculation types.
    """
    _str_type = str_type.lower().strip()
    if _str_type in STR_TYPES:
        return STR_TYPES[_str_type]
    else:
        raise InvalidParameter(f"Invalid calculation type: '{str_type}'")

def get_abs_software(software):
    _software = software.lower().strip()
    for s in SYN_SOFTWARE:
        if _software in SYN_SOFTWARE[s] or _software == s:
            return s
    raise InvalidParameter(f"Unknown software package: '{software}'")

def get_abs_method(method):
    _method = method.strip().lower()
    for m in SYN_METHODS:
        if _method in SYN_METHODS[m] or _method == m:
            return m
    raise InvalidParameter(f"Unknown method: '{method}'")

def get_abs_basis_set(basis_set):
    _bs = basis_set.strip().lower()
    for bs in SYN_BASIS_SETS:
        if _bs.lower() in SYN_BASIS_SETS[bs] or _bs.lower() == bs:
            return bs
    raise InvalidParameter(f"Unknown basis set: '{basis_set}'")

def get_abs_solvent(solvent):
    _solvent = solvent.strip().lower()
    if _solvent in ["", "vacuum", "vac"]:
        return ""
    for solv in SYN_SOLVENTS:
        if _solvent in SYN_SOLVENTS[solv] or _solvent == solv:
            return solv
    raise InvalidParameter(f"Unknown solvent: '{solvent}'")

def is_exchange_correlation_combination(method):
    for x in EXCHANGE_FUNCTIONALS:
        if method[:len(x)] == x:
            if method[len(x):] in CORRELATION_FUNCTIONALS:
                return EXCHANGE_FUNCTIONALS[method[:len(x)]] + \
                       CORRELATION_FUNCTIONALS[method[len(x):]]
    return False

def get_method(method, software):
    try:
        abs_method = get_abs_method(method)
    except InvalidParameter:
        if software == 'gaussian':
            # As far as I know, this kind of specification does not apply to ORCA
            if method.lower()[0] in ['u', 'r']:
                try:
                    abs_method = get_abs_method(method[1:])
                except InvalidParameter:
                    pass
                else:
                    if abs_method not in SOFTWARE_METHODS[software]:
                        warn(f"Unknown method '{method}'")
                        return method

                    return method[0].upper() + SOFTWARE_METHODS[software][abs_method]

            xc_check = is_exchange_correlation_combination(method.lower())
            if isinstance(xc_check, str):
                return xc_check
        warn(f"Unknown method '{method}'")
        return method
    else:
        if abs_method not in SOFTWARE_METHODS[software]:
            warn(f"Unknown method '{method}'")
            return method

        return SOFTWARE_METHODS[software][abs_method]

def get_basis_set(basis_set, software):
    try:
        abs_basis_set = get_abs_basis_set(basis_set)
    except InvalidParameter:
        warn(f"Unknown basis set '{basis_set}'")
        return basis_set

    return SOFTWARE_BASIS_SETS[software][abs_basis_set]

def get_solvent(solvent, software, solvation_model="smd"):
    try:
        abs_solvent = get_abs_solvent(solvent)
    except InvalidParameter:
        warn(f"Unknown solvent '{solvent}'")
        return solvent

    if software == "orca" and abs_solvent == "n-octanol":
        # Weird exception in ORCA
        if solvation_model == "smd":
            return "1-octanol"
        elif solvation_model == "cpcm":
            return "octanol"
        # Note that ch2cl2 is a valid keyword for SMD, although not listed in the manual

    return SOFTWARE_SOLVENTS[software][abs_solvent]

def has_dispersion_parameters(method, version="d3"):
    try:
        _method = get_abs_method(method)
    except InvalidParameter:
        warn("Unknown method, could not verify if dispersion parameters are available for it")
        return True # Don't print a second warning

    if _method in FUNCTIONALS_WITH_DISPERSION_PARAMETERS[version]:
        return True
    return False

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

def get_npxyz(str_xyz):
    lines = [i + '\n' for i in clean_xyz(str_xyz).split('\n')]

    xyz = []
    for line in lines:
        if line.strip() != '':
            a, x, y, z = line.split()
            xyz.append([a, np.array([float(x), float(y), float(z)])])

    return xyz

def get_coord(xyz, ids):
    if len(ids) == 2:
        return get_distance(xyz, *ids)
    elif len(ids) == 3:
        return get_angle(xyz, *ids)
    elif len(ids) == 4:
        return get_dihedral(xyz, *ids)
    else:
        raise InvalidParameter(f"Invalid number of atoms: {len(ids)}")



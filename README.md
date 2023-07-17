# Computational Chemistry Input Generator

[![codecov](https://codecov.io/gh/cyllab/ccinput/branch/main/graph/badge.svg?token=ox4smJs0vh)](https://codecov.io/gh/cyllab/ccinput)
[![PyPI version](https://badge.fury.io/py/ccinput.svg)](https://badge.fury.io/py/ccinput)
[![Downloads](https://pepy.tech/badge/ccinput)](https://pepy.tech/project/ccinput)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5907043.svg)](https://doi.org/10.5281/zenodo.5907043)

`ccinput` is an application to generate input files for computational chemistry software.

Example usage:
```
$ ccinput gaussian opt M062X -bs def2tzvp -f ethanol.xyz -n 8 --mem 32G
%chk=calc.chk
%nproc=8
%mem=32000MB
#p opt M062X/Def2TZVP

File created by ccinput

0 1
C   -1.31970000  -0.64380000   0.00000000
H   -0.96310000  -1.65260000   0.00000000
H   -0.96310000  -0.13940000  -0.87370000
H   -2.38970000  -0.64380000   0.00000000
C   -0.80640000   0.08220000   1.25740000
H   -1.16150000   1.09160000   1.25640000
H   -1.16470000  -0.42110000   2.13110000
O    0.62360000   0.07990000   1.25870000
H    0.94410000   0.53240000   2.04240000

$ ccinput orca sp DLPNO-CCSDT -bs cc-pVTZ -f *.xyz -o .inp -n 32 --mem 250GB --parse_name
*** No C auxiliary basis set specified, using cc-pVTZ ***
*** No C auxiliary basis set specified, using cc-pVTZ ***
*** No C auxiliary basis set specified, using cc-pVTZ ***
*** No C auxiliary basis set specified, using cc-pVTZ ***
Input file written to carbo_cation.inp
Input file written to CH4.inp
Input file written to ethanol.inp
Input file written to H2.inp

$ cat H2.inp
!SP DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C
*xyz 0 1
H   -1.43860000  -0.26820000   0.00000000
H   -2.03860000  -0.26820000   0.00000000
*
%MaxCore 7812
%pal
nprocs 32
end

$ cat carbo_cation.inp
!SP DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C
*xyz 1 1
C   -2.90770000   1.66310000   0.00000000
H   -2.55100000   2.16750000   0.87370000
H   -3.97770000   1.66310000   0.00000000
H   -2.55100000   2.16750000  -0.87370000
*
%MaxCore 7812
%pal
nprocs 32
end
```

`ccinput` can also be used as python library:

```
>>> from ccinput.wrapper import gen_input
>>> inp = gen_input(software="orca", type="ts", method="PBEh-3c", file="ethanol.xyz", nproc=16, solvent="ethanol", solvation_model="SMD")
*** No solvation radii specified; using default radii ***
>>> print(inp)
!OPTTS PBEh-3c
*xyz 0 1
C   -1.31970000  -0.64380000   0.00000000
H   -0.96310000  -1.65260000   0.00000000
H   -0.96310000  -0.13940000  -0.87370000
H   -2.38970000  -0.64380000   0.00000000
C   -0.80640000   0.08220000   1.25740000
H   -1.16150000   1.09160000   1.25640000
H   -1.16470000  -0.42110000   2.13110000
O    0.62360000   0.07990000   1.25870000
H    0.94410000   0.53240000   2.04240000
*
%pal
nprocs 16
end
%cpcm
smd true
SMDsolvent "ethanol"
end
>>>
```

## Installation
### From PyPI
```
pip install ccinput
```

### From Github
You can install the bleeding-edge version of `ccinput` from Github:
```
pip install git+https://github.com/cyllab/ccinput
```

## Usage
`ccinput` supports a wide range of options, including different solvation radii, density fitting and multiple basis sets. As of now, Gaussian 16, ORCA 5 and xtb are supported, and more packages will be added in the future.

### From the command line
Simply use the `ccinput` command with the desired parameters:
```

usage: ccinput [-h] [--basis_set BASIS_SET] [--solvent SOLVENT] [--solvation_model SOLVATION_MODEL] 
               [--solvation_radii SOLVATION_RADII] [--custom_solvation_radii CUSTOM_SOLVATION_RADII] 
               [--specifications SPECIFICATIONS] [--density_fitting DENSITY_FITTING] [--custom_basis_sets CUSTOM_BASIS_SETS] 
               [--xyz XYZ] [--file FILE [FILE ...]] [--output OUTPUT] [--constraints CONSTRAINTS] [--freeze ATOM [ATOM ...]] 
               [--scan ATOM [ATOM ...]] [--from FROM] [--to TO] [--nsteps NSTEPS] [--step STEP] [--nproc NPROC] [--mem MEM] 
               [--charge CHARGE] [--mult MULT] [--parse_name] [--trust_me] [--d3 | --d3bj] [--name NAME] 
               [--aux_name AUX_NAME] [--header HEADER] [--save SAVE] [--preset [PRESET]] [--driver {none,ORCA,pysis}] 
               [--version] [--fragments FRAGMENTS]
               [software] [type] [method]

```

More detailed information about each option can be obtained with the `ccinput -h` command.

### As Python library
The function `gen_input` returns input files as a single strings with the correct whitespace.

```
>>> from ccinput.wrapper import gen_input
>>> inp = gen_input(...)
```

The input can also be directly written to a file using `write_input`.
```
>>> from ccinput.wrapper import write_input
>>> write_input(filename, ...)
```

See the [documentation](https://ccinput.readthedocs.io/en/latest/usage.html) for all options.

## Contributors

- @RaphaelRobidas: xtb, Gaussian, ORCA and Pysisyphus support
- @zarkoivkovicc: NWChem support and Counterpoise calculation in Gaussian
- @lauragranmar: Basic Psi4 support
- @berquist: Basic QChem support

## Contributing

We welcome all contributions to the project. This includes:

- Support for new packages (even if rudimentary)
- Support for new features of supported packages
- Correction or improvement of static data (like DFT functionals or basis sets)

Planned tasks are listed in the [roadmap](https://github.com/orgs/cyllab/projects/1/views/1). This can be a good place to start when looking to contribute, although do not limit yourself to what is listed there. The contribution guidelines are detailed in CONTRIBUTING.md.

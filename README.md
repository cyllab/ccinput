# Computational Chemistry Input Generator

[![codecov](https://codecov.io/gh/cyllab/ccinput/branch/main/graph/badge.svg?token=ox4smJs0vh)](https://codecov.io/gh/cyllab/ccinput)
[![Downloads](https://pepy.tech/badge/ccinput)](https://pepy.tech/project/ccinput)

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

```

`ccinput` can also be used as python library:

```
>>> from ccinput.wrapper import gen_input
>>> inp = gen_input(software="orca", type="ts", method="PBEh-3c", in_file="ethanol.xyz", nproc=16, solvent="ethanol", solvation_model="SMD")
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
`ccinput` supports a wide range of options, including different solvation radii, density fitting and multiple basis sets. As of now, Gaussian 16 and ORCA 5 are supported, and more packages will be added in the future.
### From the command line
Simply use the `ccinput` command with the desired parameters:
```
usage: ccinput [-h] [--basis_set BASIS_SET] [--solvent SOLVENT] [--solvation_model SOLVATION_MODEL]
               [--solvation_radii SOLVATION_RADII] [--specifications SPECIFICATIONS]
               [--density_fitting DENSITY_FITTING] [--custom_basis_sets CUSTOM_BASIS_SETS] [--xyz XYZ]
               [--file FILE] [--output OUTPUT] [--constraints CONSTRAINTS] [--nproc NPROC] [--mem MEM]
               [--charge CHARGE] [--mult MULT] [--name NAME] [--header HEADER]
               software type method
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

The desired options are passed as additional named parameters:
```
gen_input(software=None, # "gaussian" or "orca"
	type=None, # Type of calculation ("sp", "opt", "freq", ...)
	method="", # Computational method ("HF", "AM1", "B3LYP", ...)
	basis_set="", # If required by the method
	solvent="", # If implicit solvation is desired
	solvation_model="", # SMD, PCM, CPCM
	solvation_radii="", # Blank or "Default" gives the default radii, other options are available ("SMD18", "Bondi", "UAKS", ...)
	specifications="", # Any other custom keywords
	density_fitting="", # If desired
	custom_basis_sets="", # Basis sets for specific atoms, uses the format "<Element>=<Basis set label>;..." (e.g. "I=Def2-TZVPD;")
	xyz="", # XYZ structure as a string (replaces 'in_file')
	in_file="", # Filename of XYZ structure (replaces 'xyz')
	constraints="", # Constrained optimisation constraints (pending documentation)
	nproc=0, # Number of CPU cores to use
	mem="", # Amount of memory to use (specify units, otherwise MB are assumed)
	charge=0, 
	multiplicity=1, 
	name="calc", # If used by the software package, name of the input/scratch files
	header="File created by ccinput", # If used by the software package, header in the input file
	)
```

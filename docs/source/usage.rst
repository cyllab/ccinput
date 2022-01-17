Usage
=====

Command Line Usage
------------------

.. argparse::
   :ref: ccinput.wrapper.get_parser
   :prog: ccinput

Python Library Usage
--------------------

When using ``ccinput`` as a Python library, the desired wrapper must be imported as such:

.. code-block:: python

        >>> from ccinput.wrapper import gen_input
        >>> inp = gen_input(...)

The generated input can be directly written to files:

.. code-block:: python

        >>> from ccinput.wrapper import write_input
        >>> inp = write_input(filename, ...)

The parameters to these two wrappers will be the same as in the command line usage, using their long name as named parameters. For example:

.. code-block:: python

        >>> from ccinput.wrapper import gen_input
        >>> inp = gen_input(software="orca", type="ts", method="PBEh-3c", in_file="ethanol.xyz", nproc=16, solvent="ethanol", solvation_model="SMD")

Parameter Details
------------------

Basis sets are not case sensitive and ``ccinput`` recognizes several synonyms for each of them. For example, "Def2-TZVP" can be written in any of these forms:

- def2-TZVP
- DEF2TZVP
- def2tzvp


Software
^^^^^^^^

Software package for which the input should be generated. This can influence the availability of certain options.

As of version 1.0, the following packages are supported for most features:

.. exec::
   from ccinput.constants import SYN_SOFTWARE
   from ccinput.documentation import format_dict_str
   print(format_dict_str(SYN_SOFTWARE, "Software id"))

Type
^^^^

Type of calculation (*e.g.* single-point energy, geometrical optimization, transition state optimization...)

.. exec::
   from ccinput.constants import INV_STR_TYPES
   from ccinput.documentation import format_dict_enum
   print(format_dict_enum(INV_STR_TYPES, "Calculation id"))

|

Method
^^^^^^

Calculation method to use.

.. note::

   Coupled cluster methods (CCSD...) and MP2 are not currently available
        
.. collapse:: All computational methods

        .. exec::
           from ccinput.constants import SYN_METHODS
           from ccinput.documentation import format_dict_str
           print(format_dict_str(SYN_METHODS, "Method"))

|

Basis set
^^^^^^^^^

Basis set used for the calculation. Required for most methods, except "-3c" methods (*e.g.* HF-3c, PBEh-3c) and semi-empirical methods (*e.g.* AM1, PM3, ...), including tight-binding methods (*e.g.* GFN2-xTB).

.. collapse:: All basis sets

        .. exec::
           from ccinput.constants import SYN_BASIS_SETS
           from ccinput.documentation import format_dict_str
           print(format_dict_str(SYN_BASIS_SETS, "Basis set"))

|

Custom basis sets
^^^^^^^^^^^^^^^^^

Specific basis sets for specific atoms. Uses the format "<Element>=<Basis set label>;..." (*e.g.* "I=Def2-TZVPD;"). To see all the valid basis set labels per element, visit the `Basis Set Exchange <https://www.basissetexchange.org/>`_.

Density fitting
^^^^^^^^^^^^^^^

Basis set to use for density fitting.

`Gaussian documentation about density fitting <https://gaussian.com/basissets/>`_

.. note::

   Only available for Gaussian 16 for the moment

Solvent
^^^^^^^

Solvent to model using implicit solvation.

.. collapse:: All solvents

        .. exec::
           from ccinput.constants import SYN_SOLVENTS
           from ccinput.documentation import format_dict_str
           print(format_dict_str(SYN_SOLVENTS, "Solvent"))

|

Solvation model
^^^^^^^^^^^^^^^

Model used for implicit solvation.

========== ========
Software   Models
========== ========
Gaussian   SMD

           PCM

           CPCM

ORCA       SMD

           CPCM
========== ========

Solvation radii
^^^^^^^^^^^^^^^

Set of element radii to use in the solvation model.

.. note:

   Only the default radii and SMD18 radii are implemented in ORCA; the other radii can only be used in Gaussian


======= =================
Model   Sets of radii
======= =================
SMD     Default

        SMD18 [SMD18]_

All     UFF (g16 default)

        UA0 

        UAHF

        UAKS

        Pauling

        Bondi

======= =================

.. [SMD18] E. Engelage, N. Schulz, F. Heinen, S. M. Huber, D. G. Truhlar, C. J. Cramer, *Chem. Eur. J.* **2018**, *24*, 15983-15987.


Specifications
^^^^^^^^^^^^^^

Custom keywords to add to the command of the input.

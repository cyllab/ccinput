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
        >>> inp = gen_input(software="orca", type="ts", method="PBEh-3c", file="ethanol.xyz", nproc=16, solvent="ethanol", solvation_model="SMD")

Parameter Details
------------------

The parameters not case sensitive and ``ccinput`` recognizes several synonyms for each of them. For example, "Def2-TZVP" can be written in any of these forms:

- def2-TZVP
- DEF2TZVP
- def2tzvp

Spaces and symbols other than ``+`` and ``*``  will be ignored. The otherwise different synonyms will be listed for each parameter. If you believe that a common or useful synonym is missing, feel free to send a pull request.

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
   from ccinput.documentation import format_calc_types
   print(format_calc_types())

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

Specific basis sets for specific atoms. Uses the format "<Element>=<Basis set label>;..." (*e.g.* "I=Def2-TZVPD;"). To see all the valid basis set labels per element, visit the `Basis Set Exchange <https://www.basissetexchange.org/>`_. Nonetheless, ``ccinput`` will also detect close synonyms of the requested basis set.

If applicable, the effective core potential (ECP) corresponding to the requested basis set will also be added to the input file.

Density fitting
^^^^^^^^^^^^^^^

Basis set to use for density fitting.

`Gaussian documentation about density fitting <https://gaussian.com/basissets/>`_

.. note::

   Only available for Gaussian 16 for the moment

Structure files
^^^^^^^^^^^^^^^

Structure file(s) to use in the input. Only XYZ files are currently supported.

Multiple files can be specified at once when using from the command line:

.. code-block:: console

        $ ccinput [...] -f struct1.xyz struct2.xyz

If no output pattern is specified, each input file will be printed to the console sequentially separated by a header. With an output pattern, the files will be created in the chosen directory with the given prefix and extension. When specifying only one file, the exact output path will be used.

.. code-block:: console

        $ ccinput [...] -f struct1.xyz struct2.xyz -o calc_dir/sp.inp
        Input file written to calc_dir/sp_struct1.inp
        Input file written to calc_dir/sp_struct2.inp

        $ ccinput [...] -f struct1.xyz struct2.xyz -o .com
        Input file written to struct1.com
        Input file written to struct2.com

        $ ccinput [...] -f struct1.xyz -o my_struct.com
        Input file written to my_struct.com


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

xtb        GBSA

           ALPB

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

Custom solvation radii
^^^^^^^^^^^^^^^^^^^^^^

Specific solvation radii to use for some elements.

The format to use is "<ELEMENT>=<RADIUS>;...", for example: "H=1.00;Li=1.70;".

Specifications
^^^^^^^^^^^^^^

Custom keywords to add to the command of the input. 

.. code-block:: console

        $ ccinput orca opt HF --xyz "Cl 0 0 0" -c -1 -bs Def2SVP
        !OPT HF Def2-SVP
        [...]

        $ ccinput orca opt HF --xyz "Cl 0 0 0" -c -1 -bs Def2SVP --specifications "TIGHTSCF"
        !OPT HF Def2-SVP tightscf
        [...]

With Gaussian 16, this can also add parameters to the calculation keyword:

.. code-block:: console

        $ ccinput g16 opt HF --xyz "Cl 0 0 0" -c -1 -bs Def2SVP
        [...]
        #p opt HF/Def2SVP
        [...]

        $ ccinput g16 opt HF --xyz "Cl 0 0 0" -c -1 -bs Def2SVP --specifications "opt(maxstep=5)"
        [...]
        #p opt(maxstep=5) HF/Def2SVP
        [...]

        $ ccinput g16 opt HF --xyz "Cl 0 0 0" -c -1 -bs Def2SVP --specifications "opt(maxstep=5) SCF(restart)"
        [...]
        #p opt(maxstep=5) HF/Def2SVP scf(restart)
        [...]

Note that the specifications are not checked for validity beyond simple syntax checks. This allows you to use all valid keywords of the software, but can also lead to invalid inputs.

Constraints
^^^^^^^^^^^

Constraints (freeze or scan) can be specified either as a string or as separate parameters.

The string representation allows to specified all the constraints at once in a relatively readable fashion:

.. code-block:: console

   # Freeze the bond between atoms 1 and 2 (starting at 1)
   $ ccinput [...] --constraints "freeze/1_2;"

   # Freeze the angle between atoms 1, 2 and 3
   $ ccinput [...] --constraints "freeze/1_2_3;"

   # Freeze the dihedral angle between atoms 1, 2, 3 and 4
   $ ccinput [...] --constraints "freeze/1_2_3_4;"

   # Different constraints are delimited by semi-colons
   $ ccinput [...] --constraints "freeze/1_2;freeze/3_4_5;"

   # Scans the bond 1-2 from 2 A to 1 A in 10 steps
   # Note that Gaussian does not allow starting values other than those of the current structure.
   # As such, the starting value is ignored with Gaussian.
   $ ccinput [...] --constraints "scan_2_1_10/1_2;"

   # Scans the angle 1-2-3 from 90 degrees to 0 degrees in 10 steps
   $ ccinput [...] --constraints "scan_90_0_10/1_2_3;"

   # Scans the angle 1-2-3 from its current value to 0 degrees in 10 steps
   $ ccinput [...] --constraints "scan_auto_0_10/1_2_3;"

   # Scans the dihedral angle 1-2-3-4 from 90 degrees to 0 degrees in 10 steps
   $ ccinput [...] --constraints "scan_90_0_10/1_2_3_4;"

   # Different types of constraint can be combined
   $ ccinput [...] --constraints "scan_90_0_10/1_2_3;freeze/4_5;"

The library usage uses an identical syntax:

.. code-block:: python

        >>> from ccinput.wrapper import gen_input
        >>> inp = gen_input([...], constraints="scan_90_0_10/1_2_3;freeze/4_5;")

With the string, scans always require the starting value, the final value as well as the number of steps. However, more convenient options are available using ``--scan``, ``--from``, ``--to``, ``--nsteps`` and ``--step``:

.. code-block:: console

   # Also scans the bond 1-2 from 2 A to 1 A in 10 steps
   $ ccinput [...] --scan 1 2 --from 2 --to 1 --nsteps 10

   # Scans the bond 1-2 from its current value to 1 A in 10 steps
   $ ccinput [...] --scan 1 2 --to 1 --nsteps 10

   # Scans the bond 1-2 from its current value to 1 A in steps of 0.1 A
   $ ccinput [...] --scan 1 2 --to 1 --step 0.1

   # Equivalent to the above
   $ ccinput [...] --scan 1 2 --to 1 --step -0.1

   # Scans the bonds 1-2 and 3-4 from their current values to 1 A in steps of 0.1 A
   $ ccinput [...] --scan 1 2 --to 1 --step 0.1 --scan 3 4 --to 1 --step 0.1

   # The order of the sets of parameters matters:
   # Scans the bond 1-2 from 1.9 A to 0.5 A by step of 0.15 and 
   # the bond 3-4 from 1.5 A to 1.0 A in step of 0.1 A
   $ ccinput [...] --scan 1 2 --from 1.9 --to 0.5 --step 0.15 --scan 3 4 --from 1.5 --to 1 --step 0.1

   # Scans the bond 1-2 from 1.5 A to 1.0 A in step of 0.1 A and 
   # the bond 3-4 from 1.9 A to 0.5 A by step of 0.15
   $ ccinput [...] --scan 1 2 --from 1.5 --to 1 --step 0.1 --scan 3 4 --from 1.9 --to 0.5 --step 0.15 

   # However, the exact ordering of each different parameter does not matter
   # Equivalent to the above
   $ ccinput [...] --scan 1 2 --to 1 --from 1.5 --step 0.1 --scan 3 4 --step 0.15 --to 0.5 --from 1.9 

   # Moreover, coordinates can be frozen in a similar fashion
   # Freezes the angle between atoms 1, 2 and 3
   $ ccinput [...] --freeze 1 2 3

   # Freezes the angle between atoms 1, 2 and 3 and 
   # scans the bond 1-2 from its current value to 1 A in steps of 0.1 A
   $ ccinput [...] --freeze 1 2 3 --scan 1 2 --to 1 --step 0.1

   # Equivalent to the above
   $ ccinput [...] --scan 1 2 --to 1 --step 0.1 --freeze 1 2 3 

   # Also equivalent to the above (although confusing)
   $ ccinput [...] --scan 1 2 --to 1 --freeze 1 2 3 --step 0.1


The library usage employs arrays for each parameter, with ``freeze`` and ``scan`` being arrays of arrays (for multiple constraints of multiple atoms each). The scanning parameters are prefixed by the letter "s" (for "scan") due to the name clash with the Python keyword ``from``.

.. code-block:: python

        >>> from ccinput.wrapper import gen_input
        >>> inp = gen_input([...], freeze=[[1, 2]])

        >>> inp = gen_input([...], scan=[[2, 3]], sfrom=[1.0], sto=[1.5], snsteps=[5])

        >>> inp = gen_input([...], scan=[[2, 3]], sfrom=[1.0], sto=[1.5], sstep=[0.1])

        >>> inp = gen_input([...], scan=[[2, 3], [1, 2, 3, 4]], sfrom=[1.0, 120], sto=[1.5, 160], sstep=[0.1, 5])

Presets
-------
Presets offer a convenient way to save sets of parameters and reuse them easily in the command line:

.. code-block:: console

        $ ccinput gaussian opt m062x -bs def2tzvp --specifications "opt(maxstep=5) 5d nosymm" -sm smd -sr smd18 -s methanol --save my_preset
        --- Saved preset 'my_preset'
        software                      gaussian                          
        type                          opt                           
        method                        m062x                         
        basis_set                     def2tzvp                      
        version                       <Y.X.Z>        
        solvent                       methanol                      
        solvation_model               smd                           
        solvation_radii               smd18                         
        specifications                opt(maxstep=5) 5d nosymm

        $ ccinput --preset my_preset --xyz "Cl 0 0 0" -c -1
        %chk=calc.chk
        %nproc=1
        %mem=1000MB
        #p opt(maxstep=5) M062X/Def2TZVP 5d nosymm SCRF(SMD, Solvent=methanol, Read)

        File created by ccinput

        -1 1
        Cl   0.00000000   0.00000000   0.00000000

        modifysph

        Br 2.60
        I 2.74

All parameters can be saved in the preset, except the calculation name and the XYZ structure. To create a preset, simply enter all the desired parameters exactly like when creating an input file and append ``--preset <preset_name>`` to the command. This will not generate any input file and will instead same the parameters as JSON in an accessible user directory (generally ``~/.local/share/ccinput`` or ``C:\\Users\\<username>\\AppData\\Local\\CYLlab\\ccinput``).

Parameters in the preset file can be overwritten when creating the input file: 

.. code-block:: console

        $ ccinput --preset my_preset --xyz "Cl 0 0 0" -c -1
        %chk=calc.chk
        %nproc=1
        %mem=1000MB
        #p opt(maxstep=5) M062X/Def2TZVP 5d nosymm SCRF(SMD, Solvent=methanol, Read)

        File created by ccinput

        -1 1
        Cl   0.00000000   0.00000000   0.00000000

        modifysph

        Br 2.60
        I 2.74


        $ ccinput --preset my_preset --xyz "Cl 0 0 0" -c -1 -s vacuum
        %chk=calc.chk
        %nproc=1
        %mem=1000MB
        #p opt(maxstep=5) M062X/Def2TZVP 5d nosymm

        File created by ccinput

        -1 1
        Cl   0.00000000   0.00000000   0.00000000

Specifications are the only parameters that are not overwritten, but combined:

.. code-block:: console

        $ ccinput --preset my_preset --xyz "Cl 0 0 0" -c -1 -s vacuum --specifications "Int(Ultrafinegrid)"

        %chk=calc.chk
        %nproc=1
        %mem=1000MB
        #p opt(maxstep=5) M062X/Def2TZVP 5d nosymm int(ultrafinegrid)

        File created by ccinput

        -1 1
        Cl   0.00000000   0.00000000   0.00000000

Presets can be permanently modified by specifying new parameters and saving them in the same preset. Unspecified options will not be modified:

.. code-block:: console

        $ ccinput -s vacuum --save my_preset
        --- Saved preset 'my_preset'
        software                      gaussian                      
        type                          opt                           
        method                        m062x                         
        basis_set                     def2tzvp                      
        version                       1.2.2+7.g8e241e2.dirty        
        solvent                       vacuum                        
        solvation_model               smd                           
        solvation_radii               smd18                         
        specifications                opt(maxstep=5) 5d nosymm

        $ ccinput --preset my_preset --xyz "Cl 0 0 0" -c -1
        %chk=calc.chk
        %nproc=1
        %mem=1000MB
        #p opt(maxstep=5) M062X/Def2TZVP 5d nosymm

        File created by ccinput

        -1 1
        Cl   0.00000000   0.00000000   0.00000000


Note that the parameters are not validated on preset creation.

Overview
========

``ccinput`` is an application to generate input files for computational chemistry software.

Example usage:

.. code-block:: console

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

``ccinput`` can also be used as python library:

.. code-block:: python

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


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

Supported features
------------------

Below is a short summary of the different features that can be requested for the supported packages. Important features (*e.g.* a whole different calculation type) that are not mentioned can be assumed to not be supported. Minor features (*e.g.* additional printout) that are not mentioned may or may not be supported; refer to the :doc:`detailed usage<usage>` for more information.

=============================== ============ ===========
Calculation Type                Gaussian 16    ORCA 5
=============================== ============ ===========
Essential calculations [1]_        yes         yes
Common calculations [2]_           yes         yes
Minimum energy path                n.a.        yes
NMR prediction                     yes         yes
TD-DFT                             yes         no
Molecular orbital visualisation    no          yes
=============================== ============ ===========


=============================== ============ ===========
Level of theory                 Gaussian 16    ORCA 5
=============================== ============ ===========
Tight-binding methods              no          yes [3]_
Semi-empirical methods             yes         yes
Hartree-fock                       yes         yes
Density Functional Theory          yes         yes
Grimme's "3c" methods              n.a.        yes
Møller-Plesset                     no          no
Coupled Cluster                    n.a.        no
=============================== ============ ===========


=============================== ============ ===========
Feature                         Gaussian 16    ORCA 5
=============================== ============ ===========
Implicit solvation                 yes         yes
Choice of solvation radii          yes         yes
Custom basis sets                  yes         yes
Density fitting                    yes         no
Custom additional keywords         yes         yes
=============================== ============ ===========


.. [1] Single-point energy calculation, geometrical optimisation, frequency calculation

.. [2] Transition state optimisation, constrained optimisation

.. [3] Requires the `xtb package <https://github.com/grimme-lab/xtb>`__; the supported methods are GFN2-xTB, GFN1-xTB, GFN0-xTB and GFN-FF

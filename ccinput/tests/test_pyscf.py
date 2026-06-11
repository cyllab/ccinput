from ccinput.tests.testing_utilities import InputTests
from ccinput.exceptions import InvalidParameter, ImpossibleCalculation


class PySCFTests(InputTests):
    def test_sp(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "PySCF",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = '''
        import os, json, pyscf
        from pyscf.lib import num_threads
			from pyscf import dft

        os.environ["OMP_NUM_THREADS"] = "8"
        os.environ["MKL_NUM_THREADS"] = "8"
        os.environ["OPENBLAS_NUM_THREADS"] = "8"
        num_threads(8)

        mol = pyscf.gto.M(
        atom="""Cl 0.0 0.0 0.0""",
        basis="def2svp",
        charge=-1,
        spin=0,
        verbose=5,
        )

        mol.max_memory = 1000
        mf = dft.rks.RKS(mol, xc="M062X")
        mf.kernel()
        print(json.dumps({
            "SCF": mf.e_tot,
        }))'''

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_df(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "PySCF",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "def2svp/JK",
        }

        inp = self.generate_calculation(**params)

        REF = '''
        import os, json, pyscf
        from pyscf.lib import num_threads
        from pyscf import dft

        os.environ["OMP_NUM_THREADS"] = "8"
        os.environ["MKL_NUM_THREADS"] = "8"
        os.environ["OPENBLAS_NUM_THREADS"] = "8"
        num_threads(8)

        mol = pyscf.gto.M(
        atom="""Cl 0.0 0.0 0.0""",
        basis="def2svp",
        charge=-1,
        spin=0,
        verbose=5,
        )

        mol.max_memory = 1000
        mf = dft.rks.RKS(mol, xc="M062X")
        mf = mf.density_fit(auxbasis="def2svp")
        mf.kernel()
        print(json.dumps({
            "SCF": mf.e_tot,
        }))'''

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_grid(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "PySCF",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "specifications": "grid(5)",
        }

        inp = self.generate_calculation(**params)

        REF = '''
        import os, json, pyscf
        from pyscf.lib import num_threads
        from pyscf import dft

        os.environ["OMP_NUM_THREADS"] = "8"
        os.environ["MKL_NUM_THREADS"] = "8"
        os.environ["OPENBLAS_NUM_THREADS"] = "8"
        num_threads(8)

        mol = pyscf.gto.M(
        atom="""Cl 0.0 0.0 0.0""",
        basis="def2svp",
        charge=-1,
        spin=0,
        verbose=5,
        grid=5,
        )

        mol.max_memory = 1000
        mf = dft.rks.RKS(mol, xc="M062X")
        mf.kernel()
        print(json.dumps({
            "SCF": mf.e_tot,
        }))'''

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_functional(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "PySCF",
            "method": "PBE0",
            "basis_set": "Def2-SVP",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = '''
        import os, json, pyscf
        from pyscf.lib import num_threads
        from pyscf import dft

        os.environ["OMP_NUM_THREADS"] = "8"
        os.environ["MKL_NUM_THREADS"] = "8"
        os.environ["OPENBLAS_NUM_THREADS"] = "8"
        num_threads(8)

        mol = pyscf.gto.M(
        atom="""Cl 0.0 0.0 0.0""",
        basis="def2svp",
        charge=-1,
        spin=0,
        verbose=5,
        )

        mol.max_memory = 1000
        mf = dft.rks.RKS(mol, xc="PBE0")
        mf.kernel()
        print(json.dumps({
            "SCF": mf.e_tot,
        }))'''

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_nproc(self):
        params = {
            "nproc": 4,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "PySCF",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
        }

        inp = self.generate_calculation(**params)

        REF = '''
        import os, json, pyscf
        from pyscf.lib import num_threads
        from pyscf import dft

        os.environ["OMP_NUM_THREADS"] = "4"
        os.environ["MKL_NUM_THREADS"] = "4"
        os.environ["OPENBLAS_NUM_THREADS"] = "4"
        num_threads(4)

        mol = pyscf.gto.M(
        atom="""Cl 0.0 0.0 0.0""",
        basis="def2svp",
        charge=-1,
        spin=0,
        verbose=5,
        )

        mol.max_memory = 1000
        mf = dft.rks.RKS(mol, xc="M062X")
        mf.kernel()
        print(json.dumps({
            "SCF": mf.e_tot,
        }))'''

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_mem(self):
        params = {
            "nproc": 8,
            "type": "Single-Point Energy",
            "file": "Cl.xyz",
            "software": "PySCF",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "charge": "-1",
            "mem": "10G",
        }

        inp = self.generate_calculation(**params)

        REF = '''
        import os, json, pyscf
        from pyscf.lib import num_threads
        from pyscf import dft

        os.environ["OMP_NUM_THREADS"] = "8"
        os.environ["MKL_NUM_THREADS"] = "8"
        os.environ["OPENBLAS_NUM_THREADS"] = "8"
        num_threads(8)

        mol = pyscf.gto.M(
        atom="""Cl 0.0 0.0 0.0""",
        basis="def2svp",
        charge=-1,
        spin=0,
        verbose=5,
        )

        mol.max_memory = 10000
        mf = dft.rks.RKS(mol, xc="M062X")
        mf.kernel()
        print(json.dumps({
            "SCF": mf.e_tot,
        }))'''

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_freq(self):
        params = {
            "nproc": 8,
            "type": "Frequency Calculation",
            "file": "benzene.xyz",
            "software": "PySCF",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
        }

        inp = self.generate_calculation(**params)

        REF = '''
        import os, json, pyscf
        from pyscf.lib import num_threads
        from pyscf.hessian.thermo import thermo, harmonic_analysis
        from pyscf import dft

        os.environ["OMP_NUM_THREADS"] = "8"
        os.environ["MKL_NUM_THREADS"] = "8"
        os.environ["OPENBLAS_NUM_THREADS"] = "8"
        num_threads(8)

        mol = pyscf.gto.M(
        atom="""C         -0.99790       -1.58780       -0.02080
                C          0.39730       -1.58780       -0.02080
                C          1.09480       -0.38010       -0.02080
                C          0.39720        0.82840       -0.02200
                C         -0.99760        0.82830       -0.02250
                C         -1.69520       -0.37990       -0.02150
                H         -1.54760       -2.54020       -0.02030
                H          0.94680       -2.54040       -0.01950
                H          2.19450       -0.38000       -0.02020
                H          0.94740        1.78060       -0.02200
                H         -1.54780        1.78060       -0.02340
                H         -2.79480       -0.37970       -0.02160""",
        basis="def2svp",
        charge=0,
        spin=0,
        verbose=5,
        )

        mol.max_memory = 1000
        mf = dft.rks.RKS(mol, xc="M062X")
        mf.kernel()
        ohess = mf.Hessian()
        hess = ohess.kernel()
        res_vib = harmonic_analysis(mf.mol, hess, imaginary_freq=False)
        res_thermo = thermo(mf, res_vib["freq_au"])
        print(json.dumps({
            "SCF": mf.e_tot,
            "H": res_thermo["H_tot"][0],
            "G": res_thermo["G_tot"][0],
            "hessian": hess.tolist(),
            "freqs": res_vib["freq_wavenumber"].tolist(),
            "modes": res_vib["norm_mode"].tolist(),
        }))'''

        self.assertTrue(self.is_equivalent(REF, inp.input_file))

    def test_ts(self):
        params = {
            "nproc": 8,
            "type": "TS Optimization",
            "file": "small_ts.xyz",
            "software": "PySCF",
            "method": "M06-2X",
            "basis_set": "Def2-SVP",
            "specifications": "verbose=4",
        }

        inp = self.generate_calculation(**params)

        REF = '''
            import os, json, pyscf
            from pyscf.lib import num_threads
            from pyscf import dft

            os.environ["OMP_NUM_THREADS"] = "8"
            os.environ["MKL_NUM_THREADS"] = "8"
            os.environ["OPENBLAS_NUM_THREADS"] = "8"
            num_threads(8)

            mol = pyscf.gto.M(
            atom="""C   -1.30199679200336      0.23770149604323     -0.21903830900536
            C   -0.65452549088797      1.37677004974294      0.23784083906211
            H   -0.89567737188508      2.31097104509050     -0.24568023160592
            H   -0.51131365003345      1.45158760940331      1.30393146403304
            O   -0.98734706501576     -0.93896491644235      0.24464647161953
            C   0.57938721313073     -1.21542290559438      0.60370912581272
            C   1.30373709071931     -0.40352018279841     -0.31679957603999
            C   1.21326659501348      0.97578532623407     -0.11796726468775
            H   1.54268192054322      1.35508505207923      0.83816378823672
            H   1.50364201896930      1.60996020761999     -0.94605411105798
            H   1.33180817109859     -0.76458664777933     -1.33953946285788
            H   0.61967054666194     -2.29023171829519      0.47066760098394
            H   0.64159071923251     -0.92315930061295      1.64911403264797
            H   -1.77094330174335      0.23668457701018     -1.20684978624530""",
            basis="def2svp",
            charge=0,
            spin=0,
            verbose=4,
            )

            mol.max_memory = 1000
            mf = dft.rks.RKS(mol, xc="M062X")
            mf.kernel()
            params = {'transition': True, 'trust': 0.02, 'tmax': 0.06, "hessian": "file:hessian.txt"}
            mol_ts = mf.Gradients().optimizer(solver='geomeTRIC').kernel(params)
            print(json.dumps({
                "xyz": mol_ts.tostring(),
            }))
            '''
        self.assertTrue(self.is_equivalent(REF, inp.input_file))


# grid level

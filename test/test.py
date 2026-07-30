import unittest as ut
import subprocess as sb
import h5py
import numpy as np
from pathlib import Path
import shutil
import tempfile

class TestFundamentalExecution(ut.TestCase):
    # testing the execution of rixs-pathway-serial for diamond example
    def test_exec_pathway_diamond(self):
        # test whether the code runs without error
        proc_=sb.Popen(r'../../../../bin/rixs-pathway-serial', 
                cwd=r'./data/diamond/pathway')
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        data_=h5py.File('./data/diamond/pathway/data.h5','r')
        ref_=h5py.File('./data/diamond/pathway/data_ref.h5','r')

        #test the eigenvalues
        self.assertEqual(data_['evals'].shape[0],
                ref_['evals'].shape[0])
        np.testing.assert_array_equal(data_['evals'], ref_['evals'])        
        #test the t(1) matrix
        np.testing.assert_array_equal(data_['t(1)'].shape,
                ref_['t(1)'].shape)
        np.testing.assert_array_equal(data_['t(1)'], ref_['t(1)'])        
        #test the t(2) matrix
        np.testing.assert_array_equal(data_['t(2)'].shape,
                ref_['t(2)'].shape)
        np.testing.assert_array_equal(data_['t(2)'], ref_['t(2)'])        

    # testing the execution of rixs-pathway-serial and
    # rixs-oscstr-serial for diamond example
    def test_exec_oscstr_diamond(self):
        # test whether both pathway and oscstr run without error
        proc1_=sb.Popen(r'../../../../bin/rixs-pathway-serial', 
                cwd=r'./data/diamond/pathway')
        out1, err1=proc1_.communicate()
        self.assertEqual(proc1_.returncode,0)
        proc2_=sb.Popen(r'../../../../bin/rixs-oscstr-serial', 
                cwd=r'./data/diamond/pathway')
        out, err=proc2_.communicate()
        self.assertEqual(proc2_.returncode,0)

        #open rixs.h5 and rixs_ref.h5
        rixs_=h5py.File('./data/diamond/pathway/rixs.h5')
        ref_=h5py.File('./data/diamond/pathway/rixs_ref.h5')

        #test the eigenvalues
        np.testing.assert_array_equal(rixs_['cevals'].shape,
                ref_['cevals'].shape)
        np.testing.assert_array_equal(rixs_['cevals'], ref_['cevals'])
        np.testing.assert_array_equal(rixs_['vevals'].shape,
                ref_['vevals'].shape)
        np.testing.assert_array_equal(rixs_['vevals'], ref_['vevals'])        
        #test the oscstr matrices
        self.assertEqual(len(rixs_['oscstr'].keys()),
                len(ref_['oscstr'].keys()))
        for entry in rixs_['oscstr'].keys():
            np.testing.assert_allclose(rixs_['oscstr'][entry],
                    ref_['oscstr'][entry], rtol=5e-12, atol=1e-14)

    # testing the execution of rixs-pathway-serial for lif example
    def test_exec_pathway_lif(self):
        # test whether the code runs without error
        proc_=sb.Popen(r'../../../bin/rixs-pathway-serial', 
                cwd=r'./data/lif')
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        data_=h5py.File('./data/lif/data.h5','r')
        ref_=h5py.File('./data/lif/data_ref.h5','r')

        #test the eigenvalues
        self.assertEqual(data_['evals'].shape[0],
                ref_['evals'].shape[0])
        np.testing.assert_array_equal(data_['evals'], ref_['evals'])        
        #test the t(1) matrix
        np.testing.assert_array_equal(data_['t(1)'].shape,
                ref_['t(1)'].shape)
        np.testing.assert_array_equal(data_['t(1)'], ref_['t(1)'])        
        #test the t(2) matrix
        np.testing.assert_array_equal(data_['t(2)'].shape,
                ref_['t(2)'].shape)
        np.testing.assert_array_equal(data_['t(2)'], ref_['t(2)'])        
    
    # testing the execution of rixs-pathway-serial and
    # rixs-oscstr-serial for lif example
    def test_exec_oscstr_lif(self):
        # test whether both pathway and oscstr run without error
        proc1_=sb.Popen(r'../../../bin/rixs-pathway-serial', 
                cwd=r'./data/lif')
        out1, err1=proc1_.communicate()
        self.assertEqual(proc1_.returncode,0)
        proc2_=sb.Popen(r'../../../bin/rixs-oscstr-serial', 
                cwd=r'./data/lif')
        out, err=proc2_.communicate()
        self.assertEqual(proc2_.returncode,0)

        #open rixs.h5 and rixs_ref.h5
        rixs_=h5py.File('./data/lif/rixs.h5')
        ref_=h5py.File('./data/lif/rixs_ref.h5')

        #test the eigenvalues
        np.testing.assert_array_equal(rixs_['cevals'].shape,
                ref_['cevals'].shape)
        np.testing.assert_array_equal(rixs_['cevals'], ref_['cevals'])     
        np.testing.assert_array_equal(rixs_['vevals'].shape,
                ref_['vevals'].shape)
        np.testing.assert_array_equal(rixs_['vevals'], ref_['vevals'])    
        #test the oscstr matrices
        self.assertEqual(len(rixs_['oscstr'].keys()),
                len(ref_['oscstr'].keys()))
        for entry in rixs_['oscstr'].keys():
            np.testing.assert_allclose(rixs_['oscstr'][entry],
                    ref_['oscstr'][entry], rtol=5e-12, atol=1e-14)

class TestCoherenceExecution(ut.TestCase):
    # testing the execution of rixs-coherence-serial for diamond example
    def test_exec_coherence_diamond(self):
        # test whether the code runs without error
        proc_=sb.Popen(r'../../../../bin/rixs-coherence-serial', 
                cwd=r'./data/diamond/coherence/')
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        rixs_=h5py.File('./data/diamond/coherence/rixs.h5','r')
        ref_=h5py.File('./data/diamond/coherence/rixs_ref.h5','r')

        #test the eigenvalues
        self.assertEqual(rixs_['cevals'].shape[0],
                ref_['cevals'].shape[0])
        np.testing.assert_array_equal(rixs_['cevals'], ref_['cevals'])        
        self.assertEqual(rixs_['vevals'].shape[0],
                ref_['vevals'].shape[0])
        np.testing.assert_array_equal(rixs_['vevals'], ref_['vevals']) 
        
        #test the oscstr matrices
        self.assertEqual(len(rixs_['oscstr'].keys()),
                len(ref_['oscstr'].keys()))
        for entry in rixs_['oscstr'].keys():
            np.testing.assert_allclose(rixs_['oscstr'][entry]['coherent'],
                    ref_['oscstr'][entry]['coherent'], rtol=5e-12, atol=1e-14)
            np.testing.assert_allclose(rixs_['oscstr'][entry]['incoherent'],
                    ref_['oscstr'][entry]['incoherent'], rtol=5e-12, atol=1e-14)


class TestNonEquilibriumPathways(ut.TestCase):
    def test_constant_occupation_factors_scale_pathways(self):
        """Both BSE amplitudes are weighted in their own transition spaces."""
        repository = Path(__file__).resolve().parents[1]
        source = repository / 'test' / 'data' / 'diamond' / 'pathway'
        executable = repository / 'bin' / 'rixs-pathway-serial'
        core_factor = 0.5
        optical_factor = 0.25

        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / 'pathway'
            shutil.copytree(source, workdir)

            with h5py.File(workdir / 'core_output.h5', 'r') as core:
                ncore = core['eigvec-singlet-TDA-BAR-full/0001/parameters/smap'].shape[0]
            with h5py.File(workdir / 'optical_output.h5', 'r') as optical:
                noptical = optical['eigvec-singlet-TDA-BAR-full/0001/parameters/smap'].shape[0]

            for filename, size, factor in (
                    ('occupations_core.h5', ncore, core_factor),
                    ('occupations_optical.h5', noptical, optical_factor)):
                with h5py.File(workdir / filename, 'w') as occupations:
                    group = occupations.require_group('IQMT_000001/transitions')
                    group.create_dataset('occupation_factors',
                                         data=np.full(size, factor, dtype=np.float64))

            config_path = workdir / 'input.cfg'
            config_text = config_path.read_text(encoding='utf-8')
            config_path.write_text(
                config_text.replace('nblocks=1', 'nblocks=2')
                + '\nnon_equilibrium=true\n',
                encoding='utf-8')

            proc = sb.run([executable], cwd=workdir, check=False)
            self.assertEqual(proc.returncode, 0)

            with h5py.File(workdir / 'data.h5', 'r') as data, \
                    h5py.File(workdir / 'data_ref.h5', 'r') as reference:
                np.testing.assert_allclose(
                    data['t(1)'][:], core_factor * reference['t(1)'][:],
                    rtol=2e-14, atol=1e-14)
                np.testing.assert_allclose(
                    data['t(2)'][:],
                    core_factor * optical_factor * reference['t(2)'][:],
                    rtol=2e-14, atol=1e-14)

    def test_transition_resolved_factors_match_weighted_eigenvectors(self):
        """Arbitrary factors act component-wise, before pathway contractions."""
        repository = Path(__file__).resolve().parents[1]
        source = repository / 'test' / 'data' / 'diamond' / 'pathway'
        executable = repository / 'bin' / 'rixs-pathway-serial'
        eigenvector_group = 'eigvec-singlet-TDA-BAR-full/0001/rvec'

        with tempfile.TemporaryDirectory() as tmp:
            nonequilibrium = Path(tmp) / 'nonequilibrium'
            weighted_reference = Path(tmp) / 'weighted_reference'
            shutil.copytree(source, nonequilibrium)
            shutil.copytree(source, weighted_reference)

            factors = {}
            for kind in ('core', 'optical'):
                output_name = f'{kind}_output.h5'
                with h5py.File(nonequilibrium / output_name, 'r') as output:
                    ntransitions = output[
                        'eigvec-singlet-TDA-BAR-full/0001/parameters/smap'
                    ].shape[0]
                factors[kind] = np.linspace(
                    0.2, 1.0, ntransitions, dtype=np.float64)
                with h5py.File(
                        nonequilibrium / f'occupations_{kind}.h5', 'w') as occupations:
                    group = occupations.require_group('IQMT_000001/transitions')
                    group.create_dataset('occupation_factors', data=factors[kind])

                with h5py.File(weighted_reference / output_name, 'r+') as output:
                    for eigenvector in output[eigenvector_group].values():
                        eigenvector[...] = eigenvector[...] * factors[kind][:, None]

            config_path = nonequilibrium / 'input.cfg'
            config_path.write_text(
                config_path.read_text(encoding='utf-8')
                + '\nnon_equilibrium=true\n',
                encoding='utf-8')

            result_nonequilibrium = sb.run(
                [executable], cwd=nonequilibrium, check=False)
            result_reference = sb.run(
                [executable], cwd=weighted_reference, check=False)
            self.assertEqual(result_nonequilibrium.returncode, 0)
            self.assertEqual(result_reference.returncode, 0)

            with h5py.File(nonequilibrium / 'data.h5', 'r') as data, \
                    h5py.File(weighted_reference / 'data.h5', 'r') as reference:
                np.testing.assert_allclose(
                    data['t(1)'][:], reference['t(1)'][:],
                    rtol=2e-14, atol=1e-14)
                np.testing.assert_allclose(
                    data['t(2)'][:], reference['t(2)'][:],
                    rtol=2e-14, atol=1e-14)

    def test_coherence_uses_both_occupation_factor_sets(self):
        repository = Path(__file__).resolve().parents[1]
        source = repository / 'test' / 'data' / 'diamond' / 'coherence'
        executable = repository / 'bin' / 'rixs-coherence-serial'
        core_factor = 0.5
        optical_factor = 0.25
        oscillator_factor = core_factor**2 * optical_factor

        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / 'coherence'
            shutil.copytree(source, workdir)

            for kind, factor in (
                    ('core', core_factor), ('optical', optical_factor)):
                with h5py.File(workdir / f'{kind}_output.h5', 'r') as output:
                    ntransitions = output[
                        'eigvec-singlet-TDA-BAR-full/0001/parameters/smap'
                    ].shape[0]
                with h5py.File(
                        workdir / f'occupations_{kind}.h5', 'w') as occupations:
                    group = occupations.require_group('IQMT_000001/transitions')
                    group.create_dataset(
                        'occupation_factors',
                        data=np.full(ntransitions, factor, dtype=np.float64))

            config_path = workdir / 'input.cfg'
            config_path.write_text(
                config_path.read_text(encoding='utf-8')
                + '\nnon_equilibrium=true\n',
                encoding='utf-8')

            proc = sb.run([executable], cwd=workdir, check=False)
            self.assertEqual(proc.returncode, 0)

            with h5py.File(workdir / 'rixs.h5', 'r') as result, \
                    h5py.File(workdir / 'rixs_ref.h5', 'r') as reference:
                for frequency in result['oscstr']:
                    for contribution in ('coherent', 'incoherent'):
                        np.testing.assert_allclose(
                            result['oscstr'][frequency][contribution][:],
                            oscillator_factor
                            * reference['oscstr'][frequency][contribution][:],
                            rtol=5e-12, atol=1e-14)


if __name__ == '__main__':
    ut.main()

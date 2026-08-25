import os
import unittest as ut
import subprocess as sb
import h5py
import numpy as np
from pathlib import Path
import shutil
import tempfile

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
# Overridden by test/CMakeLists.txt to point at the build tree; falls back
# to the installed BRIXS/bin/ for a manual `python test.py` run.
BIN_DIR = os.environ.get('BRIXS_BIN_DIR',
                          os.path.normpath(os.path.join(TEST_DIR, '..', 'bin')))


def binary(name):
    return os.path.join(BIN_DIR, name)


def data_dir(*parts):
    return os.path.join(TEST_DIR, 'data', *parts)


def assert_close(actual, desired):
    # Reference data was generated with ifort+MKL; other BLAS/compiler
    # combinations reorder floating-point sums, giving noise even
    # on a correct build.
    np.testing.assert_allclose(actual, desired, rtol=1e-12, atol=1e-14)


class TestFundamentalExecution(ut.TestCase):
    # testing the execution of rixs_pathway for diamond example
    def test_exec_pathway_diamond(self):
        # test whether the code runs without error
        proc_=sb.Popen(binary('rixs_pathway'),
                cwd=data_dir('diamond', 'pathway'))
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        data_=h5py.File(data_dir('diamond', 'pathway', 'data.h5'),'r')
        ref_=h5py.File(data_dir('diamond', 'pathway', 'data_ref.h5'),'r')

        #test the eigenvalues
        self.assertEqual(data_['evals'].shape[0],
                ref_['evals'].shape[0])
        np.testing.assert_array_equal(data_['evals'], ref_['evals'])
        #test the t(1) matrix
        np.testing.assert_array_equal(data_['t(1)'].shape,
                ref_['t(1)'].shape)
        assert_close(data_['t(1)'], ref_['t(1)'])
        #test the t(2) matrix
        np.testing.assert_array_equal(data_['t(2)'].shape,
                ref_['t(2)'].shape)
        assert_close(data_['t(2)'], ref_['t(2)'])

    # testing the execution of rixs_pathway and
    # rixs_oscstr for diamond example
    def test_exec_oscstr_diamond(self):
        # test whether both pathway and oscstr run without error
        proc1_=sb.Popen(binary('rixs_pathway'),
                cwd=data_dir('diamond', 'pathway'))
        out1, err1=proc1_.communicate()
        self.assertEqual(proc1_.returncode,0)
        proc2_=sb.Popen(binary('rixs_oscstr'),
                cwd=data_dir('diamond', 'pathway'))
        out, err=proc2_.communicate()
        self.assertEqual(proc2_.returncode,0)

        #open rixs.h5 and rixs_ref.h5
        rixs_=h5py.File(data_dir('diamond', 'pathway', 'rixs.h5'))
        ref_=h5py.File(data_dir('diamond', 'pathway', 'rixs_ref.h5'))

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
            assert_close(rixs_['oscstr'][entry], ref_['oscstr'][entry])

    # testing the execution of rixs_pathway for lif example
    def test_exec_pathway_lif(self):
        # test whether the code runs without error
        proc_=sb.Popen(binary('rixs_pathway'),
                cwd=data_dir('lif'))
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        data_=h5py.File(data_dir('lif', 'data.h5'),'r')
        ref_=h5py.File(data_dir('lif', 'data_ref.h5'),'r')

        #test the eigenvalues
        self.assertEqual(data_['evals'].shape[0],
                ref_['evals'].shape[0])
        np.testing.assert_array_equal(data_['evals'], ref_['evals'])
        #test the t(1) matrix
        np.testing.assert_array_equal(data_['t(1)'].shape,
                ref_['t(1)'].shape)
        assert_close(data_['t(1)'], ref_['t(1)'])
        #test the t(2) matrix
        np.testing.assert_array_equal(data_['t(2)'].shape,
                ref_['t(2)'].shape)
        assert_close(data_['t(2)'], ref_['t(2)'])

    # testing the execution of rixs_pathway and
    # rixs_oscstr for lif example
    def test_exec_oscstr_lif(self):
        # test whether both pathway and oscstr run without error
        proc1_=sb.Popen(binary('rixs_pathway'),
                cwd=data_dir('lif'))
        out1, err1=proc1_.communicate()
        self.assertEqual(proc1_.returncode,0)
        proc2_=sb.Popen(binary('rixs_oscstr'),
                cwd=data_dir('lif'))
        out, err=proc2_.communicate()
        self.assertEqual(proc2_.returncode,0)

        #open rixs.h5 and rixs_ref.h5
        rixs_=h5py.File(data_dir('lif', 'rixs.h5'))
        ref_=h5py.File(data_dir('lif', 'rixs_ref.h5'))

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
            assert_close(rixs_['oscstr'][entry], ref_['oscstr'][entry])

class TestCoherenceExecution(ut.TestCase):
    # testing the execution of rixs_coherence for diamond example
    def test_exec_coherence_diamond(self):
        # test whether the code runs without error
        proc_=sb.Popen(binary('rixs_coherence'),
                cwd=data_dir('diamond', 'coherence'))
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        rixs_=h5py.File(data_dir('diamond', 'coherence', 'rixs.h5'),'r')
        ref_=h5py.File(data_dir('diamond', 'coherence', 'rixs_ref.h5'),'r')

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
            assert_close(rixs_['oscstr'][entry]['coherent'],
                    ref_['oscstr'][entry]['coherent'])
            assert_close(rixs_['oscstr'][entry]['incoherent'],
                    ref_['oscstr'][entry]['incoherent'])


class TestNonEquilibriumPathways(ut.TestCase):
    def test_constant_occupation_factors_scale_pathways(self):
        """Both BSE amplitudes are weighted in their own transition spaces."""
        source = Path(data_dir('diamond', 'pathway'))
        executable = Path(binary('rixs-pathway-serial'))
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
                output_name = filename.replace('occupations_', '').replace('.h5', '_output.h5')
                with h5py.File(workdir / output_name, 'r') as output:
                    smap = output[
                        'eigvec-singlet-TDA-BAR-full/0001/parameters/smap'
                    ][:]
                with h5py.File(workdir / filename, 'w') as occupations:
                    group = occupations.require_group('IQMT_000001/transitions')
                    group.create_dataset('occupation_factors',
                                         data=np.full(size, factor, dtype=np.float64))
                    group.create_dataset('smap', data=smap)

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
                assert_close(data['t(1)'][:], core_factor * reference['t(1)'][:])
                assert_close(
                    data['t(2)'][:],
                    core_factor * optical_factor * reference['t(2)'][:])

    def test_transition_resolved_factors_match_weighted_eigenvectors(self):
        """Arbitrary factors act component-wise, before pathway contractions."""
        source = Path(data_dir('diamond', 'pathway'))
        executable = Path(binary('rixs-pathway-serial'))
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
                    with h5py.File(nonequilibrium / output_name, 'r') as output:
                        group.create_dataset(
                            'smap',
                            data=output[
                                'eigvec-singlet-TDA-BAR-full/0001/parameters/smap'
                            ][:])

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
                assert_close(data['t(1)'][:], reference['t(1)'][:])
                assert_close(data['t(2)'][:], reference['t(2)'][:])

    def test_coherence_uses_both_occupation_factor_sets(self):
        source = Path(data_dir('diamond', 'coherence'))
        executable = Path(binary('rixs-coherence-serial'))
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
                    with h5py.File(workdir / f'{kind}_output.h5', 'r') as output:
                        group.create_dataset(
                            'smap',
                            data=output[
                                'eigvec-singlet-TDA-BAR-full/0001/parameters/smap'
                            ][:])

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
                        assert_close(
                            result['oscstr'][frequency][contribution][:],
                            oscillator_factor
                            * reference['oscstr'][frequency][contribution][:])

    def test_non_equilibrium_requires_matching_transition_map(self):
        """Missing or reordered transition maps are rejected before calculation."""
        source = Path(data_dir('diamond', 'pathway'))
        executable = Path(binary('rixs-pathway-serial'))

        for mode in ('missing', 'reordered'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                workdir = Path(tmp) / 'pathway'
                shutil.copytree(source, workdir)

                for kind in ('core', 'optical'):
                    with h5py.File(workdir / f'{kind}_output.h5', 'r') as output:
                        smap = output[
                            'eigvec-singlet-TDA-BAR-full/0001/parameters/smap'
                        ][:]
                    with h5py.File(
                            workdir / f'occupations_{kind}.h5', 'w') as occupations:
                        group = occupations.require_group('IQMT_000001/transitions')
                        group.create_dataset(
                            'occupation_factors',
                            data=np.ones(smap.shape[0], dtype=np.float64))
                        if mode == 'reordered' or kind == 'optical':
                            if mode == 'reordered' and kind == 'core':
                                smap = smap.copy()
                                smap[[0, 1]] = smap[[1, 0]]
                            group.create_dataset('smap', data=smap)

                config_path = workdir / 'input.cfg'
                config_path.write_text(
                    config_path.read_text(encoding='utf-8')
                    + '\nnon_equilibrium=true\n',
                    encoding='utf-8')

                proc = sb.run(
                    [executable], cwd=workdir, check=False,
                    capture_output=True, text=True)
                self.assertNotEqual(proc.returncode, 0)
                output = proc.stdout + proc.stderr
                if mode == 'missing':
                    self.assertIn('required dataset', output)
                else:
                    self.assertIn('transition-map mismatch', output)


if __name__ == '__main__':
    ut.main()

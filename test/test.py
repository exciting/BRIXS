import os
import unittest as ut
import subprocess as sb
import h5py
import numpy as np

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
    np.testing.assert_allclose(actual, desired, rtol=1e-6, atol=1e-8)


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
if __name__ == '__main__':
    ut.main()

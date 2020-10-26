import unittest as ut
import subprocess as sb
import h5py
import numpy as np

class TestFundamentalExecution(ut.TestCase):
    # testing the execution of rixs-pathway-serial for diamond example
    def test_exec_pathway_diamond(self):
        # test whether the code runs without error
        proc_=sb.Popen(r'../../../bin/rixs-pathway-serial', 
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
        proc1_=sb.Popen(r'../../../bin/rixs-pathway-serial', 
                cwd=r'./data/diamond/pathway')
        out1, err1=proc1_.communicate()
        self.assertEqual(proc1_.returncode,0)
        proc2_=sb.Popen(r'../../../bin/rixs-oscstr-serial', 
                cwd=r'./data/diamond/pathway')
        out, err=proc2_.communicate()
        self.assertEqual(proc2_.returncode,0)

        #open rixs.h5 and rixs_ref.h5
        rixs_=h5py.File('./data/diamond/pathway/rixs.h5')
        ref_=h5py.File('./data/diamond/pathway/rixs_ref.h5')

        #test the eigenvalues
        np.testing.assert_array_equal(rixs_['evals'].shape,
                ref_['evals'].shape)
        np.testing.assert_array_equal(rixs_['evals'], ref_['evals'])        
        #test the oscstr matrices
        self.assertEqual(len(rixs_['oscstr'].keys()),
                len(ref_['oscstr'].keys()))
        for entry in rixs_['oscstr'].keys():
            np.testing.assert_array_equal(rixs_['oscstr'][entry],
                    ref_['oscstr'][entry])

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
        np.testing.assert_array_equal(rixs_['evals'].shape,
                ref_['evals'].shape)
        np.testing.assert_array_equal(rixs_['evals'], ref_['evals'])        
        #test the oscstr matrices
        self.assertEqual(len(rixs_['oscstr'].keys()),
                len(ref_['oscstr'].keys()))
        for entry in rixs_['oscstr'].keys():
            np.testing.assert_array_equal(rixs_['oscstr'][entry],
                    ref_['oscstr'][entry])

class TestCoherenceExecution(ut.TestCase):
    # testing the execution of rixs-coherence-serial for diamond example
    def test_exec_coherence_diamond(self):
        # test whether the code runs without error
        proc_=sb.Popen(r'../../../bin/rixs-coherence-serial', 
                cwd=r'./data/diamond/coherence/')
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        rixs_=h5py.File('./data/diamond/coherence/rixs.h5','r')
        ref_=h5py.File('./data/diamond/rixs_ref.h5','r')

        #test the eigenvalues
        self.assertEqual(data_['evals'].shape[0],
                ref_['evals'].shape[0])
        np.testing.assert_array_equal(data_['evals'], ref_['evals'])        
        
        #test the oscstr matrices
        self.assertEqual(len(rixs_['oscstr'].keys()),
                len(ref_['oscstr'].keys()))
        for entry in rixs_['oscstr'].keys():
            np.testing.assert_array_equal(rixs_['oscstr'][entry]['coherent'],
                    ref_['oscstr'][entry]['coherent'])
            np.testing.assert_array_equal(rixs_['oscstr'][entry]['incoherent'],
                    ref_['oscstr'][entry]['incoherent'])
if __name__ == '__main__':
    ut.main()

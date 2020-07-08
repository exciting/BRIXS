import unittest as ut
import subprocess as sb
import h5py
import numpy as np

class TestFundamentalExecution(ut.TestCase):
    # testing the execution of rixs-pathway-serial for diamond example
    def test_exec_pathwayi_diamond(self):
        proc_=sb.Popen(r'../../../bin/rixs-pathway-serial', 
                cwd=r'./data/Diamond')
        out, err=proc_.communicate()
        self.assertEqual(proc_.returncode,0)

        #open data.h5 and data_ref.h5
        data_=h5py.File('./data/Diamond/data.h5','r')
        ref_=h5py.File('./data/Diamond/data_ref.h5','r')

        #test the eigenvalues
        self.assertEqual(data_['evals'].shape[0],
                ref_['evals'].shape[0])
        np.testing.assert_array_equal(data_['evals'], ref_['evals'])        
        #test the t(1) matrix
        np.testing.assert_array_equal(data_['t(1)'].shape,
                ref_['t(1)'].shape)
        np.testing.assert_array_equal(data_['t(1)'], ref_['t(1)'])        
        #test the t(2) matrix
        np.testing.assert_array_equal(data_['t(1)'].shape,
                ref_['t(1)'].shape)
        np.testing.assert_array_equal(data_['t(1)'], ref_['t(1)'])        

    # testing the execution of rixs-pathway-serial and
    # rixs-oscstr-serial for diamond example
    def test_exec_oscstr(self):
        proc1_=sb.Popen(r'../../../bin/rixs-pathway-serial', 
                cwd=r'./data/Diamond')
        out1, err1=proc1_.communicate()
        self.assertEqual(proc1_.returncode,0)
        proc2_=sb.Popen(r'../../../bin/rixs-oscstr-serial', 
                cwd=r'./data/Diamond')
        out, err=proc2_.communicate()
        self.assertEqual(proc2_.returncode,0)
if __name__ == '__main__':
    ut.main()

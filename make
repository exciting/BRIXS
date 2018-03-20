#... with full debug options
#h5fc -g -traceback -check all -mkl -warn unused  mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_rixs.f90 rixs.f90 -o rixs
h5fc -g -traceback -check all -check bounds -mkl -warn unused  mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_blocks.f90 mod_rixs.f90 rixs_b.f90 -o rixs_b
#h5fc -g -traceback -check all -check bounds -mkl -warn unused  mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_blocks.f90 mod_rixs.f90 debug_b.f90 -o debug_b
#h5fc -g -traceback -check all -check bounds -mkl -warn unused  mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_blocks.f90 mod_rixs.f90 scaling_b.f90 -o scaling_b
#... production options
#h5fc  -mkl mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_rixs.f90 rixs.f90 -o rixs
#h5fc  -mkl mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_blocks.f90 mod_rixs.f90 rixs_b.f90 -o rixs_b
#h5fc -mkl  mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_blocks.f90 mod_rixs.f90 scaling_b.f90 -o scaling_b

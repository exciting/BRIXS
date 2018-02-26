#... with full debug options
#h5fc -g -traceback -check all -mkl -warn unused  mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_rixs.f90 rixs.f90 -o rixs
#... production options
h5fc  -mkl mod_hdf5.f90 mod_matmul.f90 mod_io.f90 mod_rixs.f90 rixs.f90 -o rixs

PHDF5=/users/stud/vorwerk/bin/phdf5-1.10.2/bin/h5pfc
HDF5=h5fc
debug:
	$(HDF5) -g -traceback -check all -check bounds -mkl -warn unused mod_config.f90 mod_hdf5.F90 mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o rixs_pathway
	$(HDF5) -g -traceback -check all -check bounds -mkl -warn unused mod_config.f90 mod_hdf5.F90 mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o rixs_osctr 

dune:
	$(PHDF5) -mkl -DMPI m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o rixs_pathway
	$(PHDF5) -mkl -DMPI m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o rixs_oscstr



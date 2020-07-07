#serial HDF5 (gfortran)
HDF5=/home1/bin/hdf5-1.10.1/bin/h5fc -L/usr/lib -llapack -L/usr/lib -lblas
#serial HDF5 (ifort+MKL)
#HDF5=/home1/bin/hdf5-1.10.1/bin/h5fc -mkl
#parallel HDF5 on DUNE
PHDF5_D=/users/stud/vorwerk/bin/phdf5-1.10.2/bin/h5pfc
#parallel HDF5 on HLRN
PHDF5=h5pfc -Wl,-rpath,/sw/dataformats/szip/intel.18/2.1.1/skl/lib/

serial:
	$(HDF5) m_config.f90 modmpi.F90 mod_phdf5.F90  mod_hdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o rixs_pathway
	$(HDF5)  m_config.f90 modmpi.F90 mod_phdf5.F90  mod_hdf5.F90 mod_matmul.F90 \
		mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o rixs_oscstr

debug:
	$(HDF5) -g -traceback -check all m_config.f90 modmpi.F90 mod_phdf5.F90 \
	 	mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o  \
		rixs_pathway
	$(HDF5) -g -traceback -check all  m_config.f90 modmpi.F90 mod_phdf5.F90 \
	 	mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o   \
	 	rixs_oscstr

dune:
	$(PHDF5_D) -mkl -DMPI m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o rixs_pathway
	$(PHDF5_D) -mkl -DMPI m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o rixs_oscstr

hlrn:
	$(PHDF5) -mkl -DMPI m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o rixs_pathway
	$(PHDF5) -mkl -DMPI m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o rixs_oscstr

clean:
	rm *.o *.mod rixs_pathway rixs_oscstr	

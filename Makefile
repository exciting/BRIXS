#serial HDF5
HDF5=h5fc
#parallel HDF5 on DUNE
PHDF5_D=/users/stud/vorwerk/bin/phdf5-1.10.2/bin/h5pfc
#parallel HDF5 on HLRN
PHDF5=h5pfc -Wl,-rpath,/sw/dataformats/szip/intel.18/2.1.1/skl/lib/

serial:
	$(HDF5) -mkl  m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o rixs_pathway
	$(HDF5) -mkl  m_config.f90 modmpi.F90 mod_phdf5.F90 mod_matmul.F90 \
		mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o rixs_oscstr

debug:
	$(HDF5) -g -traceback -check all -mkl  m_config.f90 modmpi.F90 mod_phdf5.F90 \
	 	mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o  \
		rixs_pathway
	$(HDF5) -g -traceback -check all -mkl  m_config.f90 modmpi.F90 mod_phdf5.F90 \
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



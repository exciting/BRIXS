HDF5_ROOT := /opt/software/pc2/EB-SW/software/HDF5/1.14.6-gompi-2025b

HDF5 := $(HDF5_ROOT)/bin/h5fc \
         -I$(HDF5_ROOT)/include \
         -fallow-argument-mismatch -w

PHDF5 := $(HDF5_ROOT)/bin/h5pfc \
          -I$(HDF5_ROOT)/include \
          -fallow-argument-mismatch -w

LIBS := -lflexiblas

serial:
	$(HDF5) m_config.f90 mod_mpi.F90 mod_phdf5.F90  mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 rixs_pathway.F90 -o rixs_pathway $(LIBS)
	$(HDF5)  m_config.f90 mod_mpi.F90 mod_phdf5.F90  mod_matmul.F90 \
		mod_io.F90 mod_blocks.F90 rixs_oscstr.F90 -o rixs_oscstr $(LIBS)
	$(HDF5)  m_config.f90 mod_mpi.F90 mod_phdf5.F90  mod_matmul.F90 \
		mod_io.F90 mod_blocks.F90 mod_blocks_k.F90 rixs_coherence.F90 -o \
	 	rixs_coherence $(LIBS)
	cp rixs_pathway ../bin/rixs-pathway-serial
	cp rixs_oscstr ../bin/rixs-oscstr-serial
	cp rixs_coherence ../bin/rixs-coherence-serial
	rm rixs_pathway rixs_oscstr	rixs_coherence

debug:
	$(HDF5) -g -traceback -check all m_config.f90 mod_mpi.F90 mod_phdf5.F90 \
	 	mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o  \
		rixs_pathway
	$(HDF5) -g -traceback -check all  m_config.f90 mod_mpi.F90 mod_phdf5.F90 \
	 	mod_matmul.F90 mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o   \
	 	rixs_oscstr
	cp rixs_pathway ../bin/rixs-pathway-debug
	cp rixs_oscstr ../bin/rixs-oscstr-debug
	rm rixs_pathway rixs_oscstr	

parallel:
	$(PHDF5) -mkl -DMPI m_config.f90 mod_mpi.F90 mod_phdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_pathway.F90 -o rixs_pathway $(LIBS)
	$(PHDF5) -mkl -DMPI m_config.f90 mod_mpi.F90 mod_phdf5.F90 mod_matmul.F90 \
	 	mod_io.F90 mod_blocks.F90 mod_rixs.F90 rixs_oscstr.F90 -o rixs_oscstr $(LIBS)
	cp rixs_pathway ../bin/rixs-pathway-mpi
	cp rixs_oscstr ../bin/rixs-oscstr-mpi
	rm rixs_pathway rixs_oscstr	

clean:
	rm *.o *.mod 



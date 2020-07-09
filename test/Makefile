# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

test:
	$(CONDA_ACTIVATE) pybrixs_test
	python test.py
clean:
	rm ./data/Diamond/rixs.h5 ./data/Diamond/data.h5
	rm ./data/LiF/rixs.h5 ./data/LiF/data.h5

# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
.PHONY: test
test:
	$(CONDA_ACTIVATE) pybrixs_test
	python test.py
clean:
	rm ./data/diamond/pathway/rixs.h5 ./data/diamond/pathway/data.h5
	rm ./data/lif/rixs.h5 ./data/lif/data.h5

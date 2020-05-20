# this must be run from the project root directory

dummy_file_message = "This is a dummy file used by make to track\
 if the files in this folder are up to date!!!\nManual modifications\
 of the files in this folder will be undetected by make!!!"

dummy_file_name = .dummy_makefile.txt
	
.PHONY: clean all test

all: rath_data

# Change this to point to your env bin or alias
PYTHON3 = python # I used python 3.7.5
JULIA = julia    # I used julia 1.1.0 (2019-01-21)

# RathData
# This produce all the files, later used in the scripts, that contains
# the experimental data (see README.md).
# The julia pakage will use it too
rath_data_dummy_file = data/processed/rath2017___data/${dummy_file_name}
$(rath_data_dummy_file): scripts/RathData/1_convert_rath_data.jl
	$(JULIA) --project scripts/RathData/1_convert_rath_data.jl
	-@touch ${rath_data_dummy_file}
	-@echo ${dummy_file_message} > ${rath_data_dummy_file}

rath_data: ${rath_data_dummy_file}
clear_rath_data:
	rm -fr data/processed/rath2017___data

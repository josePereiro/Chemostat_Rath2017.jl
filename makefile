# this must be run from the project root directory

dummy_file_message = "This is a dummy file used by make to track\
 if the files in this folder are up to date!!!\nManual modifications\
 of the files in this folder will be undetected by make!!!"

dummy_file_name = .dummy_makefile.txt
	
.PHONY: clean all 
all: \
	rath_data \
	MODEL1105100000

clear: \
	rath_data_clear \
	MODEL1105100000_clear

# Change this to point to your env bin or alias
PYTHON3 = python3 # I used python 3.7.5
JULIA = julia  --project  # I used julia 1.1.0 (2019-01-21) (do not forget the --project)

##############################################
# RathData
##############################################
# This produce all the files, later used in the scripts, that contains
# the experimental data (see README.md).
# The julia pakage will use it too
rath_data_dummy_file = data/processed/rath2017___data/${dummy_file_name}
$(rath_data_dummy_file): scripts/RathData/1_convert_rath_data.jl
	$(JULIA) scripts/RathData/1_convert_rath_data.jl
	-@touch ${rath_data_dummy_file}
	-@echo ${dummy_file_message} > ${rath_data_dummy_file}

rath_data: ${rath_data_dummy_file}
rath_data_clear:
	rm -fr data/processed/rath2017___data


###############################################
# MODEL1105100000
###############################################
# model mat file
data/processed/MODEL1105100000/MODEL1105100000_url.mat: scripts/MODEL1105100000/0_make_mat_file.py
	$(PYTHON3) $^

MODEL1105100000_mat_file : data/processed/MODEL1105100000/MODEL1105100000_url.mat



# mets_map
data/processed/MODEL1105100000/mets_map.csv: scripts/MODEL1105100000/1_mets_map.jl
	$(JULIA) $^

MODEL1105100000_mets_map: data/processed/MODEL1105100000/mets_map.csv



# niklas biomass
data/processed/MODEL1105100000/niklas_biomass.csv: scripts/MODEL1105100000/1_niklas_biomass.jl
	$(JULIA) $^

MODEL1105100000_niklas_biomass: data/processed/MODEL1105100000/niklas_biomass.csv



# preparing base models
data/processed/MODEL1105100000/fva_preprocessed_base_model.jls: \
	scripts/MODEL1105100000/2_prepare_base_models.jl
	$(JULIA) $^

MODEL1105100000_preparing_base_model: data/processed/MODEL1105100000/fva_preprocessed_base_model.jls



MODEL1105100000: \
	MODEL1105100000_mat_file \
	MODEL1105100000_mets_map \
	MODEL1105100000_niklas_biomass \
	MODEL1105100000_preparing_base_model

MODEL1105100000_clear:
	rm -fr data/processed/MODEL1105100000

# -*- coding: utf-8 -*-
import os
import errno
import cobra

# TODO find a way to do this in COBRA.jl

# +
# Here we take MODEL1105100000_url model and make a mat_file from it.
# -

# ### Meta

model_name = 'MODEL1105100000'


# ### Loading SMBL

raw_gems_dir = f"./data/raw/{model_name}"
if not os.path.exists(raw_gems_dir):
    raise FileNotFoundError(
        errno.ENOENT, os.strerror(errno.ENOENT), raw_gems_dir)

print(f"Model {model_name} from Shlomi et al., Genome-Scale Metabolic Modeling Elucidates the Role of Proliferative Adaptation in Causing the Warburg Effect.")
# TODO search model download link
model_file = os.path.join(raw_gems_dir, f"{model_name}_url.xml")
print('Loading ', model_file, "...")
model = cobra.io.read_sbml_model(model_file)
print('Model loaded', (len(model.metabolites), len(model.reactions)))


# ### Saving mat

processed_gems_dir = f"./data/processed/{model_name}"
if not os.path.exists(processed_gems_dir):
    os.makedirs(processed_gems_dir)
    print(f"created {processed_gems_dir}")

print("Saving mat file ...")
mat_file = os.path.join(processed_gems_dir, f"{model_name}_url.mat")
cobra.io.save_matlab_model(model, mat_file, varname = 'model')
print(f"created {mat_file}")

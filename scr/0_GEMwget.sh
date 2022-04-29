#!/bin/bash
# Get metabolic models from HMA
# johaGL 2022
cd $HOME
mkdir GEMdwld
cd GEMdwld

echo "Getting metabolic models from HMA (GEMwget.sh)"

# brain cancer cell line 
wget https://metabolicatlas.org/api/v2/repository/models/homo_sapiens/cell-line_specific_models/brain/U-251MG.xml-91080a939b86d903928cb7e2c321c2ff.zip

# healtly brain
wget https://metabolicatlas.org/api/v2/repository/models/homo_sapiens/tissue-specific_models/brain/brain.xml-0cc137f6eff1e44d85b0813c28a482f7.zip

unzip brain*.zip
unzip U-251MG*.zip

mv 'U-251 MG.xml' U-251MG.xml

echo "ok, finished GEMwget.sh"

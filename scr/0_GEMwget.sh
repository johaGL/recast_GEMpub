#!/bin/bash
# Get metabolic models from HMA
# in local:
#   ./0_GEMwget.sh $HOME/datapub
# in bb8:
#    ./0_GEMwget.sh /scratch/CBIB/jgalvis/datapub
# johaGL 2022

echo $1


if [ -d $1 ];then
	if [ !  $1 != "" ]; then echo "as argument empty, output to HOME"; fi
	
	cd $1	
	if [ ! -d GEM/ ]; then mkdir GEM/; fi
	cd GEM/

	echo "Getting metabolic models from HMA (0_GEMwget.sh)"

	# brain cancer cell line 
	wget https://metabolicatlas.org/api/v2/repository/models/homo_sapiens/cell-line_specific_models/brain/U-251MG.xml-91080a939b86d903928cb7e2c321c2ff.zip

	# healthy brain
	wget https://metabolicatlas.org/api/v2/repository/models/homo_sapiens/tissue-specific_models/brain/brain.xml-0cc137f6eff1e44d85b0813c28a482f7.zip

	unzip brain*.zip
	unzip U-251MG*.zip

	mv 'U-251 MG.xml' U-251MG.xml

	echo "ok, finished 0_GEMwget.sh"
	
else
	echo "ERROR"
	printf "Example : in local:\n \
	   ./0_GEMwget.sh $HOME/datapub \n\
	  in bb8: \n\
	   ./0_GEMwget.sh /scratch/CBIB/jgalvis/datapub \n"
fi




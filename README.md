###Requirements###
The code is developed with [Python](https://www.python.org/) 3.5 
and [RDKit](https://github.com/rdkit/rdkit/tree/Release_2015_03_1).

###Use cases###
This section contains information about basic use cases. For more information,
or alternative usage, please refer to the scripts help. The help is accessible
under the -h option of every script.

#####Extract fragments and compute descriptors####
cd molecule_preparation

python extract_fragments.py -i {directory with SDF files} -o {JSON output}

python compute_descriptors.py -f -i {JSON output} -o {output CSV}
-p {[PaDEL](http://www.yapcwsoft.com/dd/padeldescriptor/PaDEL-Descriptor.zip) directory}

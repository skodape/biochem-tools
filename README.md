The code is developed with [Python](https://www.python.org/) 3.5 
and [RDKit](https://github.com/rdkit/rdkit/tree/Release_2015_03_1).

### Use cases
This section contains information about basic use cases. For more information,
or alternative usage, please refer to the scripts documentation. 

#### Extract fragments and compute descriptor
The [PaDEL](http://www.yapcwsoft.com/dd/padeldescriptor/PaDEL-Descriptor.zip) 
software is required.

```
cd molecule_preparation
python extract_fragments.py -i {directory with SDF files} -o {JSON output} [-f {fragments to generate}]
python compute_descriptors.py -f -i {JSON output} -o {output CSV} -p {PaDEL directory}
```

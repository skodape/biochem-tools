#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Load molecules from a SDF file.

"""

from rdkit import Chem

__author__ = 'Petr Škoda'
__credits__ = ['Petr Škoda']
__license__ = 'X11'
__version__ = '1.0.0'
__maintainer__ = 'Petr Škoda'
__status__ = 'Development'


def _main():
    sdf_file_path = input('Path to a SDF file:')
    for molecule in Chem.SDMolSupplier(sdf_file_path):
        # RDKit may fail to read some molecules - skip such molecules.
        if molecule is None:
            print('Invalid molecule detected!')
            continue
        print(molecule.GetProp('_Name'), ':', molecule.GetNumAtoms())


if __name__ == '__main__':
    _main()

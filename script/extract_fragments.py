#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extract and save fragments from given molecules.

"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals
import os
import argparse
import logging
import json
import rdkit
import rdkit.Chem

__author__ = 'Petr Škoda'
__credits__ = ['Petr Škoda']
__license__ = 'X11'
__version__ = '1.0.0'
__maintainer__ = 'Petr Škoda'
__status__ = 'Development'


def extract_fragments(molecule):
    """Return fragments in given molecule.

    :param molecule:
    :return:
    """
    pattern = rdkit.Chem.MolFromSmarts('*~*~*')
    output = []
    for atoms in molecule.GetSubstructMatches(pattern):
        smiles = rdkit.Chem.MolFragmentToSmiles(molecule,
                                                atomsToUse=list(atoms),
                                                kekuleSmiles=True)
        output.append({
            'smiles': smiles,
            'atoms': atoms
        })
    return output


def _read_settings():
    """Get and return application settings.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='extract molecular fragments')
    parser.add_argument('-i', type=str, dest='input',
                        help='input SDF file', required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='output JSON file', required=True)
    args = vars(parser.parse_args())
    #
    output = {
        'source': args['input'],
        'target': args['output']
    }
    return output


def _main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    setting = _read_settings()
    output = []
    for molecule in rdkit.Chem.SDMolSupplier(setting['source']):
        output.append({
            'name': molecule.GetProp('_Name'),
            'fragments': extract_fragments(molecule)
        })
    #
    if not os.path.exists(os.path.dirname(setting['target'])):
        os.makedirs(os.path.dirname(setting['target']))

    with open(setting['target'], 'w') as output_stream:
        json.dump(output, output_stream)


if __name__ == '__main__':
    _main()

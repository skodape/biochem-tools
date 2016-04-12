#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extract fragments for molecules in given SDF files and save then into JSON.

Usage:
    python extract_fragments.py
        -i {dir with sdf files}
        -o {path to output}
"""

import os
import argparse
import logging
import json
import rdkit
import rdkit.Chem

__author__ = 'Petr Škoda'
__license__ = 'X11'


def extract_fragments(molecule):
    """Return fragments for given molecule.

    :param molecule:
    :return:
    """
    pattern = rdkit.Chem.MolFromSmarts('*~*~*')
    output = []
    for atoms in molecule.GetSubstructMatches(pattern):
        smiles = rdkit.Chem.MolFragmentToSmiles(
            molecule, atomsToUse=list(atoms), kekuleSmiles=True)
        output.append({
            'smiles': smiles
        })
    return output


def read_configuration():
    """Get and return application settings.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='extract molecular fragments')
    parser.add_argument('-i', type=str, dest='input',
                        help='path to the input',
                        required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='output JSON file', required=True)
    parser.add_argument('-r', dest='recursive',
                        help='recursive scan', action='store_true',
                        required=False)

    return vars(parser.parse_args())


def load_sdf(path):
    """Generate molecules from SDF file.

    :param path:
    """
    logging.info('Loading: %s' % path)
    for molecule in rdkit.Chem.SDMolSupplier(path):
        if molecule is None:
            logging.error('Invalid molecule detected.')
            continue
        yield ({
            'name': molecule.GetProp('_Name'),
            'smiles': rdkit.Chem.MolToSmiles(molecule),
            'fragments': extract_fragments(molecule)
        })


def recursive_scan_for_sdf(path, recursive):
    """Perform recursive scan with call callback on all SDF files.

    :param path:
    :param recursive
    :return:
    """
    result = []
    for file_name in os.listdir(path):
        file_path = path + '/' + file_name
        if os.path.isdir(file_path):
            if recursive:
                result.extend(recursive_scan_for_sdf(file_path, recursive))
        elif os.path.isfile(file_path) and file_name.lower().endswith('.sdf'):
            result.append(file_path)
    return result


def write_molecule_json(output_stream, molecules, holder):
    """Write given molecule as a JSON into stream.

    Optionaly put separator before the record based on 'holder'.
    :param output_stream:
    :param molecules:
    :param holder:
    :return:
    """
    for molecule in molecules:
        if holder['first']:
            holder['first'] = False
        else:
            output_stream.write(',')
        json.dump(molecule, output_stream)


def create_parent_directory(path):
    """Create directory if it does not exists.

    :param path:
    :return:
    """
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    configuration = read_configuration()
    if os.path.isdir(configuration['input']):
        input_files = recursive_scan_for_sdf(configuration['input'],
                                             configuration['recursive'])
    else:
        input_files = [configuration['input']]
    #
    create_parent_directory(configuration['output'])
    with open(configuration['output'], 'w') as output_stream:
        holder = {'first': True}
        #
        output_stream.write('[')
        for path in input_files:
            write_molecule_json(
                output_stream,
                load_sdf(path),
                holder)
        output_stream.write(']')


if __name__ == '__main__':
    main()

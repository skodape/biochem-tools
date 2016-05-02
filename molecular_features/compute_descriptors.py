#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute PaDEL descriptors for molecules/fragments.

Usage:
    python compute_descriptors.py
        -i {path to JSON with molecules, output of extract_fragments}
        -o {path to output csv file}
        -p {path to the PaDEL directory that contains PaDEL-Descriptor.jar}
"""

import os
import argparse
import logging
import json
import subprocess

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'


def create_parent_directory(path):
    """Create directory if it does not exists.

    :param path:
    :return:
    """
    dir_name = os.path.dirname(path)
    if not os.path.exists(dir_name) and not dir_name == "":
        os.makedirs(dir_name)


def read_configuration():
    """Get and return application settings.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Compute PaDEL descriptors for given'
                    'molecules/fragments.')
    parser.add_argument('-i', type=str, dest='input',
                        help='input JSON file',
                        required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='output CSV file', required=True)
    parser.add_argument('-p', type=str, dest='padel',
                        help='PaDEL directory', required=True)
    parser.add_argument('-f', dest='fragments',
                        help='use fragments instead of molecules',
                        action='store_true', required=False)

    return vars(parser.parse_args())


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    configuration = read_configuration()
    with open(configuration['input'], 'r') as stream:
        data = json.load(stream)
    create_parent_directory(configuration['output'])
    # Gather data
    smiles_set = set()
    if 'fragments' in configuration and configuration['fragments']:
        for molecule in data:
            for fragment in molecule['fragments']:
                if not fragment['smiles'] in smiles_set:
                    smiles_set.add(fragment['smiles'])
    else:
        for molecule in data:
            if not molecule['smiles'] in smiles_set:
                smiles_set.add(molecule['smiles'])
    # Prepare data for PaDEL.
    padel_input = os.path.dirname(configuration['output']) + '/PaDEL-temp.smi'
    with open(padel_input, 'w') as stream:
        for smiles in smiles_set:
            stream.write(smiles)
            stream.write('\t')
            stream.write(smiles)
            stream.write('\n')
    # Execute PaDEL
    logging.info('Executing PaDEL ...')
    thread = subprocess.Popen(
        ['java', '-jar',
         configuration['padel'] + '/PaDEL-Descriptor.jar',
         '-maxruntime', '5000',
         '-threads', '2',
         '-2d',
         '-dir', padel_input,
         '-file', configuration['output']],
        shell=True)
    thread.wait()
    logging.info('Executing PaDEL ... done')
    os.remove(padel_input)


if __name__ == '__main__':
    main()

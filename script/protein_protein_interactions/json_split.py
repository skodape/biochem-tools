#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Split JSON file with array.

The root entity must be an array. This script reads objects in the array
and split objects into a multiple JSON files.

The JSON file must not contains control characters ('{', '}') outside
the control positions (ie. in text).

"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals
import argparse
import os
import json
import logging

__author__ = 'Petr Škoda'
__license__ = 'X11'


def _read_json_array_stream(stream):
    """Read JSON objects from array.

    The text data must not contains control characters '{' or '}' in string.
    :param stream:
    :return:
    """
    # Read till the beginning of the array ie. first '['
    while True:
        char = stream.read(1)
        if char == '[':
            break
        if char == '':
            return
    #
    buffer = ''
    braces = 0
    while True:
        char = stream.read(1)
        if char == '':
            return
        # Skip separators between objects.
        if char == ',' and braces == 0:
            continue
        # Store characters.
        buffer += char
        if char == '{':
            braces += 1
        elif char == '}':
            braces -= 1
            if braces == 0:
                yield json.loads(buffer)
                buffer = ''


def write_file(path, value):
    logging.info('Writing file: %s', path)
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))

    with open(path, 'w') as output_stream:
        json.dump(value, output_stream)


def split_file(input_path, output_size, output_path_template):
    counter = 0
    with open(input_path, 'r') as input_stream:
        buffer = []
        for graph in _read_json_array_stream(input_stream):
            buffer.append(graph)
            if len(buffer) >= output_size:
                write_file(output_path_template.replace('{}', str(counter)),
                           buffer)
                buffer = []
                counter += 1
    # Write remaining data.
    if len(buffer) > 0:
        write_file(output_path_template.replace('{}', str(counter)),
                   buffer)


def _read_settings():
    """Get and return application settings.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='extract molecular fragments')
    parser.add_argument('-i', type=str, dest='input',
                        help='input JSON file', required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='template to JSON output file, the {} is '
                             'used as placeholder for index',
                        required=True)
    parser.add_argument('-s', type=int, dest='size',
                        help='number of object per output file', required=True)
    args = vars(parser.parse_args())
    #
    output = {
        'source': args['input'],
        'output': args['output'],
        'size': args['size']
    }
    return output


def _main():
    # Initialize logging.
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    logging.info('start')
    setting = _read_settings()
    split_file(setting['source'], setting['size'], setting['output'])


if __name__ == '__main__':
    _main()

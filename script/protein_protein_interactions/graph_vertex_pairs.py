#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute vertex pair fingerprints for given graphs.

"""

import argparse
import json
import logging
import numpy

__author__ = 'Petr Škoda'
__credits__ = ['Petr Škoda']
__license__ = 'X11'
__version__ = '1.0.0'
__maintainer__ = 'Petr Škoda'
__status__ = 'Development'


def read_json_array_stream(stream):
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


def read_configuration():
    """Get and return application settings.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='computes fingerprints for given graphs')
    parser.add_argument('-i', type=str, dest='input',
                        help='input SDF file', required=True)
    parser.add_argument('-c', type=str, dest='configuration',
                        help='configuration JSON file', required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='output SDF file', required=True)
    args = vars(parser.parse_args())
    #
    output = {
        'input': args['input'],
        'configuration': args['configuration'],
        'output': args['output']
    }
    return output


def warshall(vertices, edges):
    """Compute and return distance matrix for all vertices.

    :param vertices:
    :param edges:
    :return:
    """

    def get_distance(edges, x, y, no_path):
        """

        :param edges:
        :param x:
        :param y:
        :param no_path:
        :return:
        """
        if x == y:
            return 0
        for edge in edges:
            if (edge['from'] == x and edge['to'] == y) or \
                    (edge['from'] == y and edge['to'] == x):
                return 1
        return no_path

    n = len(vertices)
    distance_matrix = [[get_distance(edges, i, j, n + 1)
                        for j in vertices]
                       for i in vertices]

    for k in range(0, n):
        for i in range(0, n):
            for j in range(0, n):
                distance_matrix[i][j] = min(
                    distance_matrix[i][j],
                    distance_matrix[i][k] + distance_matrix[k][j])

    for k in range(0, n):
        for i in range(0, n):
            if distance_matrix[k][i] == n + 1:
                distance_matrix[k][i] = None

    return distance_matrix


def vertex_code(vertex, configuration):
    """Compute code for a single vertex.

    :param vertex:
    :param configuration:
    :return:
    """
    code = 0
    for item in configuration['vertex']['properties']:
        value = vertex[item['name']]
        if 'conversion' in item:
            if value in item['conversion']:
                value = item['conversion'][value]
            else:
                value = item['default']

        else:
            value = value % item['max']
        code = (code << int(item['size'])) + int(value)
    return code % configuration['vertex']['max']


def process_graph(graph, configuration):
    vertices = [item['id'] for item in graph['Vertices']]
    edges = graph['Edges']
    # Compute global graph descriptors.
    if configuration['distance']:
        distances = warshall(vertices, edges)
    #
    fingerprint = numpy.zeros(configuration['fp_size'])
    for left in range(0, len(vertices)):
        for right in range(0, len(vertices)):
            # Add properties.
            value = 0
            if configuration['distance']:
                value = (value << configuration['distance']['size']) + \
                        distances[left][right] % configuration['distance'][
                            'max']
            value = (value << configuration['vertex']['size']) + \
                    vertex_code(graph['Vertices'][left], configuration)
            value = (value << configuration['vertex']['size']) + \
                    vertex_code(graph['Vertices'][right], configuration)
            # Store into the fingerprint.
            fingerprint[value % configuration['fp_size']] = 1
    return fingerprint


def initialize_conversion_configuration(config):
    """Prepare conversion settings for use.

    In file we use size to determine the number of bits used, but we also
    need to know max size in order to module the value properly.
    :param config:
    :return:
    """
    if 'distance' in config:
        config['distance']['max'] = 1 << config['distance']['size']
    vertex_size = 0
    for item in config['vertex']['properties']:
        item['max'] = 1 << item['size']
        vertex_size += item['size']
    config['vertex']['size'] = vertex_size
    config['vertex']['max'] = 1 << vertex_size


def main():
    # Initialize logging.
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s.%(msecs)03d [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    logging.info('start')
    configuration = read_configuration()
    #
    with open(configuration['configuration'], 'r') as input_stream:
        conversion_configuration = json.load(input_stream)
    initialize_conversion_configuration(conversion_configuration)
    #
    counter = 0
    with open(configuration['output'], 'w') as output_stream:
        with open(configuration['input'], 'r') as input_stream:
            output_stream.write('[\n')
            first = True
            for graph in read_json_array_stream(input_stream):
                #
                id = graph['ID']
                value = process_graph(graph, conversion_configuration)
                # Write output.
                if first:
                    first = False
                    output_stream.write(' {\n')
                else:
                    output_stream.write(',{\n')
                output_stream.write('  "id":"' + str(id) + '",\n')
                output_stream.write('  "value": [')
                output_stream.write(str(int(value[0])))
                for x in value[1:]:
                    output_stream.write(',')
                    output_stream.write(str(int(x)))
                output_stream.write(']\n')
                output_stream.write(' }\n')
                # Progress
                counter += 1
                if counter % 1000 == 0:
                    logging.info(counter);
            output_stream.write(']')
    logging.info('done')


if __name__ == '__main__':
    main()

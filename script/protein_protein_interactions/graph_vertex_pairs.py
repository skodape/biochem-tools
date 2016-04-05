#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute vertex pair fingerprints for given graphs.

"""

import argparse
import json
import logging
import numpy
import math

__author__ = 'Petr Škoda'
__license__ = 'X11'


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

    def get_distance(x, y, no_path):
        """

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
    distance_matrix = [[get_distance(i, j, n + 1)
                        for j in vertices]
                       for i in vertices]

    for k in range(0, n):
        for i in range(0, n):
            for j in range(0, n):
                distance_matrix[i][j] = min(
                    distance_matrix[i][j],
                    distance_matrix[i][k] + distance_matrix[k][j])

    # Convert to indexes
    ver_list = [x for x in vertices]
    result_matrix = {}
    for k in range(0, n):
        result_matrix[ver_list[k]] = {}
        for i in range(0, n):
            if distance_matrix[k][i] == n + 1:
                result_matrix[ver_list[k]][ver_list[i]] = None
            else:
                result_matrix[ver_list[k]][ver_list[i]] = distance_matrix[k][i]
    # print('Warshall :', json.dumps(result_matrix, indent=1))
    return result_matrix


def process_template(item, template):
    value = item[template['property']]
    if template['type'] == 'property':
        if 'format' not in template:
            return value
        elif template['format'] == 'gray':
            # https://en.wikipedia.org/wiki/Gray_code
            return value ^ (value >> 1)
        else:
            return value
        # Use value as it is.
        pass
    elif template['type'] == 'mapping':
        try:
            value = template['map'][str(value)]
        except:
            logging.error('Item: %s', json.dumps(item, indent=2))
            logging.error('Template: %s', json.dumps(template, indent=2))
            raise Exception('Missing mapping for value: ' + str(value))
    elif template['type'] == 'binning':
        bin_found = False
        for bin_definition in template['bins']:
            if bin_definition['from'] <= value < bin_definition['to']:
                bin_found = True
                value = bin_definition['value']
                break
        if not bin_found:
            raise Exception('Missing bin for value: ' + str(value))
    else:
        value = 0
    return value


def get_vertex_code(index, vertices, configuration, molecule):
    vertex = vertices[index]
    # print('- - - - - - - - - - - -')
    result = 0
    shift = 0
    for item in configuration['fingerprint']['vertex']:
        value = process_template(vertex, item)
        # Store the value.
        if 'name' in item:
            # print('Store value: ', value, 'as', item['name'])
            vertex[item['name']] = value
        else:
            # print('Write value: ', value, ' as ',
            #       bin(int(value) % item['max']), ' ~ ',
            #       bin(int(value) % item['max'] << shift))
            # Put to output.
            result += (int(value) % item['max']) << shift
            shift += item['size']
            # print('    -> ', bin(result))
    # print('Vertex:', json.dumps(vertex, indent=2))
    # print('Vertex value:', bin(result % configuration['vertex_max']))
    return result % configuration['vertex_max']


def find_edge(left, right, edges):
    """Find and return edge or an empty object.

    :param left:
    :param right:
    :param edges:
    :return:
    """
    if right > left:
        swap = left
        left = right
        right = swap
    for item in edges:
        if item['from'] == left and item['to'] == right:
            return item
    return {}


def get_edge_code(left_index, right_index, vertices, edges, configuration,
                  molecule):
    left = vertices[left_index]
    right = vertices[right_index]
    edge = None
    #
    result = 0
    shift = 0
    # print('- - - - - - - - - - - -')
    for item in configuration['fingerprint']['edge']:
        if item['type'] == 'distance':
            # Artificial type of a property.
            value = molecule['distance'][left_index][right_index]
        elif item['type'] == 'compute':
            # Computed property.
            if item['method'] == 'euclidean_distance':
                value = 0
                for key in item['source']:
                    value += pow(float(left[key]) - float(right[key]), 2)
                value = math.sqrt(value)
            else:
                raise Exception('Unknown method: ' + item['method'])
        else:
            if edge is None:
                edge = find_edge(left_index, right_index, edges)
            value = process_template(edge, item)
        # Store the value.
        if 'name' in item:
            if edge is None:
                edge = find_edge(left_index, right_index, edges)
            edge[item['name']] = value
        else:
            # print('Write value:', value, 'as',
            #       bin(int(value) % item['max']), ' ~ ',
            #       bin(int(value) % item['max'] << shift))
            # Put to output.
            result += (int(value) % item['max']) << shift
            shift += item['size']
            # print('    -> ', bin(result))

    # print('Edge:', json.dumps(edge, indent=2))
    # print('Edge value:', bin(result % configuration['edge_max']))
    return result % configuration['edge_max']


def process_graph(graph, configuration):
    # print(json.dumps(graph, indent=2))

    vertices = {}
    for item in graph['Vertices']:
        vertices[item['id']] = item
    edges = graph['Edges']
    #
    fingerprint_size = configuration['fingerprint']['size']
    edge_size = configuration['fingerprint']['edge_size']
    vertex_size = configuration['fingerprint']['vertex_size']
    # Compute global descriptors.
    info = {
        'distance': warshall(vertices.keys(), edges)
    }
    # COMMENT OUT THIS AND USAGE TO IMPROVE PERFORMANCE
    indexes_count = 0
    indexes = set()
    indexes_hashed = set()
    #
    fingerprint = numpy.zeros(fingerprint_size)
    for left in vertices.keys():
        for right in vertices.keys():
            if left == right:
                continue
            # Add properties.
            left_code = get_vertex_code(left, vertices, configuration, info)
            right_code = get_vertex_code(right, vertices, configuration, info)
            edge_code = get_vertex_code(right, vertices, configuration, info)
            # Construct value.
            value = (left_code << (vertex_size + edge_size)) + \
                    (edge_code << vertex_size) + right_code
            # print('Index ', value, '(', value % fingerprint_size, ')',
            #       ' ~', bin(value))
            #
            indexes_count += 1
            indexes.add(value)
            indexes_hashed.add(value % fingerprint_size)
            # Store into the fingerprint.
            fingerprint[value % fingerprint_size] = 1
    # UNCOMMENT THIS TO SEE THE NUMBER OF UNIQUE AND HASHED FRAGMENTS
    # print('> size:', indexes_count, 'uniq:', len(indexes), 'hashed:',
    #       len(indexes_hashed), 'nodes:', len(vertices.keys()))
    # print('\t', indexes)
    # print('\t', indexes_hashed)
    return fingerprint


def initialize_conversion_configuration(configuration):
    vertex_size = 0
    for item in configuration['fingerprint']['vertex']:
        if 'size' not in item:
            continue
        vertex_size += item['size']
        item['max'] = 1 << item['size']
    configuration['fingerprint']['vertex_size'] = vertex_size
    configuration['vertex_max'] = 1 << vertex_size

    edge_size = 0
    for item in configuration['fingerprint']['edge']:
        if 'size' not in item:
            continue
        edge_size += item['size']
        item['max'] = 1 << item['size']
    configuration['fingerprint']['edge_size'] = edge_size
    configuration['edge_max'] = 1 << edge_size


def main():
    # Initialize logging.
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s.%(msecs)03d [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    logging.info('start')

    configuration = read_configuration()

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
                    logging.info(counter)
                    break
            #
            output_stream.write(']')
            print(output_stream.getvalue())
            output_stream.close()

    logging.info('done %d', counter)


if __name__ == '__main__':
    main()

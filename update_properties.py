#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Process a properties file output by Psi4 and add to properties.

Update the properties dictionary in the psi4_metadata to add their
properties and method from the properties.json file output during
a run of Psi4.
"""

import json
import argparse
from psi4_step import properties


def analyze(data, properties):
    method = data['_method']
    for key, value in data.items():
        if key[0] == '_':
            continue
        if key in properties:
            metadata = properties[key]
            methods = metadata['methods']
            if method not in methods:
                methods.append(method)
                methods.sort()
        else:
            properties[key] = {
                "calculation": [
                    "energy",
                    "optimization",
                    "thermodynamics",
                    "vibrations"
                ],
                "dimensionality": "scalar",
                'methods': [method],
                "description": "",
                "type": "float",
                "units": ""
            }
            if isinstance(value, list):
                properties[key]["dimensionality"] = [len(value)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'file',
        type=str,
        nargs='+',
        help="The 'properties.json' files to process"
    )

    args = parser.parse_args()

    for filename in args.file:
        with open(filename) as fd:
            data = json.load(fd)

        analyze(data, properties)

    print(json.dumps(properties, sort_keys=True, indent=4))

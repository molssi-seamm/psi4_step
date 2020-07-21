# -*- coding: utf-8 -*-

"""Setup and run Psi4"""

import json
import logging
from pathlib import Path

import psi4_step
import seamm
import seamm.data
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('psi4')


class Initialization(seamm.Node):

    def __init__(self, flowchart=None, title='Initialization', extension=None):
        """Initialize the node"""

        logger.debug('Creating Initialization {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.parameters = psi4_step.InitializationParameters()

        self.description = ['Initialization of Psi4']

    @property
    def header(self):
        """A printable header for this section of output"""
        return (
            'Step {}: {}'.format(
                '.'.join(str(e) for e in self._id), self.title
            )
        )

    @property
    def version(self):
        """The semantic version of this module.
        """
        return psi4_step.__version__

    @property
    def git_revision(self):
        """The git version of this module.
        """
        return psi4_step.__git_revision__

    def description_text(self, P=None):
        """Prepare information about what this node will do
        """
        if not P:
            P = self.parameters.values_to_dict()

        text = (
            f"Using the basis set '{P['basis']}' and a symmetry "
            f"tolerance of {P['symmetry_tolerance']}."
        )

        return self.header + '\n' + __(text, indent=4 * ' ').__str__()

    def get_input(self):
        """Get the input for an initialization calculation for Psi4"""
        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = '{:~P}'.format(PP[key])

        self.description = []
        self.description.append(
            __(self.description_text(PP), **PP, indent=self.indent)
        )

        # Start gathering the keywords
        result = []
        result.append('')
        result.append('#' * 80)
        result.append(f'# {self.header}')
        result.append('#' * 80)

        result.append(f"set basis {P['basis']}")
        result.append('')
        result.append(f"initial.symmetrize({P['symmetry_tolerance']})")
        result.append('point_group = initial.point_group().symbol()')

        # Dump the properties to a json file
        filename = f'@{self._id[-1]}+properties.json'
        result.append('')
        result.append("variables = {'point group': point_group}")
        result.append('')
        result.append(f"with open('{filename}', 'w') as fd:")
        result.append('    json.dump(variables, fd, sort_keys=True, indent=3)')

        return '\n'.join(result)

    def analyze(self, indent='', data={}, out=[]):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        # Read in the results from json
        directory = Path(self.directory)
        json_file = directory / 'properties.json'
        if json_file.exists():
            with json_file.open() as fd:
                data = json.load(fd)

            text = 'The calculation will use the point group {point group}.'
        else:
            data = {}
            text = (
                '\nThere are no results from Psi4. Perhaps it '
                f'failed? Looking for {str(json_file)}.'
            ),

        printer.normal(__(text, **data, indent=self.indent + 4 * ' '))

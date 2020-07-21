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


class Energy(seamm.Node):

    def __init__(self, flowchart=None, title='Energy', extension=None):
        """Initialize the node"""

        logger.debug('Creating Energy {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.parameters = psi4_step.EnergyParameters()

        self.description = 'A single point energy calculation'

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

        if P['level'] == 'recommended':
            method = P['method']
        else:
            method = P['advanced_method']

        if psi4_step.methods[method]['method'] == 'dft':
            if P['level'] == 'recommended':
                functional = P['functional']
            else:
                functional = P['advanced_functional']
            text = f'Single-point energy using {method} with an '
            text += f'exchange-correlation potential of {functional}'
            if (
                len(psi4_step.dft_functionals[functional]['dispersion']) > 1
                and P['dispersion'] != 'none'
            ):
                text += f" with the {P['dispersion']} dispersion correction."
            else:
                text += " with no dispersion correction."
        else:
            text = f'Single-point energy using {method}.'

        # Spin
        if P['spin-restricted'] == 'yes':
            text += " The spin will be restricted to a pure eigenstate."
        elif P['spin-restricted']:
            text += " The spin will not be restricted and may not be a "
            text += "proper eigenstate."
        else:
            text += " Whether the spin will be restricted to a pure "
            text += "eigenstate will be determined by {P['spin-restricted']}"

        return self.header + '\n' + __(text, **P, indent=4 * ' ').__str__()

    def get_input(self):
        """Get the input for an energy calculation for Psi4"""
        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        # references = self.parent.references

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

        lines = []
        lines.append('')
        lines.append('#' * 80)
        lines.append(f'# {self.header}')
        lines.append('#' * 80)

        # Figure out what we are doing!
        if P['level'] == 'recommended':
            method_string = P['method']
        else:
            method_string = P['advanced_method']

        method = psi4_step.methods[method_string]['method']

        lines.append('set scf_type df')
        lines.append('set guess sad')

        spin_restricted = P['spin-restricted'] == 'yes'
        if method == 'dft':
            if spin_restricted:
                lines.append('set reference rks')
            else:
                lines.append('set reference uks')
        else:
            if spin_restricted:
                lines.append('set reference rhf')
            else:
                lines.append('set reference uhf')
        lines.append('')
        if method == 'dft':
            if P['level'] == 'recommended':
                functional_string = P['functional']
            else:
                functional_string = P['advanced_functional']
            functional = psi4_step.dft_functionals[functional_string]['name']
            if len(
                psi4_step.dft_functionals[functional_string]['dispersion']
            ) > 1:
                functional = functional + '-' + P['dispersion']
            lines.append(f"E, wfn = energy('{functional}', return_wfn=True)")
        else:
            lines.append(f"E, wfn = energy('{method}', return_wfn=True)")

        # Dump the properties to a json file
        filename = f'@{self._id[-1]}+properties.json'
        lines.append('')
        lines.append('oeprop(')
        lines.append('    wfn,')
        lines.append("    'MULTIPOLE(5)',")
        lines.append("    'ESP_AT_NUCLEI',")
        lines.append("    'MO_EXTENTS',")
        lines.append("    'LOWDIN_CHARGES',")
        lines.append("    'MULLIKEN_CHARGES',")
        lines.append("    'WIBERG_LOWDIN_INDICES',")
        lines.append("    'MAYER_INDICES',")
        lines.append("    'NO_OCCUPATIONS',")
        lines.append("    title='PROP'")
        lines.append(')')
        lines.append('')
        lines.append('variables = wfn.scalar_variables()')
        lines.append('arrays = wfn.array_variables()')
        lines.append('for item in arrays:')
        lines.append(
            '    variables[item] = wfn.array_variable(item).np.tolist()'
        )
        lines.append("variables['E'] = E")
        lines.append('')
        lines.append(f"with open('{filename}', 'w') as fd:")
        lines.append('    json.dump(variables, fd, sort_keys=True, indent=3)')

        return '\n'.join(lines)

    def analyze(self, indent='', data={}, out=[]):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """

        # P = self.parameters.current_values_to_dict(
        #     context=seamm.flowchart_variables._data
        # )

        # Read in the results from json
        directory = Path(self.directory)
        json_file = directory / 'properties.json'
        if json_file.exists():
            with json_file.open() as fd:
                data = json.load(fd)

            # Put any requested results into variables or tables
            self.store_results(
                data=data,
                properties=psi4_step.properties,
                results=self.parameters['results'].value,
                create_tables=self.parameters['create tables'].get()
            )

            text = 'The calculated energy is {E:.6f} Ha.'
        else:
            data = {}
            text = (
                '\nThere are no results from Psi4. Perhaps it '
                f'failed? Looking for {str(json_file)}.'
            ),

        printer.normal(__(text, **data, indent=self.indent + 4 * ' '))

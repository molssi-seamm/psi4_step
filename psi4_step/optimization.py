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


class Optimization(psi4_step.Energy):

    def __init__(self, flowchart=None, title='Optimization', extension=None):
        """Initialize the node"""

        logger.debug('Creating Optimization {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.parameters = psi4_step.OptimizationParameters()

        self.description = 'A geometry optimization'

    def description_text(self, P=None):
        """Prepare information about what this node will do
        """

        if not P:
            P = self.parameters.values_to_dict()

        text = super().description_text(
            P=P, calculation_type='Geometry optimization'
        )

        added = (
            'The geometry optimization will use the {optimization method} '
        )
        if P['max geometry steps'] == 'default':
            added += 'method, using the default maximum number of steps, which'
            added += ' is based on the system size.'
        else:
            added += 'method, with no more than {max geometry steps} steps.'

        if P['geometry convergence'] == 'Custom':
            added += ' The convergence criterion is'
        else:
            added += " The convergence criterion is '{geometry convergence}'."

        if P['recalc hessian'] != 'never':
            added += " The Hessian will be recalculated every {recalc hession}"
            added += " steps. Note that calculating the second derivatives is "
            added += "quite expensive!"

        return text + '\n' + __(added, **P, indent=4 * ' ').__str__()

    def get_input(self, calculation_type='optimize'):
        """Get the input for an optimization calculation for Psi4"""
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
        lines.append('')
        lines.append('set opt_type min')
        opt_method = psi4_step.optimization_methods[P['optimization method']]
        lines.append(f'set step_type {opt_method}')
        max_steps = 50
        lines.append(f'set geom_maxiter {max_steps}')
        if P['geometry convergence'] == 'Custom':
            pass
        else:
            convergence = psi4_step.optimization_convergence[
                P['geometry convergence']]
            lines.append(f'set g_convergence {convergence}')
        lines.append('')

        # Add in the input from the energy part of things
        lines.append(super().get_input(calculation_type=calculation_type))

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

            text = 'The calculated energy is {Eelec:.6f} Ha.'
        else:
            data = {}
            text = (
                '\nThere are no results from Psi4. Perhaps it '
                f'failed? Looking for {str(json_file)}.'
            ),

        printer.normal(__(text, **data, indent=self.indent + 4 * ' '))

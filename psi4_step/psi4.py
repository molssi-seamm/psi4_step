# -*- coding: utf-8 -*-

"""Non-graphical part of the Psi4 step in a SEAMM flowchart
"""

import configargparse
import cpuinfo
import json
import logging
from pathlib import Path
import pprint

import psi4_step
import seamm
from seamm import data  # noqa: F401
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: 'job' and 'printer'. 'job' send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# 'printer' sends output to the file 'step.out' in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('Psi4')


class Psi4(seamm.Node):
    """
    The non-graphical part of a Psi4 step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : Psi4Parameters
        The control parameters for Psi4.

    See Also
    --------
    TkPsi4,
    Psi4, Psi4Parameters
    """

    def __init__(
        self,
        flowchart=None,
        title='Psi4',
        namespace='org.molssi.seamm.psi4',
        extension=None
    ):
        """A step for Psi4 in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
            flowchart: seamm.Flowchart
                The non-graphical flowchart that contains this step.

            title: str
                The name displayed in the flowchart.
            namespace: The namespace for the plugins of the subflowchart
            extension: None
                Not yet implemented
        """
        logger.debug('Creating Psi4 {}'.format(self))

        # Argument/config parsing
        self.parser = configargparse.ArgParser(
            auto_env_var_prefix='',
            default_config_files=[
                '/etcs/seamm/psi4.in',
                '/etc/seamm/psi4_step.ini',
                '/etc/seamm/seamm.ini',
                '~/.seamm/psi4.ini',
                '~/.seamm/psi4_step.ini',
                '~/.seamm/seamm.ini',
            ]
        )

        self.parser.add_argument(
            '--seamm-configfile',
            is_config_file=True,
            default=None,
            help='a configuration file to override others'
        )

        # Options for this plugin
        self.parser.add_argument(
            "--psi4-step-log-level",
            default=configargparse.SUPPRESS,
            choices=[
                'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
            ],
            type=str.upper,
            help="the logging level for the Psi4 step"
        )

        # Options for Psi4
        self.parser.add_argument(
            '--psi4-exe',
            default='psi4',
            help='the path to the Psi4 executable'
        )

        self.parser.add_argument(
            '--psi4-num-threads',
            default='default',
            help='How many threads to use in Psi4'
        )

        self.options, self.unknown = self.parser.parse_known_args()

        # Set the logging level for this module if requested
        if 'psi4_step_log_level' in self.options:
            logger.setLevel(self.options.psi4_step_log_level)
        self.subflowchart = seamm.Flowchart(
            parent=self, name='Psi4', namespace=namespace
        )

        super().__init__(
            flowchart=flowchart, title='Psi4', extension=extension
        )

        self.parameters = psi4_step.Psi4Parameters()

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

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        self._id = node_id

        # and set our subnodes
        self.subflowchart.set_ids(self._id)

        return self.next()

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
            P: dict
                An optional dictionary of the current values of the control
                parameters.
        Returns
        -------
            description : str
                A description of the current step.
        """

        self.subflowchart.root_directory = self.flowchart.root_directory

        # Get the first real node
        node = self.subflowchart.get_node('1').next()

        text = self.header + '\n\n'
        while node is not None:
            text += __(node.description_text(), indent=3 * ' ').__str__()
            text += '\n'
            node = node.next()

        return text

    def run(self):
        """Run a Psi4 step.

        Returns
        -------

        next_node : seamm.Node
            The next node object in the flowchart.

        """
        printer.important(self.header)
        printer.important('')

        if data.structure is None:
            logger.error('Psi4 run(): there is no structure!')
            raise RuntimeError('Psi4 run(): there is no structure!')

        # Access the options
        o = self.options

        # How many processors does this node have?
        info = cpuinfo.get_cpu_info()
        n_cores = info['count']
        # Account for Intel hyperthreading
        if info['arch'][0:3] == 'X86':
            n_cores = int(n_cores / 2)
        if n_cores < 1:
            n_cores = 1
        logger.info('The number of cores is {}'.format(n_cores))

        if o.psi4_num_threads == 'default':
            psi4_num_threads = n_cores
        else:
            psi4_num_threads = int(o.psi4_num_threads)
        if psi4_num_threads > n_cores:
            psi4_num_threads = n_cores
        if psi4_num_threads < 1:
            psi4_num_threads = 1
        logger.info('Psi4 will use {} threads.'.format(psi4_num_threads))

        # Work through the subflowchart to find out what to do.
        self.subflowchart.root_directory = self.flowchart.root_directory

        next_node = super().run(printer)

        # Get the first real node
        node = self.subflowchart.get_node('1').next()

        input_data = []
        input_data.append('import json')
        input_data.append('import numpy as np')
        input_data.append('')

        # Put the structure into the input
        input_data.append(self._convert_structure(name='initial'))

        while node is not None:
            text = node.get_input()
            input_data.append(text)
            node = node.next()

        # Write out the final structure
        input_data.append('')
        input_data.append('# Write the final structure to disk')
        input_data.append('molecule = get_active_molecule()')
        input_data.append('tmp = molecule.to_dict()')
        input_data.append('for item, value in tmp.items():')
        input_data.append('    if isinstance(value, np.ndarray):')
        input_data.append('        tmp[item] = value.tolist()')
        input_data.append('')
        input_data.append("with open('final_structure.json', 'w') as fd:")
        input_data.append('    json.dump(tmp, fd, sort_keys=True, indent=3)')

        files = {'input.dat': '\n'.join(input_data)}
        logger.info('input.dat:\n' + files['input.dat'])

        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)
        for filename in files:
            path = directory / filename
            with path.open(mode='w') as fd:
                fd.write(files[filename])

        return_files = ['output.dat', '*properties.json', '*structure.json']
        env = {'PSIPATH': Path(o.psi4_exe).parent}

        local = seamm.ExecLocal()
        result = local.run(
            cmd=[o.psi4_exe, f'-n {psi4_num_threads}'],
            files=files,
            return_files=return_files,
            env=env
        )  # yapf: disable

        if result is None:
            logger.error('There was an error running Psi4')
            return None

        logger.debug('\n' + pprint.pformat(result))

        logger.debug('stdout:\n' + result['stdout'])
        if 'stdout' in result and result['stdout'] != '':
            path = directory / 'stdout.txt'
            with path.open(mode='w') as fd:
                fd.write(result['stdout'])

        if result['stderr'] != '':
            logger.warning('stderr:\n' + result['stderr'])
            path = directory / 'stderr.txt'
            with path.open(mode='w') as fd:
                fd.write(result['stderr'])

        for filename in result['files']:
            if filename[0] == '@':
                subdir, fname = filename[1:].split('+')
                path = directory / subdir / fname
            else:
                path = directory / filename
            with path.open(mode='w') as fd:
                if result[filename]['data'] is not None:
                    fd.write(result[filename]['data'])
                else:
                    fd.write(result[filename]['exception'])

        # Analyze the results
        self.analyze()

        return next_node

    def analyze(self, indent='', **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        'printer'.

        Parameters
        ----------
            indent: str
                An extra indentation for the output
        """

        # Get the first real node
        node = self.subflowchart.get_node('1').next()

        # Loop over the subnodes, asking them to do their analysis
        while node is not None:
            for value in node.description:
                printer.important(value)

            node.analyze()

            printer.normal('')

            node = node.next()

        # Update the structure
        directory = Path(self.directory)
        structure_file = directory / 'final_structure.json'
        if structure_file.exists():
            with structure_file.open(mode='r') as fd:
                structure = json.load(fd)
        if 'geom' in structure:
            system = seamm.data.structure
            atoms = system['atoms']
            xyz = []
            it = iter(structure['geom'])
            for x in it:
                xyz.append([x, next(it), next(it)])
            atoms['coordinates'] = xyz
            printer.important(
                self.indent +
                '    Updated the system with the structure from Psi4',
            )
            printer.important('')

    def _convert_structure(self, name=None):
        """Convert the structure to the input for Psi4."""

        system = seamm.data.structure

        structure = []
        if name is None:
            structure.append('molecule {')
        else:
            structure.append('molecule ' + name + ' {')

        # Charge and multiplicity
        if 'extras' in system:
            extras = system['extras']

            if 'open' in extras:
                openshell = extras['open']
                if (
                    (
                        isinstance(openshell, tuple) or
                        isinstance(openshell, list)
                    ) and len(openshell) > 1
                ):
                    nopen = openshell[0]
                    norbitals = openshell[1]
                    if nopen != norbitals:
                        raise NotImplementedError(
                            f"Handling of open shell = '{openshell}'"
                        )
                else:
                    nopen = openshell
                    norbitals = nopen

                if 'net_charge' in extras:
                    structure.append(f"    {extras['net_charge']}   {nopen}")
                else:
                    structure.append(f'    0   {nopen}')
            else:
                if 'net_charge' in extras:
                    structure.append(f"    {extras['net_charge']}   0")

        elements = system['atoms']['elements']
        coordinates = system['atoms']['coordinates']

        if 'freeze' in system['atoms']:
            freeze = system['atoms']['freeze']
        else:
            freeze = [''] * len(elements)

        for element, xyz, frz in zip(elements, coordinates, freeze):
            x, y, z = xyz
            structure.append(
                f'    {element:2s} {float(x): 12.8f} {float(y): 12.8f} '
                f'{float(z): 12.8f}'
            )
        structure.append('}')

        return '\n'.join(structure) + '\n'

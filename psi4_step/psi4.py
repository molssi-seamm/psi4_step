# -*- coding: utf-8 -*-
"""Non-graphical part of the Psi4 step in a SEAMM flowchart
"""

import configargparse
import logging
import seamm
from seamm import data  # noqa: F401
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
import psi4_step

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


def upcase(string):
    """Return an uppercase version of the string.

    Used for the type argument in argparse
    """
    return string.upper()


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
    def __init__(self,
                 flowchart=None,
                 title='Psi4',
                 namespace='org.molssi.seamm.psi4_step',
                 extension=None):
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
                '/etc/seamm/psi4_step.ini',
                '/etc/seamm/seamm.ini',
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
            type=upcase,
            help="the logging level for the Psi4 step"
        )

        self.options, self.unknown = self.parser.parse_known_args()

        # Set the logging level for this module if requested
        if 'psi4_step_log_level' in self.options:
            logger.setLevel(self.options.psi4_step_log_level)
        self.sub_flowchart = seamm.Flowchart(
            parent=self, name='Psi4',
            namespace=namespace)

        super().__init__(
            flowchart=flowchart,
            title='Psi4',
            extension=extension)

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

        if not P:
            P = self.parameters.values_to_dict()

        text = ('Please replace this with a short summary of the '
                'Psi4 step, including key parameters.')

        return self.header + '\n' + __(text, **P, indent=4 * ' ').__str__()

    def run(self):
        """Run a Psi4 step.

        Returns
        -------

        next_node : seamm.Node
            The next node object in the flowchart.

        """

        next_node = super().run(printer)
        # Get the first real node
        node = self.sub_flowchart.get_node('1').next()

        input_data = []
        while node is not None:
            keywords = node.get_input()
            input_data.append(' '.join(keywords))
            node = node.next()

        files = {'molssi.dat': '\n'.join(input_data)}
        logger.info('molssi.dat:\n' + files['molssi.dat'])

        local = seamm.ExecLocal()
        result = local.run(
            cmd=['Psi4', '-in', 'molssi.dat'],
            files=files,
            return_files=[])

        if result is None:
            logger.error('There was an error running Psi4')
            return None

        logger.debug('\n' + pprint.pformat(result))

        logger.info('stdout:\n' + result['stdout'])
        if result['stderr'] != '':
            logger.warning('stderr:\n' + result['stderr'])

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
        node = self.sub_flowchart.get_node('1').next()

        # Loop over the subnodes, asking them to do their analysis
        while node is not None:
            for value in node.description:
                printer.important(value)
                printer.important(' ')

            node.analyze()

            node = node.next()

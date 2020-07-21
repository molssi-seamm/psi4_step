# -*- coding: utf-8 -*-

"""The graphical part of a Psi4 Initialization node"""

import configargparse
import logging
import tkinter as tk

import seamm
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkInitialization(seamm.TkNode):

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=120,
        y=20,
        w=200,
        h=50,
        my_logger=logger,
        keyword_metadata=None
    ):
        """Initialize the graphical Tk Psi4 initialization step

        Keyword arguments:
        """
        self.results_widgets = []

        # Argument/config parsing
        self.parser = configargparse.ArgParser(
            auto_env_var_prefix='',
            default_config_files=[
                '/etc/seamm/psi4_initialization.ini',
                '/etc/seamm/seamm.ini',
                '~/.seamm/psi4_initialization.ini',
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
            "--psi4-tk-initialization-log-level",
            default=configargparse.SUPPRESS,
            choices=[
                'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
            ],
            type=lambda string: string.upper(),
            help="the logging level for the Psi4 Tk_initialization step"
        )

        self.options, self.unknown = self.parser.parse_known_args()

        # Set the logging level for this module if requested
        if 'psi4_tk_initialization_log_level' in self.options:
            my_logger.setLevel(self.options.psi4_tk_initialization_log_level)
            my_logger.critical(
                'Set log level to {}'.format(
                    self.options.psi4_tk_initialization_log_level
                )
            )

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
            my_logger=my_logger,
            keyword_metadata=keyword_metadata
        )

    def right_click(self, event):
        """Probably need to add our dialog...
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self, title='Edit Psi4 Initialization Step'):
        """Create the dialog!"""
        self.logger.debug('Creating the dialog')
        frame = super().create_dialog(title=title)

        # Create all the widgets
        P = self.node.parameters
        for key in P:
            if key not in ('results', 'extra keywords', 'create tables'):
                self[key] = P[key].widget(frame)

        self.logger.debug('Finished creating the dialog')

    def reset_dialog(self, widget=None):
        frame = self['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        row = 0
        self['basis'].grid(row=row, column=0, sticky=tk.EW)
        widgets.append(self['basis'])
        row += 1
        self['symmetry_tolerance'].grid(row=row, column=0, sticky=tk.EW)
        widgets.append(self['symmetry_tolerance'])
        row += 1
        self['calculate gradients'].grid(row=row, column=0, sticky=tk.EW)
        widgets.append(self['calculate gradients'])
        row += 1
        sw.align_labels(widgets)

        return row

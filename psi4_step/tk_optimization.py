# -*- coding: utf-8 -*-

"""The graphical part of a Psi4 Optimization node"""

import configargparse
import logging
import pprint
import tkinter as tk
import tkinter.ttk as ttk

import psi4_step
# import seamm
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkOptimization(psi4_step.TkEnergy):

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
        """Initialize the graphical Tk Psi4 optimization step

        Keyword arguments:
        """
        self.results_widgets = []

        # Argument/config parsing
        self.parser = configargparse.ArgParser(
            auto_env_var_prefix='',
            default_config_files=[
                '/etc/seamm/psi4_optimization.ini',
                '/etc/seamm/seamm.ini',
                '~/.seamm/psi4_optimization.ini',
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
            "--psi4-tk-optimization-log-level",
            default=configargparse.SUPPRESS,
            choices=[
                'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
            ],
            type=lambda string: string.upper(),
            help="the logging level for the Psi4 Tk_optimization step"
        )

        self.options, self.unknown = self.parser.parse_known_args()

        # Set the logging level for this module if requested
        if 'psi4_tk_optimization_log_level' in self.options:
            my_logger.setLevel(self.options.psi4_tk_optimization_log_level)
            my_logger.critical(
                'Set log level to {}'.format(
                    self.options.psi4_tk_optimization_log_level
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

    def create_dialog(
        self, title='Edit Psi4 Optimization Step', calculation='optimization'
    ):
        """Create the edit dialog!

        This is reasonably complicated, so a bit of description
        is in order. The superclass Energy creates the dialog
        along with the calculation parameters in a 'calculation'
        frame..

        This method adds a second frame for controlling the optimizer.

        The layout is handled in part by the Energy superclass, which
        handles the calculation frame. Our part is handled by two
        methods:

        * reset_dialog does the general layout of the main frames.
        * reset_optimization handles the layout of the optimization
          section.
        """

        logger.debug('TkOptimization.create_dialog')

        # Let parent classes do their thing.
        super().create_dialog(title=title, calculation=calculation)

        # Shortcut for parameters
        P = self.node.parameters

        logger.debug('Parameters:\n{}'.format(pprint.pformat(P.to_dict())))

        # Frame to isolate widgets
        opt_frame = self['optimization'] = ttk.LabelFrame(
            self['frame'],
            borderwidth=4,
            relief='sunken',
            text='Geometry Optimization',
            labelanchor='n',
            padding=10
        )

        for key in psi4_step.OptimizationParameters.parameters:
            self[key] = P[key].widget(opt_frame)

        # and binding to change as needed
        self['geometry convergence'].combobox.bind(
            "<<ComboboxSelected>>", self.reset_optimization
        )
        self['geometry convergence'].combobox.bind(
            "<Return>", self.reset_optimization
        )
        self['geometry convergence'].combobox.bind(
            "<FocusOut>", self.reset_optimization
        )

        # Top level needs to call reset_dialog
        if calculation == 'optimization':
            self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets, letting our parents go first."""
        row = super().reset_dialog()

        self['optimization'].grid(row=row, column=0)
        row += 1

        self.reset_optimization()

        return row

    def reset_optimization(self, widget=None):
        convergence = self['geometry convergence'].get()

        frame = self['optimization']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        # widgets2 = []
        row = 0

        self['optimization method'].grid(
            row=row, column=0, columnspan=2, sticky=tk.EW
        )
        widgets.append(self['optimization method'])
        row += 1

        self['max geometry steps'].grid(
            row=row, column=0, columnspan=2, sticky=tk.EW
        )
        widgets.append(self['max geometry steps'])
        row += 1

        self['geometry convergence'].grid(
            row=row, column=0, columnspan=2, sticky=tk.EW
        )
        widgets.append(self['geometry convergence'])
        row += 1

        if convergence == 'Custom':
            pass

        self['recalc hessian'].grid(
            row=row, column=0, columnspan=2, sticky=tk.EW
        )
        widgets.append(self['recalc hessian'])
        row += 1

        sw.align_labels(widgets)
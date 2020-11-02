# -*- coding: utf-8 -*-

"""The graphical part of a Psi4 Energy node"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import psi4_step
import seamm
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkEnergy(seamm.TkNode):

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
        """Initialize the graphical Tk Psi4 energy step

        Keyword arguments:
        """
        self.results_widgets = []
        self._calculation = None

        # Set the logging level for this module if requested
        # if 'psi4_tk_energy_log_level' in self.options:
        #     my_logger.setLevel(self.options.psi4_tk_energy_log_level)
        #     my_logger.critical(
        #         'Set log level to {}'.format(
        #             self.options.psi4_tk_energy_log_level
        #         )
        #     )

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
        self, title='Edit Psi4 Energy Step', calculation='energy'
    ):
        """Create the dialog!"""
        self.logger.debug('Creating the dialog')
        self._calculation = calculation
        frame = super().create_dialog(
            title=title, widget='notebook', results_tab=True
        )

        # Create a frame for the calculation control
        self['calculation'] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief='sunken',
            text='Calculation',
            labelanchor='n',
            padding=10
        )

        # Create all the widgets
        P = self.node.parameters
        for key in psi4_step.EnergyParameters.parameters:
            if key not in ('results', 'extra keywords', 'create tables'):
                self[key] = P[key].widget(self['calculation'])

        # bindings...
        self['level'].combobox.bind(
            "<<ComboboxSelected>>", self.reset_calculation
        )
        self['level'].config(state='readonly')

        self['method'].combobox.bind(
            "<<ComboboxSelected>>", self.reset_calculation
        )
        self['method'].combobox.bind("<Return>", self.reset_calculation)
        self['method'].combobox.bind("<FocusOut>", self.reset_calculation)

        self['advanced_method'].combobox.bind(
            "<<ComboboxSelected>>", self.reset_calculation
        )
        self['advanced_method'].combobox.bind(
            "<Return>", self.reset_calculation
        )
        self['advanced_method'].combobox.bind(
            "<FocusOut>", self.reset_calculation
        )

        # self.setup_results(psi4_step.properties, calculation=calculation)

        # Top level needs to call reset_dialog
        if calculation == 'energy':
            self.reset_dialog()

        self.logger.debug('Finished creating the dialog')

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        frame = self['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        self['calculation'].grid(row=0, column=0)
        self.reset_calculation()
        return 1

    def reset_calculation(self, widget=None):
        level = self['level'].get()

        if level == 'recommended':
            long_method = self['method'].get()
            if self.is_expr(long_method):
                method = None
                meta = None
            else:
                method = psi4_step.methods[long_method]['method']
                meta = psi4_step.methods[long_method]
            functional = self['functional'].get()
        else:
            long_method = self['advanced_method'].get()
            if self.is_expr(long_method):
                method = None
                meta = None
            else:
                method = psi4_step.methods[long_method]['method']
                meta = psi4_step.methods[long_method]
            functional = self['advanced_functional'].get()

        # Set up the results table because it depends on the method
        self.results_widgets = []
        self.setup_results(
            psi4_step.properties, calculation=self._calculation, method=method
        )

        frame = self['calculation']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        widgets2 = []
        row = 0
        self['level'].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        row += 1
        if level == 'recommended':
            self['method'].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self['method'])
            row += 1
            if method is None or method == 'dft':
                self['functional'].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self['functional'])
                row += 1
            if meta is None or 'freeze core?' in meta and meta['freeze core?']:
                self['freeze-cores'].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self['freeze-cores'])
                row += 1
        else:
            self['advanced_method'].grid(
                row=row, column=0, columnspan=2, sticky=tk.EW
            )
            widgets.append(self['advanced_method'])
            row += 1
            if method is None or method == 'dft':
                self['advanced_functional'].grid(
                    row=row, column=1, sticky=tk.EW
                )
                widgets2.append(self['advanced_functional'])
                row += 1
            if meta is None or 'freeze core?' in meta and meta['freeze core?']:
                self['freeze-cores'].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self['freeze-cores'])
                row += 1
        if method is None or method == 'dft':
            dispersions = psi4_step.dft_functionals[functional]['dispersion']
            if len(dispersions) > 1:
                w = self['dispersion']
                w.config(values=dispersions)
                if w.get() not in dispersions:
                    w.value(dispersions[1])
                w.grid(row=row, column=1, sticky=tk.W)
                widgets2.append(self['dispersion'])
                row += 1
                sw.align_labels(widgets2)
            frame.columnconfigure(0, minsize=30)
        self['spin-restricted'].grid(
            row=row, column=0, columnspan=2, sticky=tk.EW
        )
        widgets.append(self['spin-restricted'])
        row += 1
        sw.align_labels(widgets)

        return row
